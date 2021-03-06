// **************************************************************************
//                                   lj_smooth.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for acceleration of the lj/smooth pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
#else
_texture_2d( pos_tex,int4);
#endif
#else
#define pos_tex x_
#endif

__kernel void k_lj_smooth(const __global numtyp4 *restrict x_,
                   const __global numtyp4 *restrict lj1,
                   const __global numtyp4 *restrict lj3,
                   const __global numtyp4 *restrict ljsw,
                   const int lj_types,
                   const __global numtyp *restrict sp_lj,
                   const __global int * dev_nbor,
                   const __global int * dev_packed,
                   __global acctyp4 *restrict ans,
                   __global acctyp *restrict engv,
                   const int eflag, const int vflag, const int inum,
                   const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;

    numtyp factor_lj;
    numtyp r, t, tsq, fskin;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
      
      int mtype=itype*lj_types+jtype;
      if (rsq<lj1[mtype].z) {
        numtyp force, r6inv;
        numtyp r2inv=ucl_recip(rsq);
        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          force = factor_lj*r2inv*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        }
        else {
          r = sqrt(rsq);
          t = r - lj3[mtype].w;
          tsq = t*t;
          fskin = ljsw[mtype].x + ljsw[mtype].y*t +
            ljsw[mtype].z*tsq + ljsw[mtype].w*tsq*t;
          force = factor_lj*r2inv*fskin*r;
        }

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e;
          if (rsq < lj1[mtype].w)
            e = r6inv * (lj3[mtype].x*r6inv - lj3[mtype].y) - lj3[mtype].z;
          else
            e = ljsw[mtype].x - ljsw[mtype].x*t -
              ljsw[mtype].y*tsq/2.0 - ljsw[mtype].z*tsq*t/3.0 -
              ljsw[mtype].z*tsq*tsq/4.0 - lj3[mtype].z;

          //numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_lj_smooth_fast(const __global numtyp4 *restrict x_,
                        const __global numtyp4 *restrict lj1_in,
                        const __global numtyp4 *restrict lj3_in,
                        const __global numtyp4 *restrict ljsw,
                        const __global numtyp *restrict sp_lj_in,
                        const __global int * dev_nbor,
                        const __global int * dev_packed,
                        __global acctyp4 *restrict ans,
                        __global acctyp *restrict engv,
                        const int eflag, const int vflag, const int inum,
                        const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (eflag>0)
      lj3[tid]=lj3_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    numtyp r, t, tsq, fskin;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<lj1[mtype].z) {
        numtyp force, r6inv;
        numtyp r2inv=ucl_recip(rsq);
        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          force = factor_lj*r2inv*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        }
        else {
          r = sqrt(rsq);
          t = r - lj3[mtype].w;
          //printf("%f\n", r - lj3[mtype].w);
          tsq = t*t;
          fskin = ljsw[mtype].x + ljsw[mtype].y*t +
            ljsw[mtype].z*tsq + ljsw[mtype].w*tsq*t;
          force = factor_lj*r2inv*fskin*r;
        }

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e;
          if (rsq < lj1[mtype].w)
            e = r6inv * (lj3[mtype].x*r6inv - lj3[mtype].y) - lj3[mtype].z;
          else
            e = ljsw[mtype].x - ljsw[mtype].x*t - 
              ljsw[mtype].y*tsq/2.0 - ljsw[mtype].z*tsq*t/3.0 -
              ljsw[mtype].z*tsq*tsq/4.0 - lj3[mtype].z; //???

          //numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

