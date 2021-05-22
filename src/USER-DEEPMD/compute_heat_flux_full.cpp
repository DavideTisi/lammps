/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: German Samolyuk (ORNL) and
                         Mario Pinto (Computational Research Lab, Pune, India)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "compute_heat_flux_full.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_nnp.h"
#include "group.h"
#include "error.h"
#include "comm.h"
//#include "comm_brick.h"

using namespace LAMMPS_NS;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */
//void ComputeHeatFluxFull::print_matrix(std::ofstream ofsforce, double ** f,std::string filename,int rows,int cols)
void ComputeHeatFluxFull::print_matrix(std::ofstream& ofsforce, double ** f,int rows,int cols)
{
//int nall=io->atom->nlocal+io->atom->nghost ; 
// std::ofstream ofsforce;
//  ofsforce.open(filename ,std::ofstream::out | std::ofstream::app);
 ofsforce<<" nuova TIMESTEP "<<std::endl;
 ofsforce<<std::endl;
for (int ii = 0; ii < rows; ++ii){
  ofsforce<<ii<<" ";
  for (int dd = 0; dd < cols; ++dd){
    ofsforce<<f[ii][dd]<<" ";
  } 
  ofsforce<<std::endl;
}
ofsforce<<std::endl;
//ofsforce.close();
 
}


/* ---------------------------------------------------------------------- */
//void ComputeHeatFluxFull::comunicate(CommBrick * comm)
//{
//  int n ;
//  MPI_Request request;
//  AtomVec *avec = atom->avec; 
//  double **f = force->pair->vatom;
//  double *buf;
//  int nswap = comm-> nswap ;
//  int *sendproc = comm->sendproc;  
//  int *recvproc = comm->recvproc;  
//  int *size_reverse_recv = comm->size_reverse_recv;  
//  int *size_reverse_send = comm->size_reverse_send;  
//  
//  for (int iswap = nswap-1; iswap >= 0; iswap--) {
//    if (sendproc[iswap] != comm->me) {
//        if (size_reverse_recv[iswap])
//          MPI_Irecv(buf_recv,size_reverse_recv[iswap]*3,MPI_DOUBLE,
//                    sendproc[iswap],0,world,&request);
//        if (size_reverse_send[iswap]) {
//          buf = f[firstrecv[iswap]];
//          MPI_Send(buf,size_reverse_send[iswap]*3,MPI_DOUBLE,
//                   recvproc[iswap],0,world);
//        }
//        if (size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
//    }
//  }
//
//
//}
/* ---------------------------------------------------------------------- */

ComputeHeatFluxFull::ComputeHeatFluxFull(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_ke(NULL), id_pe(NULL)
{
  if (narg != 5 && narg!= 6) error->all(FLERR,"Illegal compute heat/flux/full command");
//  peratom_flag = 1;
  timeflag = 1;
  pressatomflag = 1;
  vector_flag = 1;
  if (narg==6) {
     error->message(FLERR,"Adding some debug tests!");
     size_vector = 9+9;
     vector = new double[9+9];
     debug=true;
  } else {
     size_vector = 9;
     vector = new double[9];
     debug=false;
  }
  extvector = 1;

  // store ke/atom, pe/atom, stress/atom IDs used by heat flux computation
  // insure they are valid for these computations

  int n = strlen(arg[3]) + 1;
  id_ke = new char[n];
  strcpy(id_ke,arg[3]);

  n = strlen(arg[4]) + 1;
  id_pe = new char[n];
  strcpy(id_pe,arg[4]);


  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  if (ike < 0 || ipe < 0 )
    error->all(FLERR,"Could not find compute heat/flux compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");

  if (!force->pair || strcmp(force->pair_style,"deepmd") != 0)
    error->all(FLERR,"Compute heat/flux/full needs deempd pair style");
 
}

/* ---------------------------------------------------------------------- */

ComputeHeatFluxFull::~ComputeHeatFluxFull()
{
  delete [] id_ke;
  delete [] id_pe;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxFull::init()
{
  // error checks

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  if (ike < 0 || ipe < 0 )
    error->all(FLERR,"Could not find compute heat/flux compute ID");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxFull::compute_vector()
{
  invoked_vector = update->ntimestep;

  // invoke 3 computes if they haven't been already

  if (!(c_ke->invoked_flag & INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= INVOKED_PERATOM;
  }

  // heat flux vector = jc[3] + jv[3]
  // jc[3] = convective portion of heat flux = sum_i (ke_i + pe_i) v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included

  double *ke = c_ke->vector_atom;
  double *pe = c_pe->vector_atom;
  double **stress = force->pair->vatom;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double jc[3] = {0.0,0.0,0.0};
  double jv[3] = {0.0,0.0,0.0};
  double vt[9] = {0.0,0.0,0.0,
                  0.0,0.0,0.0,
                  0.0,0.0,0.0};
  double eng;

  //std::ofstream ofvatom;
  //std::ofstream ofvatom_dopo;
  //std::ofstream ofvel;
  //ofvatom.open("vatom_lammps.dat",std::ofstream::out | std::ofstream::app);
  //ofvatom_dopo.open("vatom_lammps_dopo.dat",std::ofstream::out | std::ofstream::app);
  //ofvel.open("vel_lammps.dat",std::ofstream::out | std::ofstream::app);
  //print_matrix(ofvatom,stress,nall,9);
  //print_matrix(ofvel,v,nall,3);
  ////comm->reverse_comm_stress();
  ////comunicate(comm);
  //print_matrix(ofvatom_dopo,stress,nall,9);
 // ofvatom<<endl;
 // ofvatom<<" nuova TIMESTEP "<<endl;
 // ofvatom<<endl;
 // ofvel<<endl;
 // ofvel<<" nuova TIMESTEP "<<endl;
 // ofvel<<endl;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      eng = pe[i] + ke[i];
      jc[0] += eng*v[i][0];
      jc[1] += eng*v[i][1];
      jc[2] += eng*v[i][2];
      jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] +
        stress[i][4]*v[i][2];
      jv[1] -= stress[i][6]*v[i][0] + stress[i][1]*v[i][1] +
        stress[i][5]*v[i][2];
      jv[2] -= stress[i][7]*v[i][0] + stress[i][8]*v[i][1] +
        stress[i][2]*v[i][2];
      if (debug) {
        for (unsigned int i0=0;i0<3;i0++) {
          for (unsigned int i1=0;i1<3;i1++) {
            vt[i0*3+i1]+=stress[i][i0*3+i1];
          }
        }
      }
    }
    //j+=1;
    //if(j>=nlocal) j=0;
 //   ofvatom<<i<<" ";
 //   ofvel<<i<<" ";
 //   ofvatom<<stress[i][0]<<" "<<stress[i][1]<<" "<<stress[i][2]<<" "<<stress[i][3]<<" "<<stress[i][4]<<" "<<stress[i][5]<<" "<<stress[i][6]<<" "<<stress[i][7]<<" "<<stress[i][8]<<" "<<endl;
 //   ofvel<<v[i][0]<<" "<<v[i][1]<<" "<<v[i][2]<<endl;
  }
 // ofvatom.close();
 // ofvel.close();
  // convert jv from stress*volume to energy units via nktv2p factor

  double nktv2p = force->nktv2p;
 // jv[0] /= nktv2p;
 // jv[1] /= nktv2p;
 // jv[2] /= nktv2p;

  // sum across all procs
  // 1st 3 terms are total heat flux
  // 2nd 3 terms are just conductive portion
  if (debug){
  cout<< std::scientific<<std::setprecision(20)<<"nktv2p "<< nktv2p<<endl;
    double data[9+9] = {jc[0]+jv[0], jc[1]+jv[1],jc[2]+jv[2],
                        jc[0],       jc[1],      jc[2],      
                        jv[0],       jv[1],      jv[2],
   //                     vt[0]*nktv2p, vt[8]*nktv2p, vt[7]*nktv2p,
   //                     vt[3]*nktv2p, vt[1]*nktv2p, vt[6]*nktv2p,
   //                     vt[4]*nktv2p, vt[5]*nktv2p, vt[2]*nktv2p
                        vt[0], vt[6], vt[7],
                        vt[3], vt[1], vt[8],
                        vt[4], vt[5], vt[2]
                         };     
    MPI_Allreduce(data,vector,9+9,MPI_DOUBLE,MPI_SUM,world); 

  } else {
    double data[9] = {jc[0]+jv[0], jc[1]+jv[1],jc[2]+jv[2],
                      jc[0],       jc[1],      jc[2],
                      jv[0],       jv[1],      jv[2],
                      };
    MPI_Allreduce(data,vector,9,MPI_DOUBLE,MPI_SUM,world);
  }
}
