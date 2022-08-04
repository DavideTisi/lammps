// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_bin_atomonly.h"

#include "atom.h"
#include "error.h"
#include "neighbor.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;
using namespace NeighConst;

/* ---------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI>
NPairBinAtomonly<HALF, NEWTON, TRI>::NPairBinAtomonly(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI>
void NPairBinAtomonly<HALF, NEWTON, TRI>::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,bin_start;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    ibin = atom2bin[i];

    for (k = 0; k < nstencil; k++) {
      bin_start = binhead[ibin+stencil[k]];
      if (stencil[k] == 0) {
        if (HALF && NEWTON && (!TRI)) {
          // Half neighbor list, newton on, orthonormal
          // loop over rest of atoms in i's bin, ghosts are at end of linked list
          bin_start = bins[i];
        }
      }

      for (j = bin_start; j >= 0; j = bins[j]) {
        if (!HALF) {
          // Full neighbor list
          // only skip i = j
          if (i == j) continue;
        } else if (!NEWTON) {
          // Half neighbor list, newton off
          // only store pair if i < j
          // stores own/own pairs only once
          // stores own/ghost pairs on both procs
          if (j <= i) continue;
        } else if (TRI) {
          // Half neighbor list, newton on, triclinic
          // pairs for atoms j "below" i are excluded
          // below = lower z or (equal z and lower y) or (equal zy and lower x)
          //         (equal zyx and j <= i)
          // latter excludes self-self interaction but allows superposed atoms
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp) {
              if (x[j][0] < xtmp) continue;
              if (x[j][0] == xtmp && j <= i) continue;
            }
          }
        } else {
          // Half neighbor list, newton on, orthonormal
          // store every pair for every bin in stencil,except for i's bin

          if (stencil[k] == 0) {
            // if j is owned atom, store it, since j is beyond i in linked list
            // if j is ghost, only store if j coords are "above and to the "right" of i
            if (j >= nlocal) {
              if (x[j][2] < ztmp) continue;
              if (x[j][2] == ztmp) {
                if (x[j][1] < ytmp) continue;
                if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
              }
            }
          }
        }

        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  if (!HALF) list->gnum = 0;
}

namespace LAMMPS_NS {
template class NPairBinAtomonly<0,1,0>;
template class NPairBinAtomonly<1,1,0>;
}
