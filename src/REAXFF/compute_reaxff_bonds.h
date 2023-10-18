/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (LANL)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(reaxff/bonds,ComputeReaxFFBonds);
// clang-format on
#else

#ifndef LMP_COMPUTE_REAXFF_BONDS_H
#define LMP_COMPUTE_REAXFF_BONDS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReaxFFBonds : public Compute {
 public:
  ComputeReaxFFBonds(class LAMMPS *, int, char **);
  ~ComputeReaxFFBonds() override;
  void init() override;
  void compute_local() override;
  void compute_peratom() override;
  virtual void compute_bonds();
  double memory_usage() override;

 protected:
  bigint invoked_bonds;     // last timestep on which compute_bonds() was invoked
  int nlocal;
  int nbonds;
  int prev_nbonds;

  tagint **neighid;
  double **abo;
  int *bondcount;
  class PairReaxFF *reaxff;

 private:
  int FindBond();
};

}    // namespace LAMMPS_NS

#endif
#endif
