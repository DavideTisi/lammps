/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/debye/dielectric,PairLJCutCoulDebyeDielectric);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_DEBYE_DIELECTRIC_H
#define LMP_PAIR_LJ_CUT_COUL_DEBYE_DIELECTRIC_H

#include "pair_lj_cut_coul_debye.h"

namespace LAMMPS_NS {

class PairLJCutCoulDebyeDielectric : public PairLJCutCoulDebye {
 public:
  PairLJCutCoulDebyeDielectric(class LAMMPS *);
  ~PairLJCutCoulDebyeDielectric() override;
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void init_style() override;

  double **efield;
  double *epot;

 protected:
  class AtomVecDielectric *avec;
  int nmax;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/cut/coul/debye/dielectric requires atom attribute q

The atom style defined does not have this attribute.

*/
