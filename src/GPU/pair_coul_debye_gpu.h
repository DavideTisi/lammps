/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(coul/debye/gpu,PairCoulDebyeGPU);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_DEBYE_GPU_H
#define LMP_PAIR_COUL_DEBYE_GPU_H

#include "pair_coul_debye.h"

namespace LAMMPS_NS {

class PairCoulDebyeGPU : public PairCoulDebye {
 public:
  PairCoulDebyeGPU(LAMMPS *lmp);
  ~PairCoulDebyeGPU() override;
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int) override;
  void init_style() override;
  void reinit() override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
};

}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Pair style coul/debye/gpu requires atom attribute q

The atom style defined does not have this attribute.

E: Cannot use newton pair with coul/debye/gpu pair style

Self-explanatory.

*/
