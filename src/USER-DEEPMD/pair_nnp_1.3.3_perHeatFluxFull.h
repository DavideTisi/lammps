/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(deepmd,PairNNP)

#else

#ifndef LMP_PAIR_NNP_H
#define LMP_PAIR_NNP_H

#include "pair.h"
#include "NNPInter.h"
#include <iostream>
#include <fstream>

#define GIT_SUMM 
#define GIT_HASH abee8b9
#define GIT_BRANCH master
#define GIT_DATE 2021-05-12 23:55:50 +0800
#ifdef HIGH_PREC
#define FLOAT_PREC double
#else
#define FLOAT_PREC float
#endif
#define DEEPMD_ROOT /home/dtisi/deepmd-1.3.3/deepmdlib
#define TensorFlow_INCLUDE_DIRS /opt/sissa/libs/libtensorflow_cc/1.14.0/include;/opt/sissa/libs/libtensorflow_cc/1.14.0/include
#define TensorFlow_LIBRARY /opt/sissa/libs/libtensorflow_cc/1.14.0/lib/libtensorflow_cc.so;/opt/sissa/libs/libtensorflow_cc/1.14.0/lib/libtensorflow_framework.so
#define DPMD_CVT_STR(x) #x
#define DPMD_CVT_ASSTR(X) DPMD_CVT_STR(X)
#define STR_GIT_SUMM DPMD_CVT_ASSTR(GIT_SUMM)
#define STR_GIT_HASH DPMD_CVT_ASSTR(GIT_HASH)
#define STR_GIT_BRANCH DPMD_CVT_ASSTR(GIT_BRANCH)
#define STR_GIT_DATE DPMD_CVT_ASSTR(GIT_DATE)
#define STR_FLOAT_PREC DPMD_CVT_ASSTR(FLOAT_PREC)
#define STR_DEEPMD_ROOT DPMD_CVT_ASSTR(DEEPMD_ROOT)
#define STR_TensorFlow_INCLUDE_DIRS DPMD_CVT_ASSTR(TensorFlow_INCLUDE_DIRS)
#define STR_TensorFlow_LIBRARY DPMD_CVT_ASSTR(TensorFlow_LIBRARY)

namespace LAMMPS_NS {

class PairNNP : public Pair {
 public:
  PairNNP(class LAMMPS *);
  virtual ~PairNNP();
  virtual void compute(int, int);
  virtual void *extract(const char *, int &);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int i, int j);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void print_summary(const string pre) const;
  double **vatomc ;// complete vatom matrix 3*3 not simmetric
  int get_node_rank();
 protected:  
  virtual void allocate();
  double **scale;
  virtual void ev_setup_vatomc( int, int alloc = 1);
  int maxvatomc ; 

private:  
  NNPInter nnp_inter;
  NNPInterModelDevi nnp_inter_model_devi;
  unsigned numb_models;
  double cutoff;
  int numb_types;
  vector<vector<double > > all_force;
  ofstream fp;
  int out_freq;
  string out_file;
  int dim_fparam;
  int dim_aparam;
  int out_each;
  int out_rel;
  bool single_model;
  bool multi_models_mod_devi;
  bool multi_models_no_mod_devi;
#ifdef HIGH_PREC
  vector<double > fparam;
  vector<double > aparam;
  double eps;
#else
  vector<float > fparam;
  vector<float > aparam;
  float eps;
#endif
  void make_ttm_aparam(
#ifdef HIGH_PREC
      vector<double > & dparam
#else
      vector<float > & dparam
#endif
      );
  bool do_ttm;
  string ttm_fix_id;
};

}

#endif
#endif
