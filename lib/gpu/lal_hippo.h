/***************************************************************************
                                  hippo.h
                             -------------------
                          Trung Dac Nguyen (Northwestern)

  Class for acceleration of the hippo pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#ifndef LAL_HIPPO_H
#define LAL_HIPPO_H

#include "lal_base_amoeba.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class Hippo : public BaseAmoeba<numtyp, acctyp> {
 public:
  Hippo();
  ~Hippo();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    *
    * Returns:
    * -  0 if successful
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, const int max_amtype, const int max_amclass,
           const double *host_pdamp, const double *host_thole,
           const double *host_dirdamp, const int *host_amtype2class,
           const double *host_special_mpole,
           const double *host_special_hal,
           const double *host_special_repel,
           const double *host_special_disp,
           const double *host_special_polar_wscale,
           const double *host_special_polar_piscale,
           const double *host_special_polar_pscale,
           const double *host_csix, const double *host_adisp,
           const double *host_pcore, const double *host_palpha,
           const int nlocal, const int nall, const int max_nbors,
           const int maxspecial, const int maxspecial15, const double cell_size,
           const double gpu_split, FILE *_screen,
           const double polar_dscale, const double polar_uscale);

  /// Reallocate per-atom arrays if needed, and build neighbor lists once, if needed
  int** precompute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, int *host_amtype,
                int *host_amgroup, double **host_rpole, double **host_uind,
                double **host_uinp, double* host_pval, double *sublo, double *subhi,
                tagint *tag, int **nspecial, tagint **special,
                int *nspecial15, tagint **special15,
                const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **&ilist, int **&numj, const double cpu_time, bool &success,
                double *charge, double *boxlo, double *prd);

  /// Compute dispersion real-space with device neighboring
  int** compute_dispersion_real(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, int *host_amtype,
                int *host_amgroup, double **host_rpole, double *sublo, double *subhi,
                tagint *tag, int **nspecial, tagint **special,
                int *nspecial15, tagint **special15,
                const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **numj, const double cpu_time, bool &success,
                const double aewald, const double off2_disp, double *charge,
                double *boxlo, double *prd);

  /// Compute multipole real-space with device neighboring
  virtual int** compute_multipole_real(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, int *host_amtype,
                int *host_amgroup, double **host_rpole, double *host_pval,
                double *sublo, double *subhi, tagint *tag, int **nspecial, tagint **special,
                int *nspecial15, tagint **special15,
                const bool eflag, const bool vflag,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **numj, const double cpu_time, bool &success,
                const double aewald, const double felec, const double off2_mpole, double *charge,
                double *boxlo, double *prd, void **tep_ptr);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// pdamp = coeff_amtype.x; thole = coeff_amtype.y;
  /// dirdamp = coeff_amtype.z; amtype2class = coeff_amtype.w
  UCL_D_Vec<numtyp4> coeff_amtype;
  /// csix = coeff_amclass.x; adisp = coeff_amclass.y;
  UCL_D_Vec<numtyp4> coeff_amclass;
  /// Special polar values [0-4]: 
  ///   sp_polar.x = special_polar_wscale
  ///   sp_polar.y special_polar_pscale,
  ///   sp_polar.z = special_polar_piscale
  ///   sp_polar.w = special_mpole
  UCL_D_Vec<numtyp4> sp_polar;
  /// Special nonpolar values [0-4]: 
  ///   sp_nonpolar.x = special_hal
  ///   sp_nonpolar.y special_repel
  ///   sp_nonpolar.z = special_disp
  UCL_D_Vec<numtyp4> sp_nonpolar;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _polar_dscale, _polar_uscale;
  numtyp _qqrd2e;

  UCL_Kernel k_dispersion;

 protected:
  bool _allocated;
  int dispersion_real(const int eflag, const int vflag);
  int multipole_real(const int eflag, const int vflag);
  int udirect2b(const int eflag, const int vflag);
  int umutual2b(const int eflag, const int vflag);
  int polar_real(const int eflag, const int vflag);
  
};

}

#endif
