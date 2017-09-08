
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2017-09-08 09:59:24.619000";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/global/u1/d/dkhatami/src/ecp-astro/chimera2castro";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux cori06 4.4.59-92.24-default #1 SMP Thu Jun 22 14:29:09 UTC 2017 (d11a83a) x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "/global/homes/d/dkhatami/src/ecp-astro/amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "6.3.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "CC";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "ftn";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = " -g -O3  -fopenmp -DNDEBUG -DBL_USE_MPI -DBL_USE_OMP -DDO_PROBLEM_POST_INIT -DAMREX_GIT_VERSION=\"17.09\" -DBL_GCC_VERSION='6.3.0' -DBL_GCC_MAJOR_VERSION=6 -DBL_GCC_MINOR_VERSION=3 -DBL_SPACEDIM=2 -DAMREX_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -DBL_Linux -DCRSEGRNDOMP -DSPONGE -DGRAVITY -DSELF_GRAVITY -DPOINTMASS -DMG_USE_FBOXLIB -DBL_USE_F_BASELIB -DBL_USE_FORTRAN_MPI -I. -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/LinearSolvers/C_to_F_MG -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/LinearSolvers/F_MG -I/opt/cray/pe/hdf5-parallel/1.10.0.3/GNU/5.1/include -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/Base -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/AmrCore -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/Amr -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/Boundary -I/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/util -I. -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/driver -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/driver/param_includes -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/hydro -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/problems -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/sources -I/global/homes/d/dkhatami/src/ecp-astro/Castro/constants -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Util/model_parser -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Source/gravity -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Microphysics/EOS -I/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/EOS/helmholtz -I/global/homes/d/dkhatami/src/ecp-astro/Castro/Microphysics/networks -I/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/networks/aprox13 -I/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/EOS -I/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/networks -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Src/F_BaseLib -I/global/homes/d/dkhatami/src/ecp-astro/amrex/Tools/C_scripts";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = " -g -O3 -ffree-line-length-none -fno-range-check -fno-second-underscore -Jtmp_build_dir/o/2d.gnu.MPI.OMP.EXE -I tmp_build_dir/o/2d.gnu.MPI.OMP.EXE -fimplicit-none  -fopenmp";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "-L/opt/cray/pe/hdf5-parallel/1.10.0.3/GNU/5.1/lib -L. -L/opt/gcc/6.3.0/snos/lib/gcc/x86_64-suse-linux/6.3.0/../../../../lib64/ -L/opt/cray/pe/hdf5-parallel/1.10.0.3/GNU/5.1/lib";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "-lmpichf90 -lhdf5 -lhdf5_fortran -lhdf5 -lz -lgfortran -lquadmath";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 2;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";
  static const char AUX1[] = "EOS";
  static const char AUX2[] = "NETWORK";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;
    case 2: return AUX2;

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";
  static const char AUX1[] = "/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/EOS/helmholtz";
  static const char AUX2[] = "/global/homes/d/dkhatami/src/ecp-astro/Microphysics-starkiller/networks/aprox13";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;
    case 2: return AUX2;

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "";
  static const char HASH2[] = "17.09";
  static const char HASH3[] = "17.09";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;
    case 3: return HASH3;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

}
