#%Module#####################################################
## Build and run module for Odin
#############################################################

module load craype/2.6.2
module load craype-ivybridge
module load PrgEnv-intel
module swap intel/19.0.5.281
module load cray-mpich/7.7.10
module load cray-libsci
module load cray-netcdf-hdf5parallel
module load cray-parallel-netcdf
module load cray-hdf5-parallel

export NETCDF=/opt/cray/pe/netcdf-hdf5parallel/4.6.3.2/INTEL/19.0

#module use -a /oldscratch/ywang/external/modulefiles
module use /oldscratch/ywang/external/NCEPLIBS_SRW/modules
module load w3nco
module load w3emc
module load sp
module load ip
module load bacio
module load sigio
module load sfcio
module load nemsio
module load nemsiogfs
module load gfsio
module load landsfcutil
module load g2
module load wgrib2

#module load esmf/8.0.0
export ESMFMKFILE=/oldscratch/ywang/external/NCEPLIBS_SRW/lib64/esmf.mk

export CMAKE_Fortran_COMPILER=ftn
export CMAKE_C_COMPILER=cc

#export WGRIB2_ROOT=/oldscratch/ywang/external/NCEPLIBS_SRW/wgrib2-2.0.8
