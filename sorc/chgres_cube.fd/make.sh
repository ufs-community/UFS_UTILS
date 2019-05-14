#!/bin/sh

set -x

if [[ ! -d ../exec ]];then
  mkdir ../exec
fi

mac=$(hostname -f)

case $mac in

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON WCOSS CRAY.
#---------------------------------------------------------------------------------

llogin? | slogin?)

  module purge
  module load modules/3.2.6.7
  module load PrgEnv-intel/5.2.56
  module rm intel
  module load intel/16.3.210
  module load cray-mpich/7.2.0
  module load craype-haswell
  module load cray-netcdf
  module load w3nco-intel/2.0.6
  module load nemsio-intel/2.2.3
  module load bacio-intel/2.0.2
  module load sp-intel/2.0.2
  module load sigio-intel/2.1.0
  module load sfcio-intel/1.0.0

# module use /gpfs/hps3/emc/nems/noscrub/emc.nemspara/soft/modulefiles
# module load esmf/7.1.0r
  export ESMFMKFILE=/gpfs/hps3/emc/global/noscrub/George.Gayno/esmf/8_0_0_bs20/lib/esmf.mk

  export FCOMP=ftn
  export FFLAGS="-O0 -g -r8 -i4 -qopenmp -convert big_endian -check bounds -assume byterecl -warn unused"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON THEIA.
#---------------------------------------------------------------------------------

tfe??)

  #source /apps/lmod/lmod/init/sh
  module purge
  module load intel/18.1.163
  module load impi/5.1.1.109
  module load netcdf/4.3.0
  module load hdf5/1.8.14
  module load wgrib2/2.0.8

  module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
  module load w3nco/v2.0.6
  module load w3emc/v2.3.0
  module load nemsio/v2.2.3
  module load bacio/v2.0.1
  module load sp/v2.0.2
  module load sfcio/v1.0.0
  module load sigio/v2.0.1

  module use /scratch4/NCEPDEV/nems/noscrub/emc.nemspara/soft/modulefiles
  module load yaml-cpp
  module load esmf/8.0.0bs01

  export FCOMP=mpiifort
  export FFLAGS="-O0 -g -traceback -r8 -i4 -qopenmp -convert big_endian -check bounds -warn unused -assume byterecl"

  export NEMSIOGFS_INC=/scratch4/NCEPDEV/global/save/glopara/svn/nceplibs/branches/nemsiogfsv2.0.1/include/nemsiogfs
  export NEMSIOGFS_LIB=/scratch4/NCEPDEV/global/save/glopara/svn/nceplibs/branches/nemsiogfsv2.0.1/libnemsiogfs.a

  export WGRIB2API_LIB="/apps/wgrib2/2.0.8/intel/18.0.3.222/lib/libwgrib2_api.a"
  export WGRIB2API_INC="/apps/wgrib2/2.0.8/intel/18.0.3.222/lib"
  export WGRIB2_LIB="/apps/wgrib2/2.0.8/intel/18.0.3.222/lib/libwgrib2.a"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON DELL.
#---------------------------------------------------------------------------------

m????.ncep.noaa.gov | v????.ncep.noaa.gov )

  module purge
  module use /usrx/local/dev/modulefiles

  module load ips/18.0.1.163
  module load impi/18.0.1
  module load NetCDF/4.5.0
# module load ESMF/7_1_0r
  export ESMFMKFILE=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/esmf_lib/8_0_0bs20/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk
  module load w3nco/2.0.6
  module load sp/2.0.2
  module load nemsio/2.2.3
  module load bacio/2.0.2
  module load sfcio/1.0.0
  module load sigio/2.1.0

  export FCOMP=mpif90
  export FFLAGS="-O0 -g -traceback -r8 -i4 -qopenmp -convert big_endian -check bounds -warn unused -assume byterecl"

  make clean
  make
  rc=$? ;;

#---------------------------------------------------------------------------------
# BUILD PROGRAM ON NSSL ODIN
#---------------------------------------------------------------------------------

nid?????)

  module use /home/larissa.reames/modulefiles

  #module switch PrgEnv-cray/6.0.4 PrgEnv-intel

  module load cray-netcdf-hdf5parallel
  module load cray-parallel-netcdf
  module load cray-hdf5-parallel
  #cray-mpich/7.6.2
  module load sp/v2.0.2
  module load nemsio/v2.2.2
  module load nemsiogfs/v2.0.1
  module load w3nco/v2.0.6
  module load bacio/v2.0.2
  module load esmf/8.0.0.bs1
  module load sigio/v2.0.1
  module load sfcio/v1.0.0
  
  export FCOMP=ftn
  export WGRIB2API_LIB="/home/larissa.reames/tmp/wgrib2-2/grib2/lib/libwgrib2_api.a"
  export WGRIB2API_INC="/home/larissa.reames/tmp/wgrib2-2/grib2/lib"
  export WGRIB2_LIB="/home/larissa.reames/tmp/wgrib2-2/grib2/lib/libwgrib2.a"

  export FFLAGS="-O0 -g -traceback -r8 -i4 -qopenmp -convert big_endian -check bounds -assume byterecl -warn unused"

  make clean
  make
  rc=$?  ;;

esac

exit
