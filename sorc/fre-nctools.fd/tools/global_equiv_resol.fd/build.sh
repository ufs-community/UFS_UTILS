#!/bin/sh

set -x

case $1 in
  "cray" )
     export FC=ftn
     export FFLAGS="-O2 -g"
     export LIBS=" "  ;;
  "wcoss" )
     export FC=ifort
     export FFLAGS="-O2 -g ${NETCDF_FFLAGS}"
     export LIBS="${NETCDF_LDFLAGS_F}" ;;
  "wcoss_dell_p3" )
     export FC=ifort
     export FFLAGS="-O2 -g ${NETCDF_FFLAGS}"
     export LIBS="${NETCDF_LDFLAGS_F} ${HDF5_LDFLAGS_F}" ;;
  "jet" )
     export FC=ifort
     export FFLAGS="-O2 -g -I$NETCDF/include" 
     export LIBS="-L$NETCDF/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5 -lhdf5_fortran" ;;
  "hera" )
     export FC=ifort
     export FFLAGS="-O2 -g -I$NETCDF/include"
     export LIBS="-L$NETCDF/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5 -lhdf5_fortran" ;;
  *)
     echo "GLOBAL_EQUIV_RESOL UTILITY BUILD NOT TESTED ON MACHINE $1"
     exit 1 ;;
esac

make clean
make
make install

if ((rc != 0)); then
  echo "ERROR BUILDING GLOBAL_EQUIV_RESOL UTILITY"
  exit $rc
else
  exit 0
fi
