#!/bin/sh

set -x

case $1 in
  "cray" )
     export FCMP=ftn
     export FFLAGS="-O2 -g"
     export LIBS=" "  ;;
  "wcoss" )
     export FCMP=ifort
     export FFLAGS="-O2 -g ${NETCDF_FFLAGS}"
     export LIBS="${NETCDF_LDFLAGS_F}" ;;
  "wcoss_dell_p3" )
     export FCMP=ifort
     export FFLAGS="-O2 -g ${NETCDF_FFLAGS}"
     export LIBS="${NETCDF_LDFLAGS_F} ${HDF5_LDFLAGS_F}" ;;
  "jet" )
     export FCMP=ifort
     export FFLAGS="-O2 -g -I$NETCDF/include" 
     export LIBS="-L$NETCDF/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5 -lhdf5_fortran" ;;
  "theia" )
     export FCMP=ifort
     export FFLAGS="-O2 -g -I$NETCDF/include"
     export LIBS="-L$NETCDF/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5 -lhdf5_fortran" ;;
  *)
     echo "REGIONAL GRID UTILITY BUILD NOT TESTED ON MACHINE $1"
     exit 1 ;;
esac

make clean
make
make install

if ((rc != 0)); then
  echo "ERROR BUILDING REGIONAL GRID UTILITY"
  exit $rc
else
  exit 0
fi
