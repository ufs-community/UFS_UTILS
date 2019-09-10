#!/bin/sh

set -eux

gfortran -o filter_topo -I${NETCDF}/include -g -O2 -fcray-pointer -ffree-line-length-none -fno-range-check -fbacktrace -fdefault-real-8 -fdefault-double-8 filter_topo.F90 -L${NETCDF}/lib -L${HDF5}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lz
