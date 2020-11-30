#!/bin/sh

# build program on Dell

set -x

module purge

module load ips/18.0.1.163

module load HDF5-serial/1.10.1
module load NetCDF/4.5.0

make clean
make
