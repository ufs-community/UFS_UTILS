#!/bin/sh

# build program on Orion.

set -x

source /apps/lmod/lmod/init/sh
module purge

module load intel/2018.4
module load impi/2018.4
module load netcdf/4.7.2
module load hdf5/1.10.5

make clean
make
