#build program on wcoss2

set -x

module reset

module load PrgEnv-intel/8.1.0
module load craype/2.7.13
module load intel/19.1.3.304
module load cray-mpich/8.1.7


module load hdf5/1.10.6
module load netcdf/4.7.4

module list

make clean
make
