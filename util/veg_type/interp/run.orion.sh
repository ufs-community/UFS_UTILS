#!/bin/sh --login

# run program on Orion.

#SBATCH -J vegt
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

set -x

module load intel/2018.4
module load impi/2018.4
module load netcdf/4.7.2
module load hdf5/1.10.5

rm -f test.nc

srun ./vegt.exe
