#!/bin/sh

# run program on Dell.

#BSUB -oo log
#BSUB -eo log
#BSUB -q debug
#BSUB -P GFS-DEV
#BSUB -J vegt
#BSUB -W 0:03
#BSUB -x                 # run not shared
#BSUB -n 1               # total tasks
#BSUB -R span[ptile=1]   # tasks per node
#BSUB -R affinity[core(1):distribute=balance]

set -x

module purge
module load ips/18.0.1.163
module load HDF5-serial/1.10.1
module load NetCDF/4.5.0

rm -f vegt.nc

 ./vegt.exe
