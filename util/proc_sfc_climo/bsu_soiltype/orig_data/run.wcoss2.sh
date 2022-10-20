#!/bin/bash

# Run script for WCOSS2

#PBS -l walltime=00:03:00
#PBS -o log
#PBS -e log
#PBS -N soil
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l select=1:ncpus=1:mem=10GB

module load PrgEnv-intel/8.1.0
module load craype/2.7.13
module load intel/19.1.3.304
module load cray-mpich/8.1.7

module load hdf5/1.10.6
module load netcdf/4.7.4

module list

WORK_DIR=/lfs/h2/emc/stmp/george.gayno/bnu.soil
rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

$PBS_O_WORKDIR/soil.exe

date
nccopy -d 1 -m 100000000 soil.nc soil.compressed.nc
date

exit
