#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle regression test on Hera.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in regression.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J cycle_reg_test
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o regression.log
#SBATCH -e regression.log
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH -q debug
#SBATCH -t 00:05:00

set -x

. /apps/lmod/lmod/init/sh
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load szip/2.1
module load hdf5/1.10.5
module use /scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles
module load netcdf_parallel/4.7.4
module list

export DATA=/scratch2/NCEPDEV/stmp1/$LOGNAME/reg_tests.cycle

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export HOMEreg=/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2

export APRUNCY="srun"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/apps/nccmp/1.8.5/intel/18.0.3.051/bin/nccmp

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
