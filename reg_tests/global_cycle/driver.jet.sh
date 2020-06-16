#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle regression test on Jet.
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
#SBATCH -A emcda
#SBATCH --open-mode=truncate
#SBATCH -o regression.log
#SBATCH -e regression.log
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --partition=xjet
#SBATCH -q windfall
#SBATCH -t 00:05:00

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
source ../../modulefiles/build.$target

export DATA=/lfs4/HFIP/emcda/$LOGNAME/stmp/reg_tests.cycle

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export HOMEreg=/lfs4/HFIP/emcda/George.Gayno/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2

export APRUNCY="srun"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/apps/nccmp/1.8.5/intel/18.0.5.274/bin/nccmp

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
