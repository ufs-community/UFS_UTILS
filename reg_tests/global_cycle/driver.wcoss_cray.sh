#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle regression test on WCOSS-Cray.
#
# Set $DATA to your working directory.  Set the project code (BSUB -P)
# and queue (BSUB -q) as appropriate.
#
# Invoke the script as follows:  cat $script | bsub
#
# Log output is placed in regression.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#BSUB -oo regression.log
#BSUB -eo regression.log
#BSUB -q debug
#BSUB -P GDAS-T2O
#BSUB -J cycle_regt
#BSUB -M 2400
#BSUB -W 00:05
#BSUB -extsched 'CRAYLINUX[]'

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export DATA="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/global-cycle"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export NODES=1

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=4

export KMP_AFFINITY=disabled

export APRUNCY="aprun -n 6 -N 6 -j 1 -d $OMP_NUM_THREADS_CY -cc depth"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
