#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle regression test on WCOSS-Dell.
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

#BSUB -W 00:05
#BSUB -n 6
#BSUB -R span[ptile=6]
#BSUB -x
#BSUB -o regression.log
#BSUB -e regression.log
#BSUB -R "affinity[core(1)]"
#BSUB -M 2400
#BSUB -J glc_regt
#BSUB -q debug
#BSUB -P GFS-DEV

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export DATA="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/global-cycle"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2

export APRUNCY="mpirun -l"

export NWPROD=$PWD/../..

export COMOUT=$DATA

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
