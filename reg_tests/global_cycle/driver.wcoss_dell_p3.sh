#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle consistency tests on WCOSS-Dell.
#
# Set $WORK_DIR to your working directory.  Set the project code 
# and queue as appropriate.
#
# Invoke the script from the command line as follows: ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log.
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

WORK_DIR="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"

PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
QUEUE="${QUEUE:-debug}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

DATA_DIR="${WORK_DIR}/reg-tests/global-cycle"

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2

export APRUNCY="mpirun -l"

export NWPROD=$PWD/../..

reg_dir=$PWD

LOG_FILE=consistency.log01
export DATA="${DATA_DIR}/test1"
export COMOUT=$DATA
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c768.fv3gfs -W 0:05 -x -n 6 \
        -M 2400 -R "span[ptile=6]" -R "affinity[core(1)]" "$PWD/C768.fv3gfs.sh"

LOG_FILE=consistency.log02
export DATA="${DATA_DIR}/test2"
export COMOUT=$DATA
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c768.lndinc -W 0:05 -x -n 6 \
        -M 2400 -R "span[ptile=6]" -R "affinity[core(1)]" "$PWD/C768.lndinc.sh"

LOG_FILE=consistency.log
bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c768.*)' "grep -a '<<<' "${LOG_FILE}*" >> summary.log"

exit
