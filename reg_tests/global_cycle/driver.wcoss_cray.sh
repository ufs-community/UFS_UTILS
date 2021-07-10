#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle consistency tests on WCOSS-Cray.
#
# Set $WORK_DIR to your working directory. Set the project code nd
# and queue as appropriate.
#
# Invoke the script as follows:  ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

WORK_DIR="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"

PROJECT_CODE="${PROJECT_CODE:-GDAS-T2O}"
QUEUE="${QUEUE:-debug}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

DATA_DIR="${WORK_DIR}/reg-tests/global-cycle"

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=4

export KMP_AFFINITY=disabled

export APRUNCY="aprun -n 6 -N 6 -j 1 -d $OMP_NUM_THREADS_CY -cc depth"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp

reg_dir=$PWD

LOG_FILE=consistency.log01
export DATA="${DATA_DIR}/test1"
export COMOUT=$DATA
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c768.fv3gfs -M 2400 -W 0:05 \
        -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/C768.fv3gfs.sh"

LOG_FILE=consistency.log02
export DATA="${DATA_DIR}/test2"
export COMOUT=$DATA
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c768.lndinc -M 2400 -W 0:05 \
        -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/C768.lndinc.sh"

LOG_FILE=consistency.log
bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c768.*)' "grep -a '<<<' "${LOG_FILE}*" >> summary.log"

exit
