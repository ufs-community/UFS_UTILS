#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation regression tests on WCOSS-Cray.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.
#
# Invoke the script with no arguments.  A series of daily-
# chained jobs will be submitted.  To check the queue, type: "bjobs".
#
# Log output from the suite will be in LOG_FILE.  Once the suite
# has completed, a summary is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg
#
#-----------------------------------------------------------------------------

source ../../sorc/machine-setup.sh > /dev/null 2>&1
source ../../modulefiles/build.$target
module list

set -x

QUEUE="debug"
PROJECT_CODE="GFS-DEV"
export WORK_DIR=/gpfs/hps3/stmp/$LOGNAME/reg_tests.grid

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

export home_dir=$PWD/../..
LOG_FILE=regression.log
SUM_FILE=summary.log
export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
export APRUN_SFC="aprun -j 1 -n 24 -N 24"
export OMP_STACKSIZE=2048m
export OMP_NUM_THREADS=6
export machine=WCOSS_C
export KMP_AFFINITY=disabled
export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp
export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/grid_gen/baseline_data

rm -fr $WORK_DIR

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.uniform -W 0:15 -M 2400 \
        -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/c96.uniform.sh"

#-----------------------------------------------------------------------------
# C96 regional grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -W 0:10 -M 2400 \
        -w 'ended(c96.uniform)' -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/c96.regional.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 -w 'ended(c96.regional)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
