#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on WCOSS-Dell.
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
module use ../../modulefiles
module load build.$target.intel
module list

set -x

QUEUE="${QUEUE:-debug}"
PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
export WORK_DIR="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"
export WORK_DIR="${WORK_DIR}/reg-tests/grid-gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export APRUN_SFC="mpirun -l"
export OMP_STACKSIZE=2048m
export machine=WCOSS_DELL_P3
export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/grid_gen/baseline_data
export OMP_NUM_THREADS=24

rm -fr $WORK_DIR

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.uniform -W 0:15 -x -n 24 \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/c96.uniform.sh"

#-----------------------------------------------------------------------------
# GFDL regional grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J gfdl.regional -W 0:10 -x -n 24 -w 'ended(c96.uniform)' \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/gfdl.regional.sh"

#-----------------------------------------------------------------------------
# ESG regional grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J esg.regional -W 0:10 -x -n 24 -w 'ended(gfdl.regional)' \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/esg.regional.sh"

#-----------------------------------------------------------------------------
# Regional GSL gravity wave drag.
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J reg.gsl.gwd -W 0:08 -x -n 24 -w 'ended(esg.regional)' \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/regional.gsl.gwd.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(reg.gsl.gwd)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"
