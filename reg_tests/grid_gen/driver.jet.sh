#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation regression tests on Jet.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.  To see which projects you are authorized to use,
# type "account_params".
#
# Invoke the script with no arguments.  A series of daily-
# chained jobs will be submitted.  To check the queue, type:
# "squeue -u USERNAME".
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

set -x

QUEUE="windfall"
PROJECT_CODE="emcda"
export WORK_DIR=/lfs4/HFIP/emcda/$LOGNAME/stmp/reg_tests.grid

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

LOG_FILE=regression.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export APRUN_SFC=srun
export OMP_STACKSIZE=2048m
export machine=JET
export NCCMP=/apps/nccmp/1.8.5/intel/18.0.5.274/bin/nccmp
export HOMEreg=/lfs4/HFIP/emcda/George.Gayno/reg_tests/grid_gen/baseline_data

ulimit -a
ulimit -s unlimited

rm -fr $WORK_DIR

export OMP_NUM_THREADS=24

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE ./c96.uniform.sh)

#-----------------------------------------------------------------------------
# C96 regional grid
#-----------------------------------------------------------------------------

TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST2 << EOF
#!/bin/sh
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF
