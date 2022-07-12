#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on Jet.
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
module use ../../modulefiles
module load build.$target.intel
module list

set -x

QUEUE="${QUEUE:-windfall}"
PROJECT_CODE="${PROJECT_CODE:-emcda}"
export WORK_DIR="${WORK_DIR:-/lfs4/HFIP/emcda/$LOGNAME/stmp}"
export WORK_DIR="${WORK_DIR}/reg-tests/grid-gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

LOG_FILE=consistency.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export APRUN_SFC=srun
export OMP_STACKSIZE=2048m
export HOMEreg=/lfs4/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/grid_gen/baseline_data

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
# C96 uniform grid using viirs vegetation type data.
#-----------------------------------------------------------------------------

TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.viirs.vegt \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c96.viirs.vegt.sh)

#-----------------------------------------------------------------------------
# gfdl regional grid
#-----------------------------------------------------------------------------

TEST3=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J gfdl.regional \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST2 ./gfdl.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid
#-----------------------------------------------------------------------------

TEST4=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST3 ./esg.regional.sh)

#-----------------------------------------------------------------------------
# Regional GSL gravity wave drag.
#-----------------------------------------------------------------------------

TEST5=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J reg.gsl.gwd \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST4 ./regional.gsl.gwd.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST5 << EOF
#!/bin/bash
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF
