#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on Orion.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.  To see which projects you are authorized to use,
# type "saccount_params".
#
# Invoke the script with no arguments.  A series of daily-
# chained jobs will be submitted.  To check the queue, type:
# "squeue -u $LOGNAME".
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
ulimit -s unlimited

export WORK_DIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
export WORK_DIR="${WORK_DIR}/reg-tests/grid-gen"
QUEUE="${QUEUE:-batch}"
PROJECT_CODE=${PROJECT_CODE:-fv3-cpu}

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export APRUN_SFC=srun
export OMP_STACKSIZE=2048m
export OMP_NUM_THREADS=24
export machine=ORION
export HOMEreg=/work/noaa/da/ggayno/save/ufs_utils.git/reg_tests/grid_gen/baseline_data

rm -fr $WORK_DIR

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
      -o $LOG_FILE -e $LOG_FILE ./c96.uniform.sh)

#-----------------------------------------------------------------------------
# GFDL regional grid
#-----------------------------------------------------------------------------

TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J gfdl.regional \
      --open-mode=append -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./gfdl.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid
#-----------------------------------------------------------------------------

TEST3=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional \
      --open-mode=append -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST2 ./esg.regional.sh)

#-----------------------------------------------------------------------------
# Regional grid with GSL gravity wave drag fields.
#-----------------------------------------------------------------------------

TEST4=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J reg.gsl.gwd \
      --open-mode=append -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST3 ./regional.gsl.gwd.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST4 << EOF
#!/bin/bash
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF
