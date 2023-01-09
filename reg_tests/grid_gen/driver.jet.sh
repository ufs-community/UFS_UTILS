#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on Jet.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.  To see which projects you are authorized to use,
# type "account_params".
#
# Invoke the script with no arguments.  Several tests will
# be submitted to run in parallel. To check the queue, type:
# "squeue -u USERNAME".
#
# Log output from each test will be placed in its own LOG_FILE.
# Once the suite has completed, a summary is placed in SUM_FILE.
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

QUEUE="${QUEUE:-batch}"
PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
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

LOG_FILE1=${LOG_FILE}01
TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
      --partition=xjet -o $LOG_FILE1 -e $LOG_FILE1 ./c96.uniform.sh)

#-----------------------------------------------------------------------------
# C96 uniform grid using viirs vegetation type data.
#-----------------------------------------------------------------------------

LOG_FILE2=${LOG_FILE}02
TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.viirs.vegt \
      --partition=xjet -o $LOG_FILE2 -e $LOG_FILE2 ./c96.viirs.vegt.sh)

#-----------------------------------------------------------------------------
# gfdl regional grid
#-----------------------------------------------------------------------------

LOG_FILE3=${LOG_FILE}03
TEST3=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J gfdl.regional \
      --partition=xjet -o $LOG_FILE3 -e $LOG_FILE3 ./gfdl.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid (output dominate soil/vegetation type).
#-----------------------------------------------------------------------------

LOG_FILE4=${LOG_FILE}04
TEST4=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional \
      --partition=xjet -o $LOG_FILE4 -e $LOG_FILE4 ./esg.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid (output percent of each soil and vegetation type and
# the dominate category).
#-----------------------------------------------------------------------------

LOG_FILE5=${LOG_FILE}05
TEST5=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional.pct.cat \
      --partition=xjet -o $LOG_FILE5 -e $LOG_FILE5 ./esg.regional.pct.cat.sh)

#-----------------------------------------------------------------------------
# Regional GSL gravity wave drag.
#-----------------------------------------------------------------------------

LOG_FILE6=${LOG_FILE}06
TEST6=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J reg.gsl.gwd \
      --partition=xjet -o $LOG_FILE6 -e $LOG_FILE6 ./regional.gsl.gwd.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       -q $QUEUE -d afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
EOF
