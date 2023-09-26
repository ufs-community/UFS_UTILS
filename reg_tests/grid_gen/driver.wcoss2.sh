#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on WCOSS2.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.
#
# Invoke the script with no arguments. Several tests will be started
# in parallel.  To check the queue, type: "qstat -u USERNAME".
#
# Log output from the tests will be stored in their own LOG_FILE.  
# Once the tests have completed, a summary is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# $HOMEreg.
#
#-----------------------------------------------------------------------------

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

set -x

export WORK_DIR="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export WORK_DIR="${WORK_DIR}/reg-tests/grid-gen"
QUEUE="${QUEUE:-dev}"
PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export NCCMP=/lfs/h2/emc/global/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp

LOG_FILE=consistency.log
rm -f ${LOG_FILE}
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export APRUN_SFC="mpiexec -n 30 -ppn 30 -cpu-bind core"
export OMP_STACKSIZE=2048m
export OMP_NUM_THREADS=30 # orog code uses threads
export OMP_PLACES=cores
export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/grid_gen/baseline_data
this_dir=$PWD

ulimit -a

rm -fr $WORK_DIR

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

LOG_FILE1=${LOG_FILE}01
TEST1=$(qsub -V -o $LOG_FILE1 -e $LOG_FILE1 -q $QUEUE -A $PROJECT_CODE -l walltime=00:10:00 \
        -N c96.uniform -l select=1:ncpus=30:mem=40GB $PWD/c96.uniform.sh)

#-----------------------------------------------------------------------------
# C96 uniform grid using viirs vegetation and bnu soil type data.
#-----------------------------------------------------------------------------

LOG_FILE2=${LOG_FILE}02
TEST2=$(qsub -V -o $LOG_FILE2 -e $LOG_FILE2 -q $QUEUE -A $PROJECT_CODE -l walltime=00:15:00 \
        -N c96.viirs.bnu -l select=1:ncpus=30:mem=250GB $PWD/c96.viirs.bnu.sh)

#-----------------------------------------------------------------------------
# gfdl regional grid
#-----------------------------------------------------------------------------

LOG_FILE3=${LOG_FILE}03
TEST3=$(qsub -V -o $LOG_FILE3 -e $LOG_FILE3 -q $QUEUE -A $PROJECT_CODE -l walltime=00:07:00 \
        -N gfdl.regional -l select=1:ncpus=30:mem=40GB  $PWD/gfdl.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid (output dominant soil/vegetation type).
#-----------------------------------------------------------------------------

LOG_FILE4=${LOG_FILE}04
TEST4=$(qsub -V -o $LOG_FILE4 -e $LOG_FILE4 -q $QUEUE -A $PROJECT_CODE -l walltime=00:07:00 \
        -N esg.regional -l select=1:ncpus=30:mem=40GB $PWD/esg.regional.sh)

#-----------------------------------------------------------------------------
# ESG regional grid (output percent of each soil and vegetation type and
# the dominant category).
#-----------------------------------------------------------------------------

LOG_FILE5=${LOG_FILE}05
TEST5=$(qsub -V -o $LOG_FILE5 -e $LOG_FILE5 -q $QUEUE -A $PROJECT_CODE -l walltime=00:07:00 \
        -N esg.regional.pct.cat -l select=1:ncpus=30:mem=40GB $PWD/esg.regional.pct.cat.sh)

#-----------------------------------------------------------------------------
# Regional GSL gravity wave drag test.
#-----------------------------------------------------------------------------

LOG_FILE6=${LOG_FILE}06
TEST6=$(qsub -V -o $LOG_FILE6 -e $LOG_FILE6 -q $QUEUE -A $PROJECT_CODE -l walltime=00:07:00 \
        -N rsg.gsl.gwd -l select=1:ncpus=30:mem=40GB $PWD/regional.gsl.gwd.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:02:00 \
        -N grid_summary -l select=1:ncpus=1:mem=100MB -W depend=afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6 << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' ${LOG_FILE}* | grep -v echo > $SUM_FILE
EOF
