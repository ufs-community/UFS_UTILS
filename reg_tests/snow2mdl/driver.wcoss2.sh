#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on WCOSS2.
#
# Set $DATA_ROOT to your working directory.  Set the project code 
# and queue as appropriate.
#
# Invoke the script as follows:  ./$script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module load grib_util/1.2.2
module load wgrib2/2.0.8
module list

set -x

export DATA_ROOT="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

PROJECT_CODE=${PROJECT_CODE:-"GFS-DEV"}
QUEUE=${QUEUE:-"dev"}

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..

LOG_FILE=consistency.log
SUM_FILE=summary.log

rm -fr $DATA_ROOT

#-----------------------------------------------------------------------------
# Test GFS ops snow.
#-----------------------------------------------------------------------------

export DATA=$DATA_ROOT/test.ops
TEST1=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l select=1:ncpus=1:mem=2500MB \
        -N snow.ops -l walltime=00:03:00 $PWD/snow2mdl.ops.sh)

#-----------------------------------------------------------------------------
# Test afwa global snow.
#-----------------------------------------------------------------------------

export DATA=$DATA_ROOT/test.global
TEST2=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l select=1:ncpus=1:mem=2500MB \
        -N snow.global -l walltime=00:03:00 -W depend=afterok:$TEST1 $PWD/snow2mdl.global.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

this_dir=$PWD
qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N snow_summary -l select=1:ncpus=1:mem=100MB -W depend=afterok:$TEST2 << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' $LOG_FILE | grep -v echo > $SUM_FILE
EOF
exit 0
