#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on WCOSS-Dell.
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
module load git
module list

set -x

export DATA_ROOT="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"
export DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

PROJECT_CODE=${PROJECT_CODE:-"GFS-DEV"}
QUEUE=${QUEUE:-"debug"}

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib
export WGRIB2=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib2

LOG_FILE=consistency.log
SUM_FILE=summary.log

rm -fr $DATA_ROOT

export DATA=$DATA_ROOT/test.ops
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J snow.ops -W 0:02  \
        -R "affinity[core(1)]" "$PWD/snow2mdl.ops.sh"

export DATA=$DATA_ROOT/test.global
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J snow.global -W 0:02  \
        -R "affinity[core(1)]" -w 'ended(snow.ops)' "$PWD/snow2mdl.global.sh"

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(snow.global)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
