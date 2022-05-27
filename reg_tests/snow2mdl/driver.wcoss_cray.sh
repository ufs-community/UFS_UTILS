#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on WCOSS-Cray.
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

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

DATA_ROOT="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

PROJECT_CODE=${PROJECT_CODE:-GFS-DEV}
QUEUE=${QUEUE:-dev}

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.2/exec/wgrib
export WGRIB2=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.2/exec/wgrib2

rm -fr $DATA_ROOT

LOG_FILE="consistency.log"
SUM_FILE="summary.log"

# Test the ops function of snow2mdl.

export DATA=$DATA_ROOT/test.ops
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J snow.ops -W 0:02 \
        -R "rusage[mem=2000]" "$PWD/snow2mdl.ops.sh"

# Test the afwa global snow data.

export DATA=$DATA_ROOT/test.global
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J snow.global -W 0:02 \
        -R "rusage[mem=2000]" -w 'ended(snow.ops)' "$PWD/snow2mdl.global.sh"

# Create a summary file.

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(snow.global)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
