#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on Jet.
#
# Set $DATA_ROOT to your working directory.  Set the project code and
# and queue as appropriate.
#
# Invoke the script as follows:  ./$script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------


set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module load wgrib2/2.0.8
set +x
module list
set -x

DATA_ROOT="${WORK_DIR:-/lfs4/HFIP/emcda/$LOGNAME/stmp}"
DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/lfs4/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/apps/wgrib/1.8.1.0b/bin/wgrib

rm -fr $DATA_ROOT

# This tests the OPS GFS snow processing.

export DATA="${DATA_ROOT}/test.ops"
TEST1=$(sbatch --parsable --nodes=1 --partition=xjet --time 0:02 -J snow.ops -o consistency.log \
        -e consistency.log -A $PROJECT_CODE -q $QUEUE ./snow2mdl.ops.sh)

# Test the new global afwa grib2 data.

export DATA="${DATA_ROOT}/test.global"
TEST2=$(sbatch --parsable --nodes=1 --partition=xjet --time 0:02 -J snow.global -o consistency.log \
        -e consistency.log -A $PROJECT_CODE -q $QUEUE -d afterok:$TEST1 ./snow2mdl.global.sh)

# Create summary file.

sbatch --nodes=1 --partition=xjet -t 0:01:00 -A $PROJECT_CODE -J snow.summary -o consistency.log \
       -e consistency.log --open-mode=append -q $QUEUE -d afterok:$TEST2 << EOF
#!/bin/bash
grep -a '<<<' consistency.log  > summary.log
EOF
exit 0
