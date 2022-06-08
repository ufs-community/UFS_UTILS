#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency test on S4.
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
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

DATA_ROOT="${WORK_DIR:-/scratch/short/users/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

PROJECT_CODE="${PROJECT_CODE:-star}"
QUEUE="${QUEUE:-s4}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/data/users/dhuber/save/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/wgrib
export WGRIB2=/data/prod/hpc-stack/intel-2022.1/wgrib2/2.0.8/bin/wgrib2

# The first test mimics GFS OPS.

export DATA="${DATA_ROOT}/test.ops"
TEST1=$(sbatch --parsable -J snow.ops -A ${PROJECT_CODE} -o consistency.log -e consistency.log \
      --ntasks=1 -q ${QUEUE} -t 00:03:00 ./snow2mdl.ops.sh)

# The second test is for the new AFWA global GRIB2 data.

export DATA="${DATA_ROOT}/test.global"
TEST2=$(sbatch --parsable -J snow.global -A ${PROJECT_CODE} -o consistency.log -e consistency.log \
      --ntasks=1 -q ${QUEUE} -t 00:03:00 -d afterok:$TEST1 ./snow2mdl.global.sh)

# Create summary file.

sbatch --nodes=1 -t 0:01:00 -A ${PROJECT_CODE} -J snow_summary -o consistency.log -e consistency.log \
       --open-mode=append -q ${QUEUE} -d afterok:$TEST2 << EOF
#!/bin/bash
grep -a '<<<' consistency.log  > summary.log
EOF

exit 0
