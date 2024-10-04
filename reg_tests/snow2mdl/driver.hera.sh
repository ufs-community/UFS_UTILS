#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on Hera.
#
# Set $DATA_ROOT to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script from the command line as follows:  ./$script
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
module load grib-util
module load wgrib2/3.1.1
module list

# Because of a bug in the grib-util module, need to construct this
# variable.
WGRIB=${grib_util_ROOT}/bin/wgrib

export WGRIB
export WGRIB2

DATA_ROOT="${WORK_DIR:-/scratch2/NCEPDEV/stmp1/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

rm -fr $DATA_ROOT

PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..

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
