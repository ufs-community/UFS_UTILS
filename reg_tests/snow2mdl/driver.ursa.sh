#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on Ursa.
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

compiler=${compiler:-"intelllvm"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module load grib-util
module load wgrib2
module load prod_util
module list

WGRIB=${GRIB_UTIL_ROOT}/bin/wgrib

export WGRIB
export WGRIB2

DATA_ROOT="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
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

export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..

# The first test uses hemispheric afwa/airforce data, as was done in OPS.

export DATA="${DATA_ROOT}/test.hemi"
TEST1=$(sbatch --parsable -J snow.hemi -A ${PROJECT_CODE} -o consistency.log -e consistency.log \
      --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 ./snow2mdl.hemi.sh)

# The second test mimics current OPS, which uses global afwa/airforce data.

export DATA="${DATA_ROOT}/test.global"
TEST2=$(sbatch --parsable -J snow.global -A ${PROJECT_CODE} -o consistency.log -e consistency.log \
      --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 -d afterok:$TEST1 ./snow2mdl.global.sh)

# Create summary file.

sbatch --nodes=1 -t 0:01:00 -A ${PROJECT_CODE} -J snow_summary -o consistency.log -e consistency.log \
       --open-mode=append -q ${QUEUE} -d afterok:$TEST2 << EOF
#!/bin/bash
grep -a '<<<' consistency.log  > summary.log
EOF

exit 0
