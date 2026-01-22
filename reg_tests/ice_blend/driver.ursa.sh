#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on Ursa.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
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

RT_DIR=${RT_DIR:-${PWD}/..}

if [[ ! -v PID_LIST ]]; then
  waitlocal=true
fi

if [[ -f "${RT_DIR}/rt.control" ]]; then
    source "${RT_DIR}/rt.control"
else
    echo "ERROR: Cannot find rt.control script"
    exit 1
fi

source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${HOMEUFSUTILS}/modulefiles

compiler=${compiler:-"intelllvm"}

# source ../../sorc/machine-setup.sh > /dev/null 2>&1
# module use ../../modulefiles
module load build.${MACHINE_ID}.$compiler
module load grib-util
module load wgrib2/3.6.0
module list

export CNVGRIB=${GRIB_UTIL_ROOT}/bin/cnvgrib
export COPYGB=${GRIB_UTIL_ROOT}/bin/copygb
export COPYGB2=${GRIB_UTIL_ROOT}/bin/copygb2

export DATA="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
export UPDATE_BASELINE

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/ice_blend
export HOMEgfs=$PWD/../..

./ice_blend.sh

if [[ "${waitlocal}" == "true" ]]; then
  sleep_time=0
  echo "Waiting for ice_blend test to complete..."
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       mail -s "UFS_UTILS Consistency Test ICE_BLEND timed out on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
       exit 1
    fi
  done
  mail -s "UFS_UTILS Consistency Test ICE_BLEND COMPLETED on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
fi

exit 0
