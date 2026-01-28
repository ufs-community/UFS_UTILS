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

RT_DIR=${RT_DIR:-${PWD}/..}

notlocal=${notlocal:-false}
if [[ ${notlocal} == "false" ]]; then
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

#source ../../sorc/machine-setup.sh > /dev/null 2>&1
#module use ../../modulefiles
module load build.${MACHINE_ID,,}.$compiler
module list

DATA_ROOT="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

rm -fr $DATA_ROOT

PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
export UPDATE_BASELINE

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..

# The first test uses hemispheric afwa/airforce data, as was done in OPS.

LOG_FILE=consistency.log
SUM_FILE=summary.log
if [[ -f ${LOG_FILE} ]]; then
  rm -f ${LOG_FILE}*
fi

if [[ -f ${SUM_FILE} ]]; then
  rm -f ${SUM_FILE}
fi

export DATA="${DATA_ROOT}/test.hemi"
TEST1=$(sbatch --parsable -J snow.hemi -A ${PROJECT_CODE} -o ${LOG_FILE}01 -e ${LOG_FILE}01 \
      --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 ./snow2mdl.hemi.sh)

# The second test mimics current OPS, which uses global afwa/airforce data.

export DATA="${DATA_ROOT}/test.global"
TEST2=$(sbatch --parsable -J snow.global -A ${PROJECT_CODE} -o ${LOG_FILE}02 -e ${LOG_FILE}02 \
      --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 -d afterok:$TEST1 ./snow2mdl.global.sh)

# Create summary file.

(sbatch --nodes=1 -t 0:01:00 -A ${PROJECT_CODE} -J snow_summary -o ${LOG_FILE} -e ${LOG_FILE} \
       --open-mode=append -q ${QUEUE} -d afterok:$TEST2 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}*  > summary.log
EOF
) &

if [[ "${waitlocal}" == "true" ]]; then
  sleep_time=0
  echo "Waiting for snow2mdl tests to complete..."
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       mail -s "UFS_UTILS Consistency Test SNOW2MDL timed out on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
       exit 1
    fi
  done
  mail -s "UFS_UTILS Consistency Test SNOW2MDL COMPLETED on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
fi
exit 0
