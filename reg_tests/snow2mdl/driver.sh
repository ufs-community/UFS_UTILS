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

submit_test() {
    local suffix="$1"; shift
    local ntasks_per_node="$1"; shift
    local nodes="$1"; shift
    local mem="$1"; shift
    local walltime="$1"; shift
    local partition="$1"; shift
    local exclusive="$1"; shift
    local jobname="$1"; shift
    local script="$1"; shift
    local waitonjobid="$1"; shift

    local logfile="${LOG_FILE}${suffix}"
    export OMP_NUM_THREADS_CY=2

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    export DATA="${DATA_ROOT}/test${suffix}"

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUNCY="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core --depth ${OMP_NUM_THREADS_CY}"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=1:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} "./${script}")
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        jobid=$(sbatch --parsable --partition="${partition}" --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
               -A "${PROJECT_CODE}" -q "${QUEUE}" -J "${jobname}" --open-mode=append ${exclusive_flag:+"${exclusive_flag}"} \
               ${dep_flag_slurm:+"${dep_flag_slurm}"} -o "${logfile}" -e "${logfile}" "./${script}")
    else
        echo "Error: Unsupported scheduler '${SCHEDULER}'"
        exit 1
    fi
    if [[ "${jobid}" == "" ]]; then
        echo "Error submitting job to slurm scheduler"
        exit 1
    fi
    TEST_IDS+=(":${jobid}")
    echo ${jobid}
}

# compiler=${compiler:-"intelllvm"}

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

test_name="snow2mdl"
# source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1
# module use ${HOMEUFSUTILS}/modulefiles

# #source ../../sorc/machine-setup.sh > /dev/null 2>&1
# #module use ../../modulefiles
# module load build.${MACHINE_ID,,}.$compiler
# module list

DATA_ROOT="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/${test_name}"

rm -fr $DATA_ROOT

# PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
# QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
export UPDATE_BASELINE

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

# export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
HOMEreg="${HOMEreg}/${test_name}"
export HOMEgfs=$PWD/../..

# The first test uses hemispheric afwa/airforce data, as was done in OPS.

LOG_FILE=consistency.log
SUM_FILE=summary.log

declare -a TEST_IDS=()

if [[ -f ${LOG_FILE} ]]; then
  rm -f ${LOG_FILE}*
fi

if [[ -f ${SUM_FILE} ]]; then
  rm -f ${SUM_FILE}
fi

# export DATA="${DATA_ROOT}/test.hemi"
# TEST1=$(sbatch --parsable -J snow.hemi -A ${PROJECT_CODE} -o ${LOG_FILE}01 -e ${LOG_FILE}01 \
#       --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 ./snow2mdl.hemi.sh)

# # The second test mimics current OPS, which uses global afwa/airforce data.

# export DATA="${DATA_ROOT}/test.global"
# TEST2=$(sbatch --parsable -J snow.global -A ${PROJECT_CODE} -o ${LOG_FILE}02 -e ${LOG_FILE}02 \
#       --ntasks=1 --mem=5GB -q ${QUEUE} -t 00:03:00 -d afterok:$TEST1 ./snow2mdl.global.sh)

case ${MACHINE_ID,,} in
    hercules)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 hercules false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 hercules false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    jet)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 xjet true snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 xjet true snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    orion)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 orion false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 orion false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    ursa)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 u1-compute false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 u1-compute false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    wcoss2)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 dev false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 dev false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

# Create summary file.
if [[ "${SCHEDULER}" == "pbs" ]]; then
  this_dir=$PWD
  (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
          -N snow_summary -l select=1:ncpus=1:mem=100MB -W depend="afterok$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' $LOG_FILE | grep -v echo > $SUM_FILE
EOF
  ) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
  (sbatch --nodes=1 -t 0:01:00 -A ${PROJECT_CODE} -J snow_summary -o ${LOG_FILE} -e ${LOG_FILE} \
        --open-mode=append -q ${QUEUE} -d "afterok$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}*  > summary.log
EOF
  ) &
else
    echo "Error: Unsupported scheduler '${SCHEDULER}'"
    exit 1
fi

sleep_time=0
echo "Waiting for ${test_name^^} tests to complete..."
while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
        if [[ "${waitlocal}" == "true" ]]; then
            mail -s "UFS_UTILS Consistency Test ${test_name^^} timed out on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
            exit 1
        fi
    fi
done
if [[ "${waitlocal}" == "true" ]]; then
    mail -s "UFS_UTILS Consistency Test ${test_name^^} COMPLETED on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
fi
exit 0
