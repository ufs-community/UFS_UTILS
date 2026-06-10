#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on Ursa.
#
# Set ../rt.control variables to specify the number of tasks, memory, and walltime
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
    local slurmcluster="$1"; shift
    local exclusive="$1"; shift
    local jobname="$1"; shift
    local script="$1"; shift
    local waitonjobid="$1"; shift

    local logfile="${LOG_FILE}${suffix}"
    export OMP_NUM_THREADS=1

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${slurmcluster}" != "false" ]]; then
        slurmflag="--clusters=${slurmcluster}"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    export DATA="${DATA_ROOT}/test${suffix}"

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUNCY="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core --depth ${OMP_NUM_THREADS}"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=${OMP_NUM_THREADS}:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} "./${script}")
        jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        jobid=$(sbatch --parsable --partition="${partition}" ${slurmflag:+"${slurmflag}"} --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
               -A "${PROJECT_CODE}" -q "${QUEUE}" -J "${jobname}" --open-mode=append ${exclusive_flag:+"${exclusive_flag}"} \
               ${dep_flag_slurm:+"${dep_flag_slurm}"} -o "${logfile}" -e "${logfile}" "./${script}")
        jobid=${jobid%.*}
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

RT_DIR=${RT_DIR:-${PWD}/..}

notlocal=${notlocal:-false}
waitlocal=false
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

DATA_ROOT="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
DATA_ROOT="${DATA_ROOT}/reg-tests/${test_name}"

rm -fr $DATA_ROOT

UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
export UPDATE_BASELINE

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

HOMEreg="${HOMEreg}/${test_name}"
HOMEglobal=$PWD/../..
export HOMEreg HOMEglobal

# The first test uses hemispheric afwa/airforce data, as was done in OPS.
LOG_FILE=consistency.log
SUM_FILE=summary.log

declare -a TEST_IDS=()

rm -f ${LOG_FILE}* ${SUM_FILE}

case ${MACHINE_ID,,} in
    hercules)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 hercules false false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 hercules false false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    orion)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 orion false false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 orion false false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    ursa)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 u1-compute false false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 u1-compute false false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    gaeac6)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 batch c6 false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 batch c6 false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    wcoss2)
        jobkeep=$(submit_test 01 1 1 5G 0:03:00 dev false false snow.hemi snow2mdl.hemi.sh false)
        submit_test 02 1 1 5G 0:03:00 dev false false snow.global snow2mdl.global.sh "${jobkeep}"
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

# Create summary file from logs.
this_dir=$PWD
if [[ "${SCHEDULER}" == "pbs" ]]; then
  (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
          -N snow_summary -l select=1:ncpus=1:mem=100MB -W depend="afterok$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' $LOG_FILE* | grep -v echo > $SUM_FILE
EOF
  ) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
  (sbatch --nodes=1 -t 0:01:00 -A ${PROJECT_CODE} ${slurmflag:+"${slurmflag}"} -J snow_summary -o ${LOG_FILE} -e ${LOG_FILE} \
        --open-mode=append -q ${QUEUE} -d "afterok$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd ${this_dir}
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
