#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on Ursa.
#
# Invoke the script from command line as follows:  ./$script
#
# Set ../rt.control variables to specify the number of tasks, memory, and walltime
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
    export OMP_NUM_THREADS_CY=2

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

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUNCY="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core --depth ${OMP_NUM_THREADS_CY}"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=${OMP_NUM_THREADS_CY}:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} "./${script}")
        jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        jobid=$(sbatch --parsable --partition="${partition}" --ntasks-per-node="${ntasks_per_node}" ${slurmflag:+"${slurmflag}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
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
    echo "${jobid}"
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

test_name="ice_blend"
LOG_FILE=consistency.log
SUM_FILE=summary.log

rm -f ${LOG_FILE}* ${SUM_FILE}

DATA="${WORK_DIR:-/work2/noaa/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice_blend"

HOMEreg="${HOMEreg}/${test_name}"
HOMEglobal=$PWD/../..
export HOMEreg HOMEglobal

case ${MACHINE_ID,,} in
    ursa)
        module load grib-util
        module load wgrib2/3.6.0
        ;;
    hercules)
        module load grib-util/1.4.0
        module load wgrib2/3.6.0
        ;;
    orion)
        module load grib-util/1.4.0
        module load wgrib2/3.6.0
        ;;
    gaeac6)
        module load grib-util/1.4.0
        module load wgrib2/3.6.0
        ;;
    wcoss2)
        module load grib_util/1.2.3
        module load wgrib2/2.0.8
        ;;
    *)
        echo "ERROR: Unsupported MACHINE_ID '${MACHINE_ID}'"
        exit 1
        ;;
esac

export COPYGB2=${COPYGB2:-${GRIB_UTIL_ROOT}/bin/copygb2}
export WGRIB2=${WGRIB2:-${wgrib2_ROOT}/bin/wgrib2}
export CNVGRIB=${CNVGRIB:-${GRIB_UTIL_ROOT}/bin/cnvgrib}
export COPYGB=${COPYGB:-${GRIB_UTIL_ROOT}/bin/copygb}

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 1 1 5G 0:01:00 hercules false false ice_blend ice_blend.sh false
        ;;
    orion)
        submit_test 01 1 1 5G 0:01:00 orion false false ice_blend ice_blend.sh false
        ;;
    ursa)
        submit_test 01 1 1 5G 0:01:00 u1-compute false false ice_blend ice_blend.sh false
        ;;
    gaeac6)
        submit_test 01 1 1 5G 0:01:00 normal c6 false ice_blend ice_blend.sh false
        ;;
    wcoss2)
        submit_test 01 1 1 5G 0:01:00 dev false false ice_blend ice_blend.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

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
