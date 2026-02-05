#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run weight_gen consistency test on Ursa.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'nccmp' command.  The baseline file is
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

# source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1
# module use ${HOMEUFSUTILS}/modulefiles

# compiler=${compiler:-"intelllvm"}

# source ../../sorc/machine-setup.sh > /dev/null 2>&1
# module use ../../modulefiles
# module load build.${MACHINE_ID}.$compiler
# module list

# export DATA="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
test_name="weight_gen"
export DATA="${WORK_DIR}/reg_tests/${test_name}"
LOG_FILE=consistency.log

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

# UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
# export UPDATE_BASELINE

# if [ "$UPDATE_BASELINE" = "TRUE" ]; then
#   source ../get_hash.sh
# fi

# export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/weight_gen
export HOMEreg="${HOMEreg}/${test_name}"
export HOMEufs=$PWD/../..

case ${MACHINE_ID,,} in
    hercules)
        # jobkeep=$(submit_test 01 1 1 5G 0:03:00 hercules false snow.hemi snow2mdl.hemi.sh false)
        # submit_test 02 1 1 5G 0:03:00 hercules false snow.global snow2mdl.global.sh "${jobkeep}"
        submit_test 01 1 1 5G 0:03:00 hercules false weight_gen weight_gen.sh false
        ;;
    jet)
        # jobkeep=$(submit_test 01 1 1 5G 0:03:00 xjet true snow.hemi snow2mdl.hemi.sh false)
        # submit_test 02 1 1 5G 0:03:00 xjet true snow.global snow2mdl.global.sh "${jobkeep}"
        submit_test 01 1 1 5G 0:03:00 xjet true weight_gen weight_gen.sh false
        ;;
    orion)
        # jobkeep=$(submit_test 01 1 1 5G 0:03:00 orion false snow.hemi snow2mdl.hemi.sh false)
        # submit_test 02 1 1 5G 0:03:00 orion false snow.global snow2mdl.global.sh "${jobkeep}"
        submit_test 01 1 1 5G 0:03:00 orion false weight_gen weight_gen.sh false
        ;;
    ursa)
        submit_test 01 1 1 5G 0:03:00 u1-compute false weight_gen weight_gen.sh false
        ;;
    wcoss2)
        # jobkeep=$(submit_test 01 1 1 5G 0:03:00 dev false snow.hemi snow2mdl.hemi.sh false)
        # submit_test 02 1 1 5G 0:03:00 dev false snow.global snow2mdl.global.sh "${jobkeep}"
        submit_test 01 1 1 5G 0:03:00 dev false weight_gen weight_gen.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

# ./weight_gen.sh

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
