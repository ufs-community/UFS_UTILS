#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on Ursa.
#
# Set ../rt.control variables to specify the number of tasks, memory, and walltime
#
# Invoke the script with no arguments.  A set of tests will
# be submitted to run in parallel.  To check the queue, type:
# "squeue -u USERNAME".
#
# Log output from each tests will be in its own LOG_FILE.  Once 
# the tests have completed, a summary is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg.
#
#-----------------------------------------------------------------------------

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
    OMP_NUM_THREADS=${ntasks_per_node}

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    export DATA="${DATA_ROOT}/test${suffix}"

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUNCY="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core"
        export APRUN_SFC=${APRUNCY}
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=${OMP_NUM_THREADS}:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} "./${script}")
        jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        export APRUN_SFC=${APRUNCY}
        jobid=$(sbatch --parsable --partition="${partition}" --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
               -A "${PROJECT_CODE}" -q "${QUEUE}" -J "${jobname}" --open-mode=append ${exclusive_flag:+"${exclusive_flag}"} \
               --export=ALL,OMP_NUM_THREADS=${OMP_NUM_THREADS} ${dep_flag_slurm:+"${dep_flag_slurm}"} \
               -o "${logfile}" -e "${logfile}" "./${script}")
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

set -x

test_name="grid_gen"
export WORK_DIR="${WORK_DIR}/reg-tests/${test_name}"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

LOG_FILE=consistency.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export this_dir=$PWD
export OMP_STACKSIZE=2048m
HOMEreg="${HOMEreg}/${test_name}"

ulimit -a

declare -a TEST_IDS=()
rm -f ${LOG_FILE}* ${SUM_FILE}
rm -fr "${WORK_DIR}"

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 24 1 50G 0:15:00 hercules false c96.uniform c96.uniform.sh false
        submit_test 02 15 2 300G 0:15:00 hercules false c96.viirs.bnu c96.viirs.bnu.sh false
        submit_test 03 24 1 50G 0:10:00 hercules false gfdl.regional gfdl.regional.sh false
        submit_test 04 24 1 50G 0:10:00 hercules false esg.regional esg.regional.sh false
        submit_test 05 24 1 50G 0:10:00 hercules false esg.regional.pct.cat esg.regional.pct.cat.sh false
        submit_test 06 12 1 50G 0:10:00 hercules false reg.gsl.gwd.12 regional.gsl.gwd.sh false
        submit_test 07 24 1 50G 0:10:00 hercules false reg.gsl.gwd.24 regional.gsl.gwd.sh false
        ;;
    orion)
        submit_test 01 24 1 50G 0:20:00 orion false c96.uniform c96.uniform.sh false
        submit_test 02 15 2 96G 0:20:00 orion false c96.viirs.bnu c96.viirs.bnu.sh false
        submit_test 03 24 1 50G 0:10:00 orion false gfdl.regional gfdl.regional.sh false
        submit_test 04 24 1 50G 0:10:00 orion false esg.regional esg.regional.sh false
        submit_test 05 24 1 50G 0:10:00 orion false esg.regional.pct.cat esg.regional.pct.cat.sh false
        submit_test 06 12 1 50G 0:10:00 orion false reg.gsl.gwd.12 regional.gsl.gwd.sh false
        submit_test 07 24 1 50G 0:10:00 orion false reg.gsl.gwd.24 regional.gsl.gwd.sh false
        ;;
    ursa)
        submit_test 01 24 1 50G 0:15:00 u1-compute false c96.uniform c96.uniform.sh false
        submit_test 02 12 2 300G 0:15:00 u1-compute false c96.viirs.bnu c96.viirs.bnu.sh false
        submit_test 03 24 1 50G 0:07:00 u1-compute false gfdl.regional gfdl.regional.sh false
        submit_test 04 24 1 50G 0:07:00 u1-compute false esg.regional esg.regional.sh false
        submit_test 05 24 1 50G 0:07:00 u1-compute false esg.regional.pct.cat esg.regional.pct.cat.sh false
        submit_test 06 12 1 50G 0:07:00 u1-compute false reg.gsl.gwd.12 regional.gsl.gwd.sh false
        submit_test 07 24 1 50G 0:07:00 u1-compute false reg.gsl.gwd.24 regional.gsl.gwd.sh false
        ;;
    wcoss2)
        submit_test 01 30 1 40G 0:15:00 dev false c96.uniform c96.uniform.sh false
        submit_test 02 30 1 250G 0:15:00 dev false c96.viirs.bnu c96.viirs.bnu.sh false
        submit_test 03 30 1 40G 0:07:00 dev false gfdl.regional gfdl.regional.sh false
        submit_test 04 30 1 40G 0:07:00 dev false esg.regional esg.regional.sh false
        submit_test 05 30 1 40G 0:07:00 dev false esg.regional.pct.cat esg.regional.pct.cat.sh false
        submit_test 06 15 1 40G 0:07:00 dev false reg.gsl.gwd.12 regional.gsl.gwd.sh false
        submit_test 07 30 1 40G 0:07:00 dev false reg.gsl.gwd.24 regional.gsl.gwd.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------
if [[ "${SCHEDULER}" == "pbs" ]]; then
  (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:02:00 \
        -N grid_summary -l select=1:ncpus=1:mem=100MB -W depend=afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]') << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' ${LOG_FILE}* | grep -v echo > $SUM_FILE
EOF
  ) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
  (sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]') << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
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
