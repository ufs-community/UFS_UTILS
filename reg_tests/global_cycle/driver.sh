#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle consistency test on Ursa.
#
# Set ../rt.control variables to specify the number of tasks, memory, and walltime
#
# Invoke the script from the command line as follows:  ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
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
    
    export DATA="${DATA_DIR}/test${suffix}"
    export COMOUT=$DATA

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUNCY="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core --depth ${OMP_NUM_THREADS_CY}"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=1:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"}"./${script}")
        jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        jobid=$(sbatch --parsable --partition="${partition}" --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
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

test_name="global_cycle"
HOMEreg="${HOMEreg}/${test_name}"

# EXPORTED VARIABLES
DATA_DIR="${WORK_DIR}/reg-tests/${test_name}"
OMP_NUM_THREADS_CY=2
OMP_PLACES=cores
NWPROD="${WORK_DIR}/UFS_UTILS"

export DATA_DIR OMP_NUM_THREADS_CY OMP_PLACES HOMEreg NWPROD
LOG_FILE=consistency.log
SUM_FILE=summary.log
reg_dir=$PWD

rm -f ${LOG_FILE}* ${SUM_FILE}

declare -a TEST_IDS=()

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 6 1 50G 0:05:00 hercules false C768.fv3gfs C768.fv3gfs.sh false
        #submit_test 02 6 1 50G 0:05:00 hercules false C192.gsi_lndincsoilnoahmp C192.gsi_lndincsoilnoahmp.sh false
        submit_test 03 6 1 50G 0:05:00 hercules false C768.lndincsnow C768.lndincsnow.sh false
        submit_test 04 6 1 50G 0:05:00 hercules false C48.noahmp.coupled C48.noahmp.coupled.sh false
        submit_test 05 6 1 50G 0:05:00 hercules false C192.jedi_lndincsoilnoahmp C192.jedi_lndincsoilnoahmp.sh false
        submit_test 06 6 1 50G 0:05:00 hercules false C192.gsitile_lndincsoilnoahmp C192.gsitile_lndincsoilnoahmp.sh false
        ;;
    orion)
        submit_test 01 6 1 50G 0:05:00 orion false C768.fv3gfs C768.fv3gfs.sh false
        #submit_test 02 6 1 50G 0:05:00 orion false C192.gsi_lndincsoilnoahmp C192.gsi_lndincsoilnoahmp.sh false
        submit_test 03 6 1 50G 0:05:00 orion false C768.lndincsnow C768.lndincsnow.sh false
        submit_test 04 6 1 50G 0:05:00 orion false C48.noahmp.coupled C48.noahmp.coupled.sh false
        submit_test 05 6 1 50G 0:05:00 orion false C192.jedi_lndincsoilnoahmp C192.jedi_lndincsoilnoahmp.sh false
        submit_test 06 6 1 50G 0:05:00 orion false C192.gsitile_lndincsoilnoahmp C192.gsitile_lndincsoilnoahmp.sh false
        ;;
    ursa)
        submit_test 01 6 1 50G 0:05:00 u1-compute false C768.fv3gfs C768.fv3gfs.sh false
        #submit_test 02 6 1 50G 0:05:00 u1-compute false C192.gsi_lndincsoilnoahmp C192.gsi_lndincsoilnoahmp.sh false
        submit_test 03 6 1 50G 0:05:00 u1-compute false C768.lndincsnow C768.lndincsnow.sh false
        submit_test 04 6 1 50G 0:05:00 u1-compute false C48.noahmp.coupled C48.noahmp.coupled.sh false
        submit_test 05 6 1 50G 0:05:00 u1-compute false C192.jedi_lndincsoilnoahmp C192.jedi_lndincsoilnoahmp.sh false
        submit_test 06 6 1 50G 0:05:00 u1-compute false C192.gsitile_lndincsoilnoahmp C192.gsitile_lndincsoilnoahmp.sh false
        ;;
    wcoss2)
        submit_test 01 12 1 15G 0:05:00 dev false C768.fv3gfs C768.fv3gfs.sh false
        #submit_test 02 12 1 15G 0:05:00 dev false C192.gsi_lndincsoilnoahmp C192.gsi_lndincsoilnoahmp.sh false
        submit_test 03 12 1 15G 0:05:00 dev false C768.lndincsnow C768.lndincsnow.sh false
        submit_test 04 12 1 15G 0:05:00 dev false C48.noahmp.coupled C48.noahmp.coupled.sh false
        submit_test 05 12 1 15G 0:05:00 dev false C192.jedi_lndincsoilnoahmp C192.jedi_lndincsoilnoahmp.sh false
        submit_test 06 12 1 15G 0:05:00 dev false C192.gsitile_lndincsoilnoahmp C192.gsitile_lndincsoilnoahmp.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

if [[ "${SCHEDULER}" == "pbs" ]]; then
  (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N cycle_summary -l select=1:ncpus=1:mem=100MB -W depend="afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd $reg_dir
grep -a '^<<<' ${LOG_FILE}* | grep -v echo > summary.log
EOF
  ) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
  (sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q $QUEUE \
      -d "afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd $reg_dir
grep -a '^<<<' ${LOG_FILE}*  > summary.log
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
