#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run regrid_sfc consistency tests.
#
# Set $WORK_DIR to your working directory. 
# Set the $PROJECT_CODE and $QUEUE as appropriate.
#
# Invoke the script from command line as follows:  ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility. The baseline files are
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
    export PARTITION="${partition}"

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUN_REGRID="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core --depth ${OMP_NUM_THREADS_CY}"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=1:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"}"./${script}")
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUN_REGRID="srun"
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
}

test_name="regrid_sfc"

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
# compiler=${compiler:-intelllvm}
# if [[ "${compiler}" == "intelllvm" ]]; then
#   if [[ ! -f ${HOMEUFSUTILS}/modulefiles/build.$MACHINE_ID.$compiler.lua ]];then
#      set +x
#      echo "IntelLLVM not available. Will use Intel Classic."
#      set -x
#     compiler=intel
#   fi
# fi

# module load "build.${MACHINE_ID}.${compiler}"
if [[ "${MACHINE_ID}" == "wcoss2" ]];then
  module load nccmp-D/1.9.0.1
fi
# set +x
# module list
# set -x



##UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
##export UPDATE_BASELINE

if [[ "$UPDATE_BASELINE" == "TRUE" ]]; then
  if [[ -f "${RT_DIR}/get_hash.sh" ]]; then
    source "${RT_DIR}/get_hash.sh"
  else
    echo "ERROR: Cannot find detect_machine.sh script"
    exit 1
  fi
fi

export HOMEreg="${HOMEreg}/${test_name}"

# if [[ "${MACHINE_ID}" == "jet" ]];then
#  export WORK_DIR="${WORK_DIR:-/lfs5/HFIP/emcda/$LOGNAME/stmp}"
#  PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
#  QUEUE="${QUEUE:-batch}"
#  export HOMEreg=${HOMEreg}/${test_name}
  # export APRUN_REGRID=srun
  # PARTITION="xjet"
# elif [[ "$MACHINE_ID" == "ursa" ]];then
#  WORK_DIR="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
#  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
#  QUEUE="${QUEUE:-batch}"
#  export HOMEreg=${HOMEreg}/${test_name}
  # export APRUN_REGRID=srun
  # PARTITION='u1-compute'
# elif [[ "$MACHINE_ID" == "orion" ]];then
#  WORK_DIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
#  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
#  QUEUE="${QUEUE:-batch}"
#  export HOMEreg=${HOMEreg}/${test_name}
  # export APRUN_REGRID=srun
  # PARTITION='orion'
  # ulimit -a
# elif [[ "$MACHINE_ID" == "hercules" ]];then
#  WORK_DIR="${WORK_DIR:-/work2/noaa/stmp/$LOGNAME}"
#  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
#  QUEUE="${QUEUE:-batch}"
#  export HOMEreg=${HOMEreg}/${test_name}
  # export APRUN_REGRID=srun
  # PARTITION='hercules'
# elif [[ "$MACHINE_ID" == "wcoss2" ]];then
#  WORK_DIR="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
#  PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
#  QUEUE="${QUEUE:-dev}"
#  export HOMEreg=${HOMEreg/$test_name}
  # PARTITION='dev'
  # export APRUN_REGRID="mpiexec -n 6 -ppn 6 --cpu-bind core"
# fi

DATA_DIR="${WORK_DIR}/reg-tests/${test_name}"
export NWPROD=${HOMEUFSUTILS}

# LOG_FILE=consistency.log01
export DATA="${DATA_DIR}/test1"

LOG_FILE=consistency.log
rm -f ${LOG_FILE}* summary.log

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 6 1 10G 0:05:00 hercules false gauss2fv3incr gauss2fv3incr.sh false
        ;;
    jet)
        submit_test 01 6 1 10G 0:05:00 xjet true gauss2fv3incr gauss2fv3incr.sh false
        ;;
    orion)
        submit_test 01 6 1 10G 0:05:00 orion false gauss2fv3incr gauss2fv3incr.sh false
        ;;
    ursa)
        submit_test 01 6 1 10G 0:05:00 u1-compute false gauss2fv3incr gauss2fv3incr.sh false
        ;;
    wcoss2)
        submit_test 01 6 1 10G 0:05:00 dev false gauss2fv3incr gauss2fv3incr.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

# if [[ "$MACHINE_ID" == "wcoss2" ]];then
#   TEST1=$(qsub -V -o "${LOG_FILE}" -e "${LOG_FILE}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=00:05:00 \
#         -N gauss2fv3incr -l select=1:ncpus=6:ompthreads=1:mem=10GB ./gauss2fv3incr.sh)
# else
#   TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A "${PROJECT_CODE}" -q "${QUEUE}" -J gauss2fv3incr \
#       -p "${PARTITION}" -o "${LOG_FILE}" -e "${LOG_FILE}" ./gauss2fv3incr.sh)
# fi



if [[ "${MACHINE_ID}" == "wcoss2" ]];then

  this_dir=${PWD}
  (qsub -V -o "${LOG_FILE}" -e "${LOG_FILE}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=00:01:00 \
        -N summary -l select=1:ncpus=1:mem=100MB -W "depend=afterany${TEST_IDS[*]}" << EOF
#!/bin/bash
cd ${PWD}
grep -a '<<<' ${LOG_FILE}?? | grep -v echo > ./summary.log
EOF
  ) &
  
else

  (sbatch --nodes=1  -t 0:01:00 -A "${PROJECT_CODE}" -J summary -o "${LOG_FILE}" -e "${LOG_FILE}" \
       -p "${PARTITION}" --open-mode=append -q "${QUEUE}" -d "afterany${TEST_IDS[*]}" << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}* > ./summary.log
EOF
  ) &
fi

sleep_time=0
echo "Waiting for ${test_name^^} testing to complete..."
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