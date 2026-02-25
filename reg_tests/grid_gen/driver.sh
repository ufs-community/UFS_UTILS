#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on Ursa.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.  To see which projects you are authorized to use,
# type "account_params".
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
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUNCY="srun"
        export APRUN_SFC=${APRUNCY}
        jobid=$(sbatch --parsable --partition="${partition}" --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
               -A "${PROJECT_CODE}" -q "${QUEUE}" -J "${jobname}" --open-mode=append ${exclusive_flag:+"${exclusive_flag}"} \
               --export=ALL,OMP_NUM_THREADS=${OMP_NUM_THREADS} ${dep_flag_slurm:+"${dep_flag_slurm}"} \
               -o "${logfile}" -e "${logfile}" "./${script}")
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

# source ../../sorc/machine-setup.sh > /dev/null 2>&1
# module use ../../modulefiles
# module load build.${MACHINE_ID}.$compiler
# module list

set -x

# export WORK_DIR="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
test_name="grid_gen"
export WORK_DIR="${WORK_DIR}/reg-tests/${test_name}"

# QUEUE="${QUEUE:-batch}"
# PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

# UPDATE_BASELINE="${UPDATE_BASELINE:-FALSE}"
# export UPDATE_BASELINE

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

LOG_FILE=consistency.log
SUM_FILE=summary.log
export home_dir=$PWD/../..
export APRUN=time
export this_dir=$PWD
export OMP_STACKSIZE=2048m
# export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/grid_gen
HOMEreg="${HOMEreg}/${test_name}"

ulimit -a

declare -a TEST_IDS=()

rm -f ${LOG_FILE}* ${SUM_FILE}

rm -fr "${WORK_DIR}"

# export OMP_NUM_THREADS=24

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
    jet)
        submit_test 01 24 1 50G 0:20:00 xjet true c96.uniform c96.uniform.sh false
        submit_test 02 12 4 300G 0:15:00 xjet true c96.viirs.bnu c96.viirs.bnu.sh false
        submit_test 03 24 1 50G 0:07:00 xjet true gfdl.regional gfdl.regional.sh false
        submit_test 04 24 1 50G 0:07:00 xjet true esg.regional esg.regional.sh false
        submit_test 05 24 1 50G 0:07:00 xjet true esg.regional.pct.cat esg.regional.pct.cat.sh false
        submit_test 06 12 1 50G 0:07:00 xjet true reg.gsl.gwd.12 regional.gsl.gwd.sh false
        submit_test 07 24 1 50G 0:07:00 xjet true reg.gsl.gwd.24 regional.gsl.gwd.sh false
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

# #-----------------------------------------------------------------------------
# # C96 uniform grid
# #-----------------------------------------------------------------------------

# LOG_FILE1=${LOG_FILE}01
# TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
#       -o $LOG_FILE1 -e $LOG_FILE1 ./c96.uniform.sh)

# #-----------------------------------------------------------------------------
# # C96 uniform grid using viirs vegetation and bnu soil data.
# #-----------------------------------------------------------------------------

# LOG_FILE2=${LOG_FILE}02
# TEST2=$(sbatch --parsable --ntasks-per-node=12 --nodes=2 --mem=300g -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.viirs.bnu \
#       -o $LOG_FILE2 -e $LOG_FILE2 ./c96.viirs.bnu.sh)

# #-----------------------------------------------------------------------------
# # gfdl regional grid
# #-----------------------------------------------------------------------------

# LOG_FILE3=${LOG_FILE}03
# TEST3=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J gfdl.regional \
#       -o $LOG_FILE3 -e $LOG_FILE3 ./gfdl.regional.sh)

# #-----------------------------------------------------------------------------
# # ESG regional grid (output dominant soil/vegetation type).
# #-----------------------------------------------------------------------------

# LOG_FILE4=${LOG_FILE}04
# TEST4=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional \
#       -o $LOG_FILE4 -e $LOG_FILE4 ./esg.regional.sh)

# #-----------------------------------------------------------------------------
# # ESG regional grid (output percent of each soil and vegetation type and
# # the dominant category).
# #-----------------------------------------------------------------------------

# LOG_FILE5=${LOG_FILE}05
# TEST5=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J esg.regional.pct.cat \
#       -o $LOG_FILE5 -e $LOG_FILE5 ./esg.regional.pct.cat.sh)

# #-----------------------------------------------------------------------------
# # Regional GSL gravity wave drag test. This test is run with varying
# # thread counts.
# #-----------------------------------------------------------------------------

# export nthreads=12
# LOG_FILE6=${LOG_FILE}06
# TEST6=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J reg.gsl.gwd.12 \
#       -o $LOG_FILE6 -e $LOG_FILE6 ./regional.gsl.gwd.sh)

# export nthreads=24
# LOG_FILE7=${LOG_FILE}07
# TEST7=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J reg.gsl.gwd.24 \
#       -o $LOG_FILE7 -e $LOG_FILE7 ./regional.gsl.gwd.sh)

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
grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
EOF
  ) &
else
    echo "Error: Unsupported scheduler '${SCHEDULER}'"
    exit 1
fi

# if [[ "${waitlocal}" == "true" ]]; then
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
