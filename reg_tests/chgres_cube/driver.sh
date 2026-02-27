#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on Hercules.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube.
# Set the PROJECT_CODE and QUEUE as appropriate.  To see which projects you 
# are authorized to use, type:
#
#   $ sacctmgr show associations where user-$USER format=account%20,qos%50.
#
# Invoke the script with no arguments.  A series of daily-chained
# consistency tests will be submitted.  To check the queue, type:
# "squeue -u $LOGNAME".
#
# The run output will be stored in OUTDIR.  Standard output from 
# each test will be placed in its own log file. Once the suite
# has completed, a summary of results is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg.
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
    export OMP_NUM_THREADS=1  # should match cpus-per-task

    if [[ "${exclusive}" == "true" ]]; then
        exclusive_flag="--exclusive"
    fi

    if [[ "${waitonjobid}" != "false" ]]; then
        dep_flag_slurm="--dependency=afterok:${waitonjobid}"
        dep_flag_pbs="-W depend=afterok:${waitonjobid}"
    fi

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUN="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=${OMP_NUM_THREADS}:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} ./${script})
        jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUN="srun"
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
    # shellcheck source=/dev/null
    source "${RT_DIR}/rt.control"
else
    echo "ERROR: Cannot find rt.control script"
    exit 1
fi

# shellcheck source=${HOMEUFSUTILS}/sorc/machine-setup.sh
# source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1
# module use "${HOMEUFSUTILS}/modulefiles"

# # source "${RT_DIR}/rt.control"
# # source ../../sorc/machine-setup.sh > /dev/null 2>&1
# # module use ../../modulefiles
# module load build.${MACHINE_ID,,}.$compiler
# module list

# if [ "${MACHINE_ID,,}" == "hercules" ]; then
#     ulimit -s unlimited
# fi

test_name="chgres_cube"
export OUTDIR="${WORK_DIR}/reg-tests/${test_name}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

#export HOMEreg=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/chgres_cube
export HOMEreg=${HOMEreg}/${test_name}


# BASELINE_ROOT=${HOMEreg}/${test_name}/baseline_data
BASELINE_ROOT=${HOMEreg}/baseline_data
# WEIGHTS_ROOT=${HOMEreg}/cpld_gridgen/baseline_data
WEIGHTS_ROOT=${HOMEreg}/../cpld_gridgen/baseline_data
# INPUT_ROOT=${HOMEreg}/${test_name}/input_data
INPUT_ROOT=${HOMEreg}/input_data
# STMP=${WORKDIR}
ACCOUNT=${PROJECT_CODE}

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    source ../get_hash.sh
fi

HOMEufs=$PWD/../..
export HOMEufs

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $SUM_FILE ${LOG_FILE}*

export OMP_STACKSIZE=1024M

# export APRUN=srun

# export machine=${MACHINE_ID,,}

export NCCMP=${NCCMP:-nccmp}

rm -fr "$OUTDIR"

declare -a TEST_IDS=()

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 6 1 75G 0:15:00 hercules false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 hercules false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 6 2 75G 0:10:00 hercules false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 hercules false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 12 1 75G 0:10:00 hercules false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 hercules false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 75G 0:10:00 hercules false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 12 1 75G 0:10:00 hercules false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 12 1 75G 0:10:00 hercules false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 12 1 75G 0:10:00 hercules false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 hercules false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 12 1 75G 0:10:00 hercules false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 hercules false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 12 1 75G 0:10:00 hercules false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    jet)
        submit_test 01 6 1 75G 0:15:00 xjet true c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 xjet true c192.fv3.history c192.fv3.history.sh false
        submit_test 03 6 2 75G 0:10:00 xjet true c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 xjet true c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 12 1 75G 0:10:00 xjet true 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 xjet true 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 75G 0:10:00 xjet true 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 12 1 75G 0:10:00 xjet true 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 12 1 75G 0:10:00 xjet true 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 12 1 75G 0:10:00 xjet true 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 xjet true c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 12 1 75G 0:10:00 xjet true 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 xjet true c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 12 1 75G 0:10:00 xjet true 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
     ;;
    orion)
        submit_test 01 6 1 75G 0:15:00 orion false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 orion false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 6 2 75G 0:10:00 orion false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 orion false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 12 1 75G 0:10:00 orion false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 orion false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 75G 0:10:00 orion false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 12 1 75G 0:10:00 orion false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 12 1 75G 0:10:00 orion false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 12 1 75G 0:10:00 orion false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 orion false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 12 1 75G 0:10:00 orion false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 orion false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 12 1 75G 0:10:00 orion false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    ursa)
        submit_test 01 6 1 50G 0:15:00 u1-compute false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 2 100G 0:15:00 u1-compute false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 12 1 100G 0:15:00 u1-compute false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 50G 0:05:00 u1-compute false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 6 1 50G 0:05:00 u1-compute false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 100G 0:10:00 u1-compute false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 100G 0:10:00 u1-compute false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 6 1 50G 0:05:00 u1-compute false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 6 1 100G 0:05:00 u1-compute false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 6 1 100G 0:05:00 u1-compute false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 u1-compute false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 6 1 100G 0:05:00 u1-compute false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 50G 0:05:00 u1-compute false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 6 1 100G 0:05:00 u1-compute false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    wcoss2)
        submit_test 01 6 1 75G 0:15:00 dev false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 dev false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 12 1 75G 0:10:00 dev false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 dev false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 6 1 75G 0:10:00 dev false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 dev false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 12 1 75G 0:10:00 dev false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 6 1 75G 0:10:00 dev false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 6 1 75G 0:10:00 dev false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 6 1 75G 0:10:00 dev false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:25:00 dev false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 6 1 75G 0:10:00 dev false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 dev false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 6 1 75G 0:10:00 dev false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    *)
        echo "Error: Unsupported machine '${MACHINE_ID}'"
        exit 1
        ;;
esac

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

# sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
#        --open-mode=append -q $QUEUE \
#        -d afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6:$TEST7:$TEST8:$TEST9:$TEST10:$TEST11:$TEST12:$TEST13:$TEST14 << EOF
# #!/bin/bash
# grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
# EOF
if [[ "${SCHEDULER}" == "pbs" ]]; then
    (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N chgres_summary -l select=1:ncpus=1:mem=100MB \
        -W depend="afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
grep -a '^<<<' ${LOG_FILE}* | grep -v echo > ${SUM_FILE}
EOF
) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
    (sbatch --nodes=1 -t 0:01:00 -A "${PROJECT_CODE}" -J chgres_summary -o "${LOG_FILE}" -e "${LOG_FILE}" \
       --open-mode=append -q "${QUEUE}" \
       -d "afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
grep -a '^<<<' ${LOG_FILE}*  > ${SUM_FILE}
EOF
) &
else
    echo "Error: Unsupported scheduler '${SCHEDULER}'"
    exit 1
fi

# (sbatch --nodes=1 -t 0:01:00 -A "${PROJECT_CODE}" -J chgres_summary -o "${LOG_FILE}" -e "${LOG_FILE}" \
#        --open-mode=append -q "${QUEUE}" \
#        -d "afterok$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
# #!/bin/bash
# grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
# EOF
# ) &
echo "Waiting for summary log to get generated..."
TIMEOUT_LIMIT=${TIMEOUT_LIMIT:?}  # default to 1 hour
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
