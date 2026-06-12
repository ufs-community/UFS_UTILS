#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on Hercules.
#
# Set values in ../rt.control to specify the number of tasks, memory, and walltime
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
    local slurmcluster="$1"; shift
    local exclusive="$1"; shift
    local jobname="$1"; shift
    local script="$1"; shift
    local waitonjobid="$1"; shift

    local logfile="${LOG_FILE}${suffix}"
    export OMP_NUM_THREADS=1  # should match cpus-per-task

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
        export APRUN="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core"
        jobid=$(qsub -V -o "${logfile}" -e "${logfile}" -q "${QUEUE}" -A "${PROJECT_CODE}" -l walltime=${walltime} \
                -N "${jobname}" -l select=${nodes}:ncpus=${ntasks_per_node}:ompthreads=${OMP_NUM_THREADS}:mem=${mem} \
                ${dep_flag_pbs:+"${dep_flag_pbs}"} ./${script})
        # jobid=${jobid%.*}
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        export APRUN="srun"
        jobid=$(sbatch --parsable --partition="${partition}" ${slurmflag:+"${slurmflag}"} --ntasks-per-node="${ntasks_per_node}" --nodes="${nodes}" --mem="${mem}" -t "${walltime}" \
                -A "${PROJECT_CODE}" -q "${QUEUE}" -J "${jobname}" --open-mode=append ${exclusive_flag:+"${exclusive_flag}"} \
                ${dep_flag_slurm:+"${dep_flag_slurm}"} -o "${logfile}" -e "${logfile}" "./${script}")
        # jobid=${jobid%.*}
        # jobid=${jobid%%;*}
    else
        echo "Error: Unsupported scheduler '${SCHEDULER}'"
        exit 1
    fi
    status=$?
    if [ $status -ne 0 ]; then
        echo "Error submitting job: $output"
        exit 1
    fi
    jobid=${jobid%.*}
    jobid=${jobid%%;*}
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

test_name="chgres_cube"
export OUTDIR="${WORK_DIR}/reg-tests/${test_name}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------
export HOMEreg=${HOMEreg}/${test_name}

BASELINE_ROOT=${HOMEreg}/baseline_data
WEIGHTS_ROOT=${HOMEreg}/../cpld_gridgen/baseline_data
INPUT_ROOT=${HOMEreg}/input_data
ACCOUNT=${PROJECT_CODE}

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    source ../get_hash.sh
fi

HOMEufs=$PWD/../..
export HOMEufs
this_dir=$PWD

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $SUM_FILE ${LOG_FILE}*

export OMP_STACKSIZE=1024M

export NCCMP=${NCCMP:-nccmp}

rm -fr "$OUTDIR"

declare -a TEST_IDS=()

case ${MACHINE_ID,,} in
    hercules)
        submit_test 01 6 1 75G 0:15:00 hercules false false false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 hercules false false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 6 2 75G 0:10:00 hercules false false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 hercules false false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 12 1 75G 0:10:00 hercules false false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 hercules false false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 75G 0:10:00 hercules false false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 12 1 75G 0:10:00 hercules false false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 12 1 75G 0:10:00 hercules false false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 12 1 75G 0:10:00 hercules false false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 hercules false false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 12 1 75G 0:10:00 hercules false false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 hercules false false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 12 1 75G 0:10:00 hercules false false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    orion)
        submit_test 01 6 1 75G 0:15:00 orion false false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 orion false false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 6 2 75G 0:10:00 orion false false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 orion false false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 12 1 75G 0:10:00 orion false false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 orion false false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 75G 0:10:00 orion false false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 12 1 75G 0:10:00 orion false false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 12 1 75G 0:10:00 orion false false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 12 1 75G 0:10:00 orion false false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 orion false false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 12 1 75G 0:10:00 orion false false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 orion false false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 12 1 75G 0:10:00 orion false false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    ursa)
        submit_test 01 6 1 50G 0:15:00 u1-compute false false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 2 100G 0:15:00 u1-compute false false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 12 1 100G 0:15:00 u1-compute false false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 50G 0:05:00 u1-compute false false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 6 1 50G 0:05:00 u1-compute false false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 100G 0:10:00 u1-compute false false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 100G 0:10:00 u1-compute false false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 6 1 50G 0:05:00 u1-compute false false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 6 1 100G 0:05:00 u1-compute false false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 6 1 100G 0:05:00 u1-compute false false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:15:00 u1-compute false false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 6 1 100G 0:05:00 u1-compute false false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 50G 0:05:00 u1-compute false false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 6 1 100G 0:05:00 u1-compute false false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    gaeac6)
        submit_test 01 6 1 0 0:15:00 batch c6 false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 2 0 0:15:00 batch c6 false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 12 1 0 0:15:00 batch c6 false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 0 0:05:00 batch c6 false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 6 1 0 0:05:00 batch c6 false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 0 0:10:00 batch c6 false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 6 2 0 0:10:00 batch c6 false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 6 1 0 0:05:00 batch c6 false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 6 1 0 0:05:00 batch c6 false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 6 1 0 0:05:00 batch c6 false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 0 0:15:00 batch c6 false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 6 1 0 0:05:00 batch c6 false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 0 0:05:00 batch c6 false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 6 1 0 0:05:00 batch c6 false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
        ;;
    wcoss2)
        submit_test 01 6 1 75G 0:15:00 dev false false c96.fv3.restart c96.fv3.restart.sh false
        submit_test 02 6 1 75G 0:15:00 dev false false c192.fv3.history c192.fv3.history.sh false
        submit_test 03 12 1 75G 0:10:00 dev false false c96.fv3.netcdf c96.fv3.netcdf.sh false
        submit_test 04 6 1 75G 0:05:00 dev false false c192.gfs.grib2 c192.gfs.grib2.sh false
        submit_test 05 6 1 75G 0:10:00 dev false false 25km.conus.gfs.grib2 25km.conus.gfs.grib2.sh false
        submit_test 06 6 1 75G 0:10:00 dev false false 3km.conus.hrrr.gfssdf.grib2 3km.conus.hrrr.gfssdf.grib2.sh false
        submit_test 07 12 1 75G 0:10:00 dev false false 3km.conus.hrrr.newsfc.grib2 3km.conus.hrrr.newsfc.grib2.sh false
        submit_test 08 6 1 75G 0:10:00 dev false false 13km.conus.nam.grib2 13km.conus.nam.grib2.sh false
        submit_test 09 6 1 75G 0:10:00 dev false false 13km.conus.rap.grib2 13km.conus.rap.grib2.sh false
        submit_test 10 6 1 75G 0:10:00 dev false false 13km.na.gfs.ncei.grib2 13km.na.gfs.ncei.grib2.sh false
        submit_test 11 12 1 100G 0:25:00 dev false false c96.fv3.netcdf2wam c96.fv3.netcdf2wam.sh false
        submit_test 12 6 1 75G 0:10:00 dev false false 25km.conus.gfs.pbgrib2 25km.conus.gfs.pbgrib2.sh false
        submit_test 13 6 1 75G 0:05:00 dev false false c96.gefs.grib2 c96.gefs.grib2.sh false
        submit_test 14 6 1 75G 0:10:00 dev false false 13km.conus.rap-smoke.grib2 13km.conus.rap-smoke.grib2.sh false
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
    (qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N chgres_summary -l select=1:ncpus=1:mem=100MB \
        -W depend="afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd ${this_dir}
grep -a '^<<<' ${LOG_FILE}* | grep -v echo > ${SUM_FILE}
EOF
) &
elif [[ "${SCHEDULER}" == "slurm" ]]; then
    (sbatch --nodes=1 -t 0:01:00 -A "${PROJECT_CODE}" ${slurmflag:+"${slurmflag}"} -J chgres_summary -o "${LOG_FILE}" -e "${LOG_FILE}" \
       --open-mode=append -q "${QUEUE}" \
       -d "afterany$(echo "${TEST_IDS[*]}" | tr -d '[:space:]')" << EOF
#!/bin/bash
cd ${this_dir}
grep -a '^<<<' ${LOG_FILE}*  > ${SUM_FILE}
EOF
) &
else
    echo "Error: Unsupported scheduler '${SCHEDULER}'"
    exit 1
fi

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
