#!/bin/bash

submit_test() {
    # the only echo should be the jobid of the last submitted job.
    # The caller can capture this and use it for dependencies.

    # local suffix="$1"; shift
    # local ntasks_per_node="$1"; shift
    # local nodes="$1"; shift
    # local mem="$1"; shift
    # local walltime="$1"; shift
    # local partition="$1"; shift
    # local slurmcluster="$1"; shift
    # local exclusive="$1"; shift
    # local jobname="$1"; shift
    # local script="$1"; shift
    # local waitonjobid="$1"; shift
    if [[ $1 == "end" ]]; then
        END_SUBMIT="true"
        local logsuffix=""; shift
    else
        END_SUBMIT="false"
        local logsuffix="$1"; shift
    fi

    local ntasks_per_node="$1"; shift
    local nodes="$1"; shift
    local mem="$1"; shift
    local walltime="$1"; shift
    #local partition="$1"; shift #LETS PROVIDE THIS IN THE RT.CONTROL FILE
    #local slurmcluster="$1"; shift#
    local exclusive="$1"; shift
    local jobname="$1"; shift
    local script="$1"; shift
    local remainingargs=("$@");
    local waitonjobid=$(printf '%s' "${remainingargs[@]}")

    LOG_FILE=${LOG_FILE:-"slurm_job"}
    local partition=${PARTITION:?}
    local slurmcluster=${SLURM_CLUSTER:-false}

    local logfile="${LOG_FILE}${logsuffix}"

    export OMP_NUM_THREADS=${ntasks_per_node}  # should match cpus-per-task
    export APRUN_ADDITIONS=${APRUN_ADDITIONS:-""}
    export DATA_ROOT=${DATA_ROOT:-""}
    export DATA="${DATA:-${DATA_ROOT}/test${logsuffix}}"

    if [[ "${SCHEDULER}" == "pbs" ]]; then
        export APRUN="mpiexec -n ${ntasks_per_node} -ppn ${ntasks_per_node} --cpu-bind core ${APRUN_ADDITIONS}"
        export APRUNCY="${APRUN} --depth ${OMP_NUM_THREADS}"
        export APRUN_SFC=${APRUN_SFC:-${APRUNCY}}
        export APRUN_REGRID=${APRUN_REGRID:-${APRUN}}
        local pbs_args=(-V)
        if [[ -n "${logfile}" && "${logfile}" != "false" ]]; then
            pbs_args+=(-o "${logfile}" -e "${logfile}")
        fi
        if [[ -n "${QUEUE}" && "${QUEUE}" != "false" ]]; then
            pbs_args+=(-q "${QUEUE}")
        fi
        if [[ -n "${PROJECT_CODE}" && "${PROJECT_CODE}" != "false" ]]; then
            pbs_args+=(-A "${PROJECT_CODE}")
        fi
        if [[ -n "${walltime}" && "${walltime}" != "false" ]]; then
            pbs_args+=(-l "walltime=${walltime}")
        fi
        if [[ -n "${jobname}" && "${jobname}" != "false" ]]; then
            pbs_args+=(-N "${jobname}")
        fi
        if [[ ( -n "${nodes}" && "${nodes}" != "false" ) || ( -n "${ntasks_per_node}" && "${ntasks_per_node}" != "false" ) || ( -n "${mem}" && "${mem}" != "false" ) ]]; then
            local select_spec="select"
            if [[ -n "${nodes}" && "${nodes}" != "false" ]]; then
                select_spec+="=${nodes}"
            fi
            if [[ -n "${ntasks_per_node}" && "${ntasks_per_node}" != "false" ]]; then
                if [[ "${select_spec}" == "select" ]]; then
                    select_spec+="=ncpus=${ntasks_per_node}"
                else
                    select_spec+=":ncpus=${ntasks_per_node}"
                fi
            fi
            select_spec+=":ompthreads=${OMP_NUM_THREADS}"
            if [[ -n "${mem}" && "${mem}" != "false" ]]; then
                if [[ "${select_spec}" == "select" ]]; then
                    select_spec+="=mem=${mem}"
                else
                    select_spec+=":mem=${mem}"
                fi
            fi
            pbs_args+=(-l "${select_spec}")
        fi
        if [[ -n "${waitonjobid}" && "${waitonjobid}" != "false" ]]; then
            pbs_args+=(-W "depend=afterok${waitonjobid}")
        fi
        if [[ ${END_SUBMIT} == "true" ]]; then
            pbs_args+=(-W "block=true")
        fi
        jobid=$(qsub "${pbs_args[@]}" "./${script}")
    elif [[ "${SCHEDULER}" == "slurm" ]]; then
        # export APRUN="srun --mpi=pmi2"
        export APRUN=${APRUN:-"srun ${APRUN_ADDITIONS}"}
        export APRUNCY=${APRUNCY:-"srun ${APRUN_ADDITIONS}"}
        export APRUN_SFC=${APRUN_SFC:-"srun ${APRUN_ADDITIONS}"}
        export APRUN_REGRID=${APRUN_REGRID:-"srun ${APRUN_ADDITIONS}"}
        # SLURM Items
        local slurm_args=(--parsable)
        if [[ -n "${partition}" && "${partition}" != "false" ]]; then
            slurm_args+=(--partition="${partition}")
        fi
        if [[ -n "${slurmcluster}" && "${slurmcluster}" != "false" ]]; then
            slurm_args+=(--clusters="${slurmcluster}")
        fi
        if [[ -n "${ntasks_per_node}" && "${ntasks_per_node}" != "false" ]]; then
            slurm_args+=(--ntasks-per-node="${ntasks_per_node}")
        fi
        if [[ -n "${nodes}" && "${nodes}" != "false" ]]; then
            slurm_args+=(--nodes="${nodes}")
        fi
        if [[ -n "${mem}" && "${mem}" != "false" ]]; then
            slurm_args+=(--mem="${mem}")
        fi
        if [[ -n "${walltime}" && "${walltime}" != "false" ]]; then
            slurm_args+=(-t "${walltime}")
        fi
        if [[ -n "${jobname}" && "${jobname}" != "false" ]]; then
            slurm_args+=(-J "${jobname}")
        fi
        slurm_args+=(--open-mode=append)
        if [[ -n "${exclusive}" && "${exclusive}" != "false" ]]; then
            slurm_args+=(--exclusive)
        fi
        if [[ -n "${waitonjobid}" && "${waitonjobid}" != "false" ]]; then
            slurm_args+=(--dependency="afterok${waitonjobid}")
        fi
        if [[ -n "${logfile}" && "${logfile}" != "false" ]]; then
            slurm_args+=(-o "${logfile}" -e "${logfile}")
        fi
        jobid=$(sbatch "${slurm_args[@]}" "./${script}")
    else
        echo "Error: Unsupported scheduler '${SCHEDULER}'"
        exit 1
    fi
    status=$?
    if [ $status -ne 0 ]; then
        echo "Error submitting job: ${jobid}"
        exit 1
    fi
 
    if [[ ${END_SUBMIT} != "true" ]]; then
        jobid=${jobid%.*}
        jobid=${jobid%%;*}
        if [[ "${jobid}" == "" ]]; then
            echo "Error submitting job to slurm scheduler"
            exit 1
        fi
        TEST_IDS+=(":${jobid}")
        echo ":${jobid}"
    fi
    export TEST_IDS
}