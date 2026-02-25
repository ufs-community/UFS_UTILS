#!/bin/bash

set -ux

error() {
    set +x
    echo
    echo "$@" 1>&2
    exit 1
}

usage() {
    echo
    echo "Usage: $program [-c] [-m] [-h] [-b]"
    echo
    echo "  -b build the executable"
    echo "  -c create a new baseline"
    echo "  -m compare against the new baseline"
    echo "  -h display this help and exit"
    echo "  Examples"
    echo
    echo "    './rt.sh -b'  build exe file. compare against the existing baseline"
    echo "    './rt.sh -bc' build exe file. create a new baseline"
    echo "    './rt.sh -m'  do not build exe file. compare against the new baseline"
    echo
}

usage_and_exit() {
    set +x
    usage
    exit $1
}

readonly program=$(basename $0)
# PATHRT - Path to regression tests directory
# readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
# export PATHRT
# # PATHTR - Path to the UFS UTILS directory
# readonly PATHTR="$(cd $PATHRT/../.. && pwd)"
# export PATHTR



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

readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
export PATHRT
# PATHTR - Path to the UFS UTILS directory
readonly PATHTR="${HOMEUFSUTILS}"
export PATHTR

# source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1

# source $PATHTR/sorc/machine-setup.sh >/dev/null 2>&1
set +x
echo "Machine: ${MACHINE_ID}"
set -x

MOM6_version=20250128

# Adjust STMP, ACCOUNT and QUEUE as needed.
test_name=cpld_gridgen
STMP=${WORK_DIR:-?}
# export BASELINE_ROOT=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
export BASELINE_ROOT=${HOMEreg}/${test_name}/baseline_data
export NCCMP=nccmp
ACCOUNT=${PROJECT_CODE:?}

case ${MACHINE_ID} in
  orion)
    export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
    WLCLK=120
    PARTITION=orion
    ;;
  hercules)
    export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
    WLCLK=120
    PARTITION=hercules
    ;;
  ursa)
    export MOM6_FIXDIR=/scratch3/NCEPDEV/global/role.glopara/fix/mom6/${MOM6_version}
    WLCLK=40
    PARTITION=u1-compute
    ;;
  jet)
    export MOM6_FIXDIR=/lfs5/HFIP/hfv3gfs/glopara/FIX/fix/mom6/${MOM6_version}
    WLCLK=60
    PARTITION=xjet
    ;;
  wcoss2)
    export APRUN="mpiexec -n 12 -ppn 12 --cpu-bind core"
    export MOM6_FIXDIR=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/mom6/${MOM6_version}
    WLCLK=60
    PARTITION=dev
    ;;
  *)
    error "Unknown machine ${MACHINE_ID}"
    ;;

esac

# if [[ $MACHINE_ID = ursa ]]; then
#   STMP=${STMP:-/scratch4/NCEPDEV/stmp/$USER}
#                /scratch4/NCEPDEV/stmp/${LOGNAME}/ufs_utils
#   ACCOUNT=${ACCOUNT:-fv3-cpu}
#   QUEUE=${QUEUE:-batch}
#   WLCLK=40
#   export MOM6_FIXDIR=/scratch3/NCEPDEV/global/role.glopara/fix/mom6/${MOM6_version}
#   export NCCMP=nccmp
#   export BASELINE_ROOT=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
#                        /scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests
#   PARTITION=''
# elif [[ $MACHINE_ID = orion ]]; then
#   STMP=${STMP:-/work/noaa/stmp/$USER}
#   ACCOUNT=${ACCOUNT:-fv3-cpu}
#   QUEUE=${QUEUE:-batch}
#   WLCLK=120
#   export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
#   export NCCMP=nccmp
#   export BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/cpld_gridgen/baseline_data
#   PARTITION=''
#   ulimit -a
# elif [[ $MACHINE_ID = hercules ]]; then
#   STMP=${STMP:-/work2/noaa/stmp/$USER}
#   ACCOUNT=${ACCOUNT:-fv3-cpu}
#   QUEUE=${QUEUE:-batch}
#   WLCLK=120
#   export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
#   export BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/cpld_gridgen/baseline_data
#   export NCCMP=nccmp
#   PARTITION=''
#   ulimit -s unlimited
# elif [[ $MACHINE_ID = jet ]]; then
#   STMP=${STMP:-/lfs5/HFIP/emcda/$USER/stmp}
#   ACCOUNT=${ACCOUNT:-hfv3gfs}
#   QUEUE=${QUEUE:-batch}
#   WLCLK=60
#   export MOM6_FIXDIR=/lfs5/HFIP/hfv3gfs/glopara/FIX/fix/mom6/${MOM6_version}
#   export BASELINE_ROOT=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
#   export NCCMP=nccmp
#   PARTITION="--partition=xjet"
#   ulimit -s unlimited
# elif [[  $MACHINE_ID = wcoss2 ]]; then
#   STMP=${STMP:-/lfs/h2/emc/stmp/$USER}
#   ACCOUNT=${ACCOUNT:-GFS-DEV}
#   QUEUE=${QUEUE:-dev}
#   WLCLK=60
#   export MOM6_FIXDIR=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/mom6/${MOM6_version}
#   export BASELINE_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/cpld_gridgen/baseline_data
#   export APRUN="mpiexec -n 12 -ppn 12 --cpu-bind core"
#   export NCCMP=nccmp
# fi

NEW_BASELINE_ROOT=$STMP/${test_name}/baseline_data

BUILD_EXE=false
CREATE_BASELINE=false
while getopts :bcmh opt; do
    case $opt in
        b)
            BUILD_EXE=true
            ;;
        c)
            CREATE_BASELINE=true
            ;;
        m)
            BASELINE_ROOT=$NEW_BASELINE_ROOT
            ;;
        h)
            usage_and_exit 0
            ;;
        '?')
            error "$program: invalid option"
            ;;
    esac
done

export CREATE_BASELINE
if [[ $CREATE_BASELINE = true ]]; then
    rm -rf $NEW_BASELINE_ROOT
    mkdir -p $NEW_BASELINE_ROOT
fi

# compiler=${compiler:-intelllvm}
# if [[ "$compiler" == "intelllvm" ]]; then
#   if [[ ! -f ${PATHTR}/modulefiles/build.$MACHINE_ID.$compiler.lua ]];then
#     set +x
#     echo "IntelLLVM not available. Will use Intel Classic."
#     set -x
#     compiler=intel
#   fi
# fi
# export compiler
# set +x
# echo "Compiler: $compiler"
# set -x

# Build the executable file
if [[ $BUILD_EXE = true ]]; then
    COMPILE_LOG=compile.log
    cd "$PATHTR" || exit 1; echo "Can't CD into '${PATHTR}'"
    rm -rf $COMPILE_LOG $PATHTR/build $PATHTR/exec $PATHTR/lib
    ./build_all.sh >$PATHRT/$COMPILE_LOG 2>&1 && d=$? || d=$?
    if [[ d -ne 0 ]]; then
        error "Build did not finish successfully. Check $COMPILE_LOG"
    else
        set +x
        echo "Build was successful"
        set -x
        cd $PATHRT
    fi
else
    if [[ ! -f $PATHTR/exec/cpld_gridgen ]]; then
       error "cpld_gridgen exe file is not found in $PATHTR/exe/. Try -b to build or -h for help."
    fi
fi

# module use $PATHTR/modulefiles
# module use ${HOMEUFSUTILS}/modulefiles
# module load build.${MACHINE_ID,,}.$compiler
if [[ $MACHINE_ID = wcoss2 ]]; then
  module load nccmp-D/1.9.0.1
fi
set +x
module list
set -x

RUNDIR_ROOT=$STMP/reg-tests/${test_name}/rt_$$

declare -A tests
all_tests=""

rm -f fail_test* nccmp_*.log summary.log run_*log RegressionTests_${MACHINE_ID,,}.$compiler.*.log

# Kick off all tests.

i=0
while read -r line || [ "$line" ]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line =~ \# ]] && continue

  TEST_NAME=$(echo $line | cut -d'|' -f1 | sed -e 's/^ *//' -e 's/ *$//')
  TEST_NAME=${TEST_NAME##mx}
  ATMLIST=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
  if [[ -z ${ATMLIST} ]]; then
      ATMLIST=-1
  fi

  if [[ ${TEST_NAME} == "008" ]]; then
      NTASKS=24
  else
      NTASKS=12
  fi

  export NEW_BASELINE=${NEW_BASELINE_ROOT}/$TEST_NAME
  RUNDIR=$RUNDIR_ROOT/$TEST_NAME
  mkdir -p $RUNDIR
  export RUNDIR
  export OUTDIR_PATH=$RUNDIR
  export BASELINE=$BASELINE_ROOT/$TEST_NAME
  export REGRESSIONTEST_LOG=RegressionTests_$MACHINE_ID.$compiler.${TEST_NAME}.log

  cp $PATHRT/parm/grid.nml.IN $RUNDIR
  cp $PATHTR/exec/cpld_gridgen $RUNDIR

  if [[ $MACHINE_ID = wcoss2 ]]; then
    tests[$i]=$(qsub -V -o $PATHRT/run_${TEST_NAME}.log -e $PATHRT/run_${TEST_NAME}.log -q $QUEUE  -A $ACCOUNT \
       -l walltime=00:${WLCLK}:00 -N $TEST_NAME -l select=1:ncpus=${NTASKS} -v RESNAME=$TEST_NAME,ATMLIST="'$ATMLIST'" ./cpld_gridgen.sh)

  else
    tests[$i]=$(sbatch --parsable --ntasks-per-node=${NTASKS} --nodes=1 -t 00:${WLCLK}:00 -A $ACCOUNT -q $QUEUE -J $TEST_NAME \
            --partition=$PARTITION -o run_${TEST_NAME}.log -e run_${TEST_NAME}.log ./cpld_gridgen.sh "$TEST_NAME" "$ATMLIST")
  fi

  all_tests=${all_tests}":"${tests[$i]}

  ((i=i+1))

done < ./rt.conf

# export target

# Once all the jobs are finished, this summary job will run.

if [[ $MACHINE_ID = wcoss2 ]]; then

  (qsub -V -o /dev/null -e /dev/null -q $QUEUE -A $ACCOUNT -l walltime=00:01:00 \
        -N summary -l select=1:ncpus=1:mem=100MB \
        -W depend=afterany${all_tests} ./rt.summary.sh) &
else

  (sbatch --nodes=1 -t 0:01:00 -A $ACCOUNT -J summary -o /dev/null -e /dev/null \
       --partition=$PARTITION --open-mode=append -q $QUEUE -d afterany${all_tests} ./rt.summary.sh) &

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
