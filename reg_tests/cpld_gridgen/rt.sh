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
    echo
    echo "  -c create a new baseline"
    echo
    echo "  -m compare against the new baseline"
    echo
    echo "  -h display this help and exit"
    echo
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
readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
export PATHRT
# PATHTR - Path to the UFS UTILS directory
readonly PATHTR="$(cd $PATHRT/../.. && pwd)"
export PATHTR

source $PATHTR/sorc/machine-setup.sh >/dev/null 2>&1
set +x
echo "Machine: $target"
set -x

WLCLK=20
MOM6_version=20250128

# Adjust STMP, ACCOUNT and QUEUE as needed.

if [[ $target = hera ]]; then
  STMP=${STMP:-/scratch2/NCEPDEV/stmp1/$USER}
  ACCOUNT=${ACCOUNT:-fv3-cpu}
  QUEUE=${QUEUE:-batch}
  export MOM6_FIXDIR=/scratch1/NCEPDEV/global/glopara/fix/mom6/${MOM6_version}
  export NCCMP=nccmp
  BASELINE_ROOT=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
  PARTITION=hera
elif [[ $target = orion ]]; then
  STMP=${STMP:-/work/noaa/stmp/$USER}
  ACCOUNT=${ACCOUNT:-fv3-cpu}
  QUEUE=${QUEUE:-batch}
  export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
  export NCCMP=nccmp
  BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/cpld_gridgen/baseline_data
  PARTITION=orion
  ulimit -s unlimited
elif [[ $target = hercules ]]; then
  STMP=${STMP:-/work2/noaa/stmp/$USER}
  ACCOUNT=${ACCOUNT:-fv3-cpu}
  QUEUE=${QUEUE:-batch}
  export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/${MOM6_version}
  BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/cpld_gridgen/baseline_data
  export NCCMP=nccmp
  PARTITION=hercules
  ulimit -s unlimited
elif [[ $target = jet ]]; then
  STMP=${STMP:-/lfs5/HFIP/emcda/$USER/stmp}
  ACCOUNT=${ACCOUNT:-hfv3gfs}
  QUEUE=${QUEUE:-batch}
  export MOM6_FIXDIR=/lfs5/HFIP/hfv3gfs/glopara/FIX/fix/mom6/${MOM6_version}
  BASELINE_ROOT=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
  export NCCMP=nccmp
  PARTITION=xjet
  ulimit -s unlimited
elif [[  $target = wcoss2 ]]; then
  STMP=${STMP:-/lfs/h2/emc/stmp/$USER}
  ACCOUNT=${ACCOUNT:-GFS-DEV}
  QUEUE=${QUEUE:-dev}
  export MOM6_FIXDIR=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/mom6/${MOM6_version}
  BASELINE_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/cpld_gridgen/baseline_data
  export APRUN="mpiexec -n 1 -ppn 1 --cpu-bind core"
  export NCCMP=nccmp
fi

NEW_BASELINE_ROOT=$STMP/CPLD_GRIDGEN/BASELINE

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

compiler=${compiler:-intelllvm}
if [[ "$compiler" == "intelllvm" ]]; then
  if [[ ! -f ${PATHTR}/modulefiles/build.$target.$compiler.lua ]];then
    set +x
    echo "IntelLLVM not available. Will use Intel Classic."
    set -x
    compiler=intel
  fi
fi
export compiler
set +x
echo "Compiler: $compiler"
set -x

# Build the executable file
if [[ $BUILD_EXE = true ]]; then
    COMPILE_LOG=compile.log
    cd $PATHTR
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

module use $PATHTR/modulefiles
module load build.$target.$compiler
if [[ $target = wcoss2 ]]; then
  module load netcdf
  module load nccmp
fi
set +x
module list
set -x

RUNDIR_ROOT=$STMP/CPLD_GRIDGEN/

declare -A tests
all_tests=""

rm -f fail_test* nccmp_*.log summary.log run_*log RegressionTests_$target.$compiler.*.log

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

  export NEW_BASELINE=${NEW_BASELINE_ROOT}/$TEST_NAME
  RUNDIR=$RUNDIR_ROOT/$TEST_NAME
  rm -fr $RUNDIR
  mkdir -p $RUNDIR
  export RUNDIR
  export OUTDIR_PATH=$RUNDIR
  export BASELINE=$BASELINE_ROOT/$TEST_NAME
  export REGRESSIONTEST_LOG=RegressionTests_$target.$compiler.${TEST_NAME}.log

  cp $PATHRT/parm/grid.nml.IN $RUNDIR
  cp $PATHTR/exec/cpld_gridgen $RUNDIR
  
  if [[ $target = wcoss2 ]]; then
    tests[$i]=$(qsub -V -o $PATHRT/run_${TEST_NAME}.log -e $PATHRT/run_${TEST_NAME}.log -q $QUEUE  -A $ACCOUNT \
       -l walltime=00:${WLCLK}:00 -N $TEST_NAME -l select=1:ncpus=1:mem=12GB -v RESNAME=$TEST_NAME,ATMLIST="'$ATMLIST'" ./cpld_gridgen.sh)

  else
    tests[$i]=$(sbatch --parsable --ntasks-per-node=1 --nodes=1 -t 00:${WLCLK}:00 -A $ACCOUNT -q $QUEUE -J $TEST_NAME \
            --partition=$PARTITION -o run_${TEST_NAME}.log -e run_${TEST_NAME}.log ./cpld_gridgen.sh "$TEST_NAME" "$ATMLIST")
  fi

  all_tests=${all_tests}":"${tests[$i]}

  ((i=i+1))

done < ./rt.conf

export target

# Once all the jobs are finished, this summary job will run.

if [[ $target = wcoss2 ]]; then

  qsub -V -o /dev/null -e /dev/null -q $QUEUE -A $ACCOUNT -l walltime=00:01:00 \
        -N summary -l select=1:ncpus=1:mem=100MB \
        -W depend=afterok${all_tests} ./rt.summary.sh
else

  sbatch --nodes=1 -t 0:01:00 -A $ACCOUNT -J summary -o /dev/null -e /dev/null \
       --partition=$PARTITION --open-mode=append -q $QUEUE -d afterok${all_tests} ./rt.summary.sh

fi

exit
