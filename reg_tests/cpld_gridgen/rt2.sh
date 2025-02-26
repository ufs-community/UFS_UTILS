#!/bin/bash

set -ux

error() {
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
echo "Machine: $target"

MOM6_version=20250128

if [[ $target = hera ]]; then
  export MOM6_FIXDIR=/scratch1/NCEPDEV/global/glopara/fix/mom6/${MOM6_version}
  STMP=${STMP:-/scratch2/NCEPDEV/stmp1/$USER}
  export NCCMP=nccmp
  BASELINE_ROOT=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
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
     echo "IntelLLVM not available. Will use Intel Classic."
    compiler=intel
  fi
fi
export compiler
echo "Compiler: $compiler"

# Build the executable file
if [[ $BUILD_EXE = true ]]; then
    COMPILE_LOG=compile.log
    cd $PATHTR
    rm -rf $COMPILE_LOG $PATHTR/build $PATHTR/exec $PATHTR/lib
    ./build_all.sh >$PATHRT/$COMPILE_LOG 2>&1 && d=$? || d=$?
    if [[ d -ne 0 ]]; then
        error "Build did not finish successfully. Check $COMPILE_LOG"
    else
        echo "Build was successful"
        cd $PATHRT
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

rm -f fail_test* summary.log

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
  rm -f $REGRESSIONTEST_LOG
  rm -f nccmp_*log
  rm -f run_${TEST_NAME}.log

  cp $PATHRT/parm/grid.nml.IN $RUNDIR
  cp $PATHTR/exec/cpld_gridgen $RUNDIR
  
  tests[$i]=$(sbatch --parsable --ntasks-per-node=1 --nodes=1 -t 0:10:00 -A fv3-cpu -q batch -J $TEST_NAME \
            -o run_${TEST_NAME}.log -e run_${TEST_NAME}.log $PATHTR/ush/cpld_gridgen.sh "$TEST_NAME" "$ATMLIST")

  all_tests=${all_tests}":"${tests[$i]}

  ((i=i+1))
done < ./rt.conf

export target

sbatch --nodes=1 -t 0:01:00 -A fv3-cpu -J summary -o /dev/null -e /dev/null \
       --open-mode=append -q batch -d afterok${all_tests} ./rt.summary.sh

exit
