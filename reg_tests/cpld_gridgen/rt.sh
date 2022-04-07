#!/bin/bash
set -eu

error() {
  echo
  echo "$@" 1>&2
  exit 1
}

usage() {
  echo
  echo "Usage: $program [-c] [-m] [-h]"
  echo
  echo "  -c create a new baseline"
  echo
  echo "  -m compare against the new baseline"
  echo
  echo "  -h display this help and exit"
  echo
  echo "  Examples"
  echo
  echo "    './rt.sh'    compare against the existing baseline"
  echo "    './rt.sh -c' create a new baseline"
  echo "    './rt.sh -m' compare against the new baseline"
  echo
}

usage_and_exit() {
  usage
  exit $1
}

check_results() {

  [ -o xtrace ] && set_x='set -x' || set_x='set +x'
  set +x

  local test_status=PASS
  # verification run
  if [[ $CREATE_BASELINE = false ]]; then
  
    echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo "Baseline dir = $BASELINE" | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo "Checking test $TEST_NAME results ...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

    for file in $BASELINE/*.nc; do
      printf %s "Comparing " $(basename ${file}) "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

      if [[ ! -f $RUNDIR/$(basename ${file}) ]]; then
        echo "....MISSING file" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        test_status=FAIL
      else
        nccmp -dmfqS $(basename ${file}) $file && d=$? || d=$?
        if [[ $d -ne 0 ]]; then
          echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
          test_status=FAIL
        else
          echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        fi
      fi
    done

  # baseline creation run
  else

    echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo "Moving baseline files...." | tee -a $PATHRT/$REGRESSIONTEST_LOG
    echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

    mkdir -p $NEW_BASELINE

    for file in *.nc; do
      printf %s "Moving " $file "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

      cp $file $NEW_BASELINE/$file && d=$? || d=$?
      if [[ $d -ne 0 ]]; then
        echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        test_status=FAIL
      else
        echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
      fi
    done

  fi

  if [[ $test_status == FAIL ]]; then
    echo "$TEST_NAME failed" >> $PATHRT/fail_test_$TEST_NAME
  fi
}

readonly program=$(basename $0)
# PATHRT - Path to regression tests directory
readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
# PATHTR - Path to the UFS UTILS directory
readonly PATHTR="$(cd $PATHRT/../.. && pwd)"
TESTS_FILE="$PATHRT/rt.conf"
export TEST_NAME=

cd $PATHRT
export compiler=${compiler:-intel}
source $PATHTR/sorc/machine-setup.sh >/dev/null 2>&1
echo "Machine: $target"
echo "Compiler: $compiler"

COMPILE_LOG=compile.log
RUNTEST_LOG=run.log
REGRESSIONTEST_LOG=RegressionTests_$target.$compiler.log
rm -f fail_test* $COMPILE_LOG

if [[ $target = hera ]]; then
  STMP=/scratch1/NCEPDEV/stmp4
  BASELINE_ROOT=$STMP/$USER/UFS_UTILS_BASELINE
  NEW_BASELINE_ROOT=$STMP/$USER/CPLD_GRIDGEN/BASELINE
  RUNDIR_ROOT=$STMP/$USER/CPLD_GRIDGEN/rt_$$
elif [[ $target = orion ]]; then
  STMP=
elif [[ $target = jet ]]; then
  STMP=
fi

CREATE_BASELINE=false
while getopts :cmh opt; do
  case $opt in
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

# Build the executable file
cd $PATHTR
rm -rf $PATHTR/build $PATHTR/exec $PATHTR/lib
./build_all.sh >$PATHRT/$COMPILE_LOG 2>&1
if [[ ! -f $PATHTR/exec/cpld_gridgen ]]; then
  error "Build did not generate the cpld_gridgen exe file. Check $PATHTR/build/"
else
  echo "Build was successful. cpld_gridgen executable file is in $PATHTR/exec/"
fi

module use $PATHTR/modulefiles
module load build.$target.$compiler
module list

if [[ $CREATE_BASELINE = true ]]; then
  rm -rf $NEW_BASELINE_ROOT
  mkdir -p $NEW_BASELINE_ROOT
fi

date > $PATHRT/$REGRESSIONTEST_LOG
echo "Start Regression test" | tee -a  $PATHRT/$REGRESSIONTEST_LOG

# Run tests specified in $TESTS_FILE
while read -r line || [ "$line" ]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line =~ \# ]] && continue

  TEST_NAME=$(echo $line | cut -d'|' -f1 | sed -e 's/^ *//' -e 's/ *$//')
  DEP_NAME=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
  TEST_NAME=${TEST_NAME##*_}
  DEP_NAME=${DEP_NAME##*_}

  cd $PATHRT
  RUNDIR=$RUNDIR_ROOT/$TEST_NAME
  BASELINE=$BASELINE_ROOT/$TEST_NAME
  NEW_BASELINE=$NEW_BASELINE_ROOT/$TEST_NAME
  DEPDIR=$RUNDIR_ROOT/$DEP_NAME
  mkdir -p $RUNDIR
  export OUTDIR_PATH=$RUNDIR

  if [[ -n $DEP_NAME ]]; then
    cp $DEPDIR/Ct.mx025_SCRIP.nc $RUNDIR >/dev/null 2>&1 && d=$? || d=$?
    if [[ $d -eq 1 ]]; then    
      error "DEPDIR does not exist. Dependency not met"
    fi
  fi

  cp $PATHTR/exec/cpld_gridgen $RUNDIR
  cp $PATHTR/ush/cpld_gridgen.sh $RUNDIR
  cp $PATHRT/parm/grid.nml.IN $RUNDIR
  cd $RUNDIR

  ./cpld_gridgen.sh $TEST_NAME >$PATHRT/$RUNTEST_LOG 2>&1

  check_results

done <$TESTS_FILE
if [[ $? -ne 0 ]]; then
  error "Run test while loop did not finish properly"
fi

cd $PATHRT
FAIL_FILES="fail_test_*"
for file in $FAIL_FILES; do
  if [[ -f "$file" ]]; then
    cat "$file" >> fail_test
  fi
done

if [[ -e fail_test ]]; then
  echo | tee -a $REGRESSIONTEST_LOG
  for file in fail_test_*; do
    cat $file >>$REGRESSIONTEST_LOG 
  done

  echo | tee -a $REGRESSIONTEST_LOG
  echo "REGRESSION TEST FAILED" | tee -a $REGRESSIONTEST_LOG
else
  echo | tee -a $REGRESSIONTEST_LOG
  echo "REGRESSION TEST WAS SUCCESSFUL" | tee -a $REGRESSIONTEST_LOG
fi
date >> $REGRESSIONTEST_LOG
