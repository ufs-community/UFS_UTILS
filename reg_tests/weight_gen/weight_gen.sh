#!/bin/bash

set -x

#-----------------------------------------------------------------------------
# Run weight consistency test.  This script called from its machine-specific
# driver script.
#-----------------------------------------------------------------------------

rm -fr $DATA
mkdir -p $DATA
cd $DATA

HOMEexec=$HOMEufs/exec
HOMEtest=$HOMEufs/reg_tests
HOMEtest_wg=$HOMEtest/weight_gen

# Run executable
${HOMEexec}/weight_gen "C48"

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< WEIGHT_GEN TEST FAILED. >>>"
  echo "<<< WEIGHT_GEN TEST FAILED. >>>"  > $HOMEtest_wg/summary.log
  exit $iret
fi

NCCMP=${NCCMP:-$(which nccmp)}

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/C48/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "*********************************"
  echo "<<< WEIGHT_GEN TEST FAILED. >>>"
  echo "*********************************"
  echo "<<< WEIGHT_GEN TEST FAILED. >>>" > $HOMEtest_wg/summary.log
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    cd $DATA
    $HOMEtest/update_baseline.sh $HOMEreg "C48" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< WEIGHT_GEN TEST PASSED. >>>"
  echo "*********************************"
  echo "<<< WEIGHT_GEN TEST PASSED. >>>" > $HOMEtest_wg/summary.log
fi

exit 0
