#!/bin/bash

#--------------------------------------------------------------------------
# Mimic GFS OPS.  This script is run from its machine-specific
# driver.
#--------------------------------------------------------------------------

set -x

export IMS_FILE=$HOMEreg/input_data/imssnow96.grb
export AFWA_NH_FILE=$HOMEreg/input_data/NPR.SNWN.SP.S1200.MESH16
export AFWA_SH_FILE=$HOMEreg/input_data/NPR.SNWS.SP.S1200.MESH16

export MODEL_LATITUDE_FILE=$HOMEgfs/fix/fix_am/global_latitudes.t1534.3072.1536.grb
export MODEL_LONGITUDE_FILE=$HOMEgfs/fix/fix_am/global_longitudes.t1534.3072.1536.grb
export MODEL_SLMASK_FILE=$HOMEgfs/fix/fix_am/global_slmask.t1534.3072.1536.grb
export GFS_LONSPERLAT_FILE=$HOMEgfs/fix/fix_am/global_lonsperlat.t1534.3072.1536.txt

export OMP_NUM_THREADS=1
export OUTPUT_GRIB2=.false.

${HOMEgfs}/ush/emcsfc_snow.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< SNOW2MDL OPS TEST FAILED. <<<"
  exit $iret
fi

test_failed=0

cmp ${DATA}/snogrb_model $HOMEreg/baseline_data/t1534.ops/snogrb_model
iret=$?
if [ $iret -ne 0 ]; then
  test_failed=1
fi

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL OPS TEST FAILED. >>>"
  echo "*********************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    cd $DATA
    $HOMEgfs/reg_tests/update_baseline.sh $HOMEreg "t1534.ops" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL OPS TEST PASSED. >>>"
  echo "*********************************"
fi

exit
