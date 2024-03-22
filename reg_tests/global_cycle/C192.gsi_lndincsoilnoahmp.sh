#!/bin/bash

#------------------------------------------------------------------
# Run global_cycle for a C192 case to test the ingest and
# application of soil moisture and temperature increments
# from the GSI, into Noah-MP restarts.
# Compare output to a baseline set of files using the 'nccmp' 
# utility.
#------------------------------------------------------------------

set -x

NCCMP=${NCCMP:-$(which nccmp)}

export MAX_TASKS_CY=6

export HOMEgfs=$NWPROD

export FIXgfs=$HOMEreg/fix

export CYCLEXEC=$HOMEgfs/exec/global_cycle

export CDATE=2019073000
export FHOUR=00
export DELTSFC=6

export CASE=C192
export OCNRES=99999

export COMIN=$HOMEreg/input_data_noahmp

export JCAP=1534
export LONB=3072
export LATB=1536

export DONST="NO"
export use_ufo=.true.

export DO_SFCCYCLE=".FALSE." 
export DO_LNDINC=".TRUE." 
export DO_SOI_INC_GSI=".true."

export VERBOSE=YES
export CYCLVARS=FSNOL=-2.,FSNOS=99999.,

$HOMEgfs/ush/global_cycle_driver.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C192 GSI based LANDINC SOIL NOAHMP CYCLE TEST FAILED. >>>"
  exit $iret
fi

test_failed=0

cd $DATA
for files in *tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c192.gsi_lndincsoilnoahmp/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "**********************************************"
  echo "<<< C192 GSI based LANDINC SOIL-NOAHMP CYCLE TEST FAILED. >>>"
  echo "**********************************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEgfs/reg_tests/update_baseline.sh $HOMEreg "c192.lndincsoilnoahmp" $commit_num
  fi
else
  echo
  echo "*****************************************"
  echo "<<< C192 GSI based LANDINC SOIL-NOAHMP CYCLE TEST PASSED. >>>"
  echo "*****************************************"
fi

exit
