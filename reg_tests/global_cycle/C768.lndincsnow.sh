#!/bin/bash

#------------------------------------------------------------------
# Run global_cycle for a C768 case to test the ingest of snow
# increments from JEDI.  Compare output to a baseline set of
# files using the 'nccmp' utility.
#------------------------------------------------------------------

set -x

NCCMP=${NCCMP:-$(which nccmp)}

export MAX_TASKS_CY=6

export HOMEgfs=$NWPROD
export BASE_GSM=$NWPROD

export CYCLEXEC=$BASE_GSM/exec/global_cycle

export CDATE=2019073000
export FHOUR=00
export DELTSFC=6

export CASE=C768
export OCNRES=99999

export COMIN=$HOMEreg/input_data
export FNTSFA=$COMIN/gdas.t00z.rtgssthr.grb
export FNSNOA=$COMIN/gdas.t00z.snogrb_t1534.3072.1536
export FNACNA=$COMIN/gdas.t00z.seaice.5min.blend.grb
export NST_FILE=$COMIN/gdas.t00z.dtfanl.nc

export DO_SNO_INC=.true. # must be lower-case.
export JCAP=1534
export LONB=3072
export LATB=1536

export OROFIX=$HOMEreg/fix/$CASE
export FIXgsm=$BASE_GSM/fix/am

export DONST="NO"
export use_ufo=.true.

export DO_SFCCYCLE=".FALSE." 
export DO_LNDINC=".TRUE." 

export VERBOSE=YES
export CYCLVARS=FSNOL=99999.,FSNOS=99999.,

$BASE_GSM/ush/global_cycle_driver.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C768 LANDINC SNOW CYCLE TEST FAILED. >>>"
  exit $iret
fi

test_failed=0

cd $DATA
for files in *tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c768.lndincsnow/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "****************************************"
  echo "<<< C768 LANDINC SNOW CYCLE TEST FAILED. >>>"
  echo "****************************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $BASE_GSM/reg_tests/update_baseline.sh $HOMEreg "c768.lndincsnow" $commit_num
  fi
else
  echo
  echo "***************************************"
  echo "<<< C768 LANDINC SNOW CYCLE TEST PASSED. >>>"
  echo "***************************************"
fi

exit
