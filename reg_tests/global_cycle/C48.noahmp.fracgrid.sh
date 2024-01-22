#!/bin/bash

#------------------------------------------------------------------
# Run global_cycle for a C48 case that tests the NOAHMP and
# fractional grid options.  
#
# Compare output to a baseline set of files using the 'nccmp'
# utility.
#------------------------------------------------------------------

set -x

NCCMP=${NCCMP:-$(which nccmp)}

export MAX_TASKS_CY=6

export HOMEgfs=$NWPROD
export BASE_GSM=$NWPROD

export CYCLEXEC=$BASE_GSM/exec/global_cycle

export CDATE=2021032406
export FHOUR=00
export DELTSFC=6

export CASE=C48
export OCNRES=500

export COMIN=$HOMEreg/input_data_c48.noahmp.frac.grid
export FNACNA=$COMIN/gdas.t06z.seaice.5min.blend.grb
export FNTSFA=" "
export FNSNOA=" "
export NST_FILE=$COMIN/gdas.t06z.dtfanl.nc

export JCAP=1534
export LONB=3072
export LATB=1536

export OROFIX=$HOMEreg/fix/$CASE

export FIXgsm=$BASE_GSM/fix/am

export FNAISC=$FIXgsm/IMS-NIC.blended.ice.monthly.clim.grb

export DONST="YES"
export use_ufo=.true.
export FRAC_GRID=.true.

export VERBOSE=YES
export CYCLVARS=FSNOL=99999.,FSNOS=99999.,

$BASE_GSM/ush/global_cycle_driver.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C48 NOAHMP FRAC GRID TEST FAILED. >>>"
  exit $iret
fi

test_failed=0

cd $DATA
for files in *tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c48.noahmp.fracgrid/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "******************************************"
  echo "<<< C48 NOAHMP FRAC GRID TEST FAILED. >>>"
  echo "******************************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $BASE_GSM/reg_tests/update_baseline.sh $HOMEreg "c48.noahmp.fracgrid" $commit_num
  fi
else
  echo
  echo "*****************************************"
  echo "<<< C48 NOAHMP FRAC GRID TEST PASSED. >>>"
  echo "*****************************************"
fi

exit
