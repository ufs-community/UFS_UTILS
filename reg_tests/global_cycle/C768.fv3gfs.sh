#!/bin/bash

#------------------------------------------------------------------
# Run global_cycle for a C768 test case.  Compare output
# to a baseline set of files using the 'nccmp' utility.
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

export COMIN=$HOMEreg/input_data
export FNTSFA=$COMIN/gdas.t00z.rtgssthr.grb
export FNSNOA=$COMIN/gdas.t00z.snogrb_t1534.3072.1536
export FNACNA=$COMIN/gdas.t00z.seaice.5min.blend.grb
export NST_FILE=$COMIN/gdas.t00z.dtfanl.nc

export JCAP=1534
export LONB=3072
export LATB=1536

export OROFIX=$BASE_GSM/fix/orog/$CASE
export OROFILE=${CASE}_oro_data

export FIXgsm=$BASE_GSM/fix/am

export FNALBC2=$FIXgsm/global_albedo4.1x1.grb
export FNALBC=$FIXgsm/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb
export FNALBC=$FIXgsm/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb
export FNTG3C=$FIXgsm/global_tg3clim.2.6x1.5.grb
export FNVETC=$FIXgsm/global_vegtype.igbp.t1534.3072.1536.rg.grb
export FNSOTC=$FIXgsm/global_soiltype.statsgo.t1534.3072.1536.rg.grb
export FNVEGC=$FIXgsm/global_vegfrac.0.144.decpercent.grb
export FNVMNC=$FIXgsm/global_shdmin.0.144x0.144.grb
export FNVMXC=$FIXgsm/global_shdmax.0.144x0.144.grb
export FNSLPC=$FIXgsm/global_slope.1x1.grb
export FNABSC=$FIXgsm/global_mxsnoalb.uariz.t1534.3072.1536.rg.grb
export FNAISC=$FIXgsm/CFSR.SEAICE.1982.2012.monthly.clim.grb

export FIXfv3=$HOMEreg/fix

export DONST="YES"
export use_ufo=.true.

export VERBOSE=YES
export CYCLVARS=FSNOL=-2.,FSNOS=99999.,

$BASE_GSM/ush/global_cycle_driver.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C768 GLOBAL CYCLE TEST FAILED. >>>"
  exit $iret
fi

test_failed=0

cd $DATA
for files in *tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c768.fv3gfs/$files
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
  echo "<<< C768 GLOBAL CYCLE TEST FAILED. >>>"
  echo "*********************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $BASE_GSM/reg_tests/update_baseline.sh $HOMEreg "c768.fv3gfs" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< C768 GLOBAL CYCLE TEST PASSED. >>>"
  echo "*********************************"
fi

exit
