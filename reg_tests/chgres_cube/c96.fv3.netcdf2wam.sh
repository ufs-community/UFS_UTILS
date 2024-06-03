#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C96 WAM coldstart files using FV3 gaussian netcdf 
# (gfs v16) as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c96_fv3_netcdf2wam
rm -fr $DATA

export CRES=96
export ocn=100
export FIXfv3=${HOMEreg}/fix/C${CRES}

export COMIN=${HOMEreg}/input_data/fv3.netcdf
export ATM_FILES_INPUT=gfs.t00z.atmf000.nc
export VCOORD_FILE=${HOMEufs}/fix/am/global_hyblev.l64.txt
export INPUT_TYPE="gaussian_netcdf"
export CONVERT_SFC=".false."
export CONVERT_NST=".false."
export WAM_PARM_FILE=${HOMEufs}/parm/msis_lib/msis21.parm

export CDATE=2020020200

# export TRACERS_INPUT='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'
export TRACERS_TARGET='"sphum","liq_wat","spo3","ice_wat","rainwat","snowwat","graupel","spo","spo2"'

export VCOORD_FILE=${HOMEufs}/fix/am/global_hyblev.l150.txt
export WAM_COLD_START=.true.

export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS:-1}

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------------
# Invoke chgres program.
#-----------------------------------------------------------------------------

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C96 FV3 GAUSSIAN NETCDF2WAM TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output from chgres to baseline set of data.
#-----------------------------------------------------------------------------

cd $DATA

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_fv3_netcdf2wam/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 FV3 GAUSSIAN NETCDF2WAM TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEufs/reg_tests/update_baseline.sh $HOMEreg "c96_fv3_netcdf2wam" $commit_num
  fi
else
  echo "<<< C96 FV3 GAUSSIAN NETCDF2WAM TEST PASSED. >>>"
fi

exit 0
