#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C96 coldstart files using GFS sigio/sfcio files
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c96_gfs_sigio
rm -fr $DATA

export FIXfv3=${HOMEreg}/fix/C96
export COMIN=${HOMEreg}/input_data/gfs.sigio
export ATM_FILES_INPUT=gdas.t00z.sanl
export SFC_FILES_INPUT=gdas.t00z.sfcanl
export CONVERT_NST='.false.'
export VCOORD_FILE=${HOMEufs}/fix/fix_am/global_hyblev.l64.txt
export INPUT_TYPE="gfs_spectral"

# dont start/end with double quotes
export TRACERS_TARGET='"sphum","o3mr","liq_wat"'
export TRACERS_INPUT='"spfh","o3mr","clwmr"'
export CDATE=2017071700
export OMP_NUM_THREADS_CY=6

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  echo "<<< C96 GFS SIGIO TEST FAILED. <<<"
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
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_gfs_sigio/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 GFS SIGIO TEST FAILED. >>>"
else
  echo "<<< C96 GFS SIGIO TEST PASSED. >>>"
fi

exit 0
