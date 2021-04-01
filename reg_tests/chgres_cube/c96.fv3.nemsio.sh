#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C96 coldstart files using FV3 gaussian nemsio files
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c96_fv3_nemsio
rm -fr $DATA

export FIXfv3=${HOMEreg}/fix/C96
export COMIN=${HOMEreg}/input_data/fv3.nemsio
export ATM_FILES_INPUT=gfs.t12z.atmf000.nemsio
export SFC_FILES_INPUT=gfs.t12z.sfcf000.nemsio
export VCOORD_FILE=${HOMEufs}/fix/fix_am/global_hyblev.l64.txt

export CDATE=2019070412

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
  echo "<<< C96 FV3 GAUSSIAN NEMSIO TEST FAILED. <<<"
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
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_fv3_nemsio/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 FV3 GAUSSIAN NEMSIO TEST FAILED. >>>"
else
  echo "<<< C96 FV3 GAUSSIAN NEMSIO TEST PASSED. >>>"
fi

exit 0
