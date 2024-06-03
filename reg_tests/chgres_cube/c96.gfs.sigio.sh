#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C96 coldstart files using GFS sigio/sfcio files
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------------

set -x

# Orion won't let me set the ulimit in the driver script.  Set it here.
machine=${machine:-NULL}
if [ $machine == 'orion' ]; then
  ulimit -s 199000000
fi

export DATA=$OUTDIR/c96_gfs_sigio
rm -fr $DATA

export CRES=96
export ocn=100
export FIXfv3=${HOMEreg}/fix/C${CRES}

export COMIN=${HOMEreg}/input_data/gfs.sigio
export ATM_FILES_INPUT=gdas.t00z.sanl
export SFC_FILES_INPUT=gdas.t00z.sfcanl
export CONVERT_NST='.false.'
export VCOORD_FILE=${HOMEufs}/fix/am/global_hyblev.l64.txt
export INPUT_TYPE="gfs_sigio"

# dont start/end with double quotes
export TRACERS_TARGET='"sphum","o3mr","liq_wat"'
export TRACERS_INPUT='"spfh","o3mr","clwmr"'
export CDATE=2017071700
export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS:-1}

NCCMP=${NCCMP:-$(which nccmp)}

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
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
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEufs/reg_tests/update_baseline.sh $HOMEreg "c96_gfs_sigio" $commit_num
  fi
else
  echo "<<< C96 GFS SIGIO TEST PASSED. >>>"
fi

exit 0
