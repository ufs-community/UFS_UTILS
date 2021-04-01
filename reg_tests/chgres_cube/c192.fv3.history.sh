#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C192 coldstart files using FV3 tiled history files
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific 
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c192_fv3_history
rm -fr $DATA

export CRES=192
export FIXfv3=${HOMEreg}/fix/C192
export COMIN=${HOMEreg}/input_data/fv3.history

# Pay attention to the quotes.  Dont start/end with double quote.
export ATM_FILES_INPUT='dynf000.tile1.nc","dynf000.tile2.nc","dynf000.tile3.nc","dynf000.tile4.nc","dynf000.tile5.nc","dynf000.tile6.nc'
export SFC_FILES_INPUT='phyf000.tile1.nc","phyf000.tile2.nc","phyf000.tile3.nc","phyf000.tile4.nc","phyf000.tile5.nc","phyf000.tile6.nc'
export VCOORD_FILE=${HOMEufs}/fix/fix_am/global_hyblev.l64.txt
export INPUT_TYPE='history'
export MOSAIC_FILE_INPUT_GRID="${HOMEreg}/fix/C96/C96_mosaic.nc"
export OROG_DIR_INPUT_GRID=${HOMEreg}/fix/C96
export OROG_FILES_INPUT_GRID='C96_oro_data.tile1.nc","C96_oro_data.tile2.nc","C96_oro_data.tile3.nc","C96_oro_data.tile4.nc","C96_oro_data.tile5.nc","C96_oro_data.tile6.nc'
export TRACERS_TARGET='"sphum","liq_wat","o3mr"'
export TRACERS_INPUT='"spfh","clwmr","o3mr"'

export CDATE=2018100300

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
  echo "<<< C192 FV3 HISTORY TEST FAILED. <<<"
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
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c192_fv3_history/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C192 FV3 HISTORY TEST FAILED. >>>"
else
  echo "<<< C192 FV3 HISTORY TEST PASSED. >>>"
fi

exit 0
