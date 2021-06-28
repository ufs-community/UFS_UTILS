#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C96 coldstart files using FV3 tiled restart files
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c96_fv3_restart
rm -fr $DATA

export FIXfv3=${HOMEreg}/fix/C96
export COMIN=${HOMEreg}/input_data/fv3.restart
export VCOORD_FILE=${HOMEufs}/fix/fix_am/global_hyblev.l64.txt
export INPUT_TYPE='restart'
export MOSAIC_FILE_INPUT_GRID="${HOMEreg}/fix/C384/C384_mosaic.nc"
export OROG_DIR_INPUT_GRID=${HOMEreg}/fix/C384

# Note: no double quotes at the beginning or end.
export OROG_FILES_INPUT_GRID='C384_oro_data.tile1.nc","C384_oro_data.tile2.nc","C384_oro_data.tile3.nc","C384_oro_data.tile4.nc","C384_oro_data.tile5.nc","C384_oro_data.tile6.nc'
export ATM_CORE_FILES_INPUT='20190706.120000.fv_core.res.tile1.nc","20190706.120000.fv_core.res.tile2.nc","20190706.120000.fv_core.res.tile3.nc","20190706.120000.fv_core.res.tile4.nc","20190706.120000.fv_core.res.tile5.nc","20190706.120000.fv_core.res.tile6.nc","20190706.120000.fv_core.res.nc'
export ATM_TRACER_FILES_INPUT='20190706.120000.fv_tracer.res.tile1.nc","20190706.120000.fv_tracer.res.tile2.nc","20190706.120000.fv_tracer.res.tile3.nc","20190706.120000.fv_tracer.res.tile4.nc","20190706.120000.fv_tracer.res.tile5.nc","20190706.120000.fv_tracer.res.tile6.nc'
export SFC_FILES_INPUT='20190706.120000.sfc_data.tile1.nc","20190706.120000.sfc_data.tile2.nc","20190706.120000.sfc_data.tile3.nc","20190706.120000.sfc_data.tile4.nc","20190706.120000.sfc_data.tile5.nc","20190706.120000.sfc_data.tile6.nc'

export TRACERS_TARGET='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'
export TRACERS_INPUT='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'

export CDATE=2019070612

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
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_fv3_restart/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 FV3 RESTART TEST FAILED. >>>"
else
  echo "<<< C96 FV3 RESTART TEST PASSED. >>>"
fi

exit 0
