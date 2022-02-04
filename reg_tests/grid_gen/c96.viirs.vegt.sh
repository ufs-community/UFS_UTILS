#!/bin/bash

#-----------------------------------------------------------------------
# Create a C96 global uniform grid using VIIRS-based vegetation type
# data.  Compare output to a set of baseline files using the
# 'nccmp' utility.  This script is run by the machine specific
# driver script.
#-----------------------------------------------------------------------

set -x

export TEMP_DIR=${WORK_DIR}/c96.viirs.vegt.work
export out_dir=${WORK_DIR}/c96.viirs.vegt

export res=96
export gtype=uniform
export veg_type_src="viirs.igbp.0.05"

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C96 VIIRS VEGT FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir/C96

test_failed=0
for files in *tile*.nc ./fix_sfc/*tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/c96.viirs.vegt/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 VIIRS VEGT TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $home_dir/reg_tests/update_baseline.sh "${HOMEreg}/.." "c96.viirs.vegt" $commit_num
  fi
else
  echo "<<< C96 VIIRS VEGT TEST PASSED. >>>"
fi

exit 0
