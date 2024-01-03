#!/bin/bash

#-----------------------------------------------------------------------
# Create a regional esg grid. Output dominant soil and vegetation
# categories and well as the percentage of each category. 
# Compare output to a set of baseline files using the 'nccmp'
# utility.  This script is run by the machine specific driver script.
#-----------------------------------------------------------------------

set -x

TEST_NAME="esg.regional.pct.cat"
export TEMP_DIR=${WORK_DIR}/${TEST_NAME}.work
export out_dir=${WORK_DIR}/${TEST_NAME}

export gtype=regional_esg
export target_lon=-97.5      # Center longitude of the highest resolution tile
export target_lat=35.5       # Center latitude of the highest resolution tile
export idim=1301             # Dimension of grid in 'i' direction
export jdim=600              # Dimension of grid in 'j' direction
export delx=0.0145           # Grid spacing in degrees in 'i' direction
export dely=0.0145           # Grid spacing in degrees in 'j' direction
export halo=3
export vegsoilt_frac=.true.  # Output dominant soil/veg categories as well
                             # as the percentage of each category.

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< ESG REGIONAL PERCENT CATEGORY TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir/C3113

test_failed=0
for files in *tile*.nc ./sfc/*tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/${TEST_NAME}/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< ESG REGIONAL PERCENT CATEGORY TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $home_dir/reg_tests/update_baseline.sh "${HOMEreg}/.." "${TEST_NAME}" $commit_num
  fi
else
  echo "<<< ESG REGIONAL PERCENT CATEGORY TEST PASSED. >>>"
fi

exit 0
