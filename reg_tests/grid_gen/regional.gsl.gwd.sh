#!/bin/bash

#-----------------------------------------------------------------------
# Create a regional esg grid with GSL gravity wave drag files.  
# Compare output to a set of baseline files using the 'nccmp' utility.
# This script is run by the machine specific driver script.
#-----------------------------------------------------------------------

set -x

export TEMP_DIR=${WORK_DIR}/regional.gsl.gwd.work
export out_dir=${WORK_DIR}/regional.gsl.gwd

export gtype=regional_esg
export make_gsl_orog=true    # Create GSL gravity wave drag fields
export target_lon=-97.5      # Center longitude of the highest resolution tile
export target_lat=35.5       # Center latitude of the highest resolution tile
export idim=301              # Dimension of grid in 'i' direction
export jdim=200              # Dimension of grid in 'j' direction
export delx=0.0585           # Grid spacing in degrees in 'i' direction
export dely=0.0585           # Grid spacing in degrees in 'j' direction
export halo=4

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< REGIONAL GSL GWD TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir/C772

test_failed=0
for files in *tile*.nc ./sfc/*tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/regional.gsl.gwd/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< REGIONAL GSL GWD TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $home_dir/reg_tests/update_baseline.sh "${HOMEreg}/.." "regional.gsl.gwd" $commit_num
  fi
else
  echo "<<< REGIONAL GSL GWD TEST PASSED. >>>"
fi

exit 0
