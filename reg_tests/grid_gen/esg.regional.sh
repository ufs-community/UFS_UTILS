#!/bin/bash

#-----------------------------------------------------------------------
# Create a regional esg grid.  Compare output to a set
# of baseline files using the 'nccmp' utility.  This script is
# run by the machine specific driver script.
#-----------------------------------------------------------------------

set -x

export TEMP_DIR=${WORK_DIR}/esg.regional.work
export out_dir=${WORK_DIR}/esg.regional

export gtype=regional_esg
export target_lon=-97.5      # Center longitude of the highest resolution tile
export target_lat=35.5       # Center latitude of the highest resolution tile
export idim=1301             # Dimension of grid in 'i' direction
export jdim=600              # Dimension of grid in 'j' direction
export delx=0.0145           # Grid spacing in degrees in 'i' direction
export dely=0.0145           # Grid spacing in degrees in 'j' direction
export halo=3

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< ESG REGIONAL TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir/C3113

test_failed=0
for files in *tile*.nc ./fix_sfc/*tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/esg.regional/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< ESG REGIONAL TEST FAILED. >>>"
else
  echo "<<< ESG REGIONAL TEST PASSED. >>>"
fi

exit 0
