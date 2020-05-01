#!/bin/bash

#-----------------------------------------------------------------------
# Create a C96 regional grid.  Compare output to a set
# of baseline files using the 'nccmp' utility.  This script is
# run by the machine specific driver script.
#-----------------------------------------------------------------------

set -x

export TMPDIR=${WORK_DIR}/c96.regional.work
export out_dir=${WORK_DIR}/c96.regional

export res=96
export gtype=regional
export stretch_fac=1.5       # Stretching factor for the grid
export target_lon=-97.5      # Center longitude of the highest resolution tile
export target_lat=35.5       # Center latitude of the highest resolution tile
export refine_ratio=3        # The refinement ratio
export istart_nest=27        # Starting i-direction index of nest grid in parent tile supergrid
export jstart_nest=37        # Starting j-direction index of nest grid in parent tile supergrid
export iend_nest=166         # Ending i-direction index of nest grid in parent tile supergrid
export jend_nest=164         # Ending j-direction index of nest grid in parent tile supergrid
export halo=3

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C96 REGIONAL TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir

test_failed=0
for files in *tile*.nc ./fix_sfc/*tile*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/c96.regional/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 REGIONAL TEST FAILED. >>>"
else
  echo "<<< C96 REGIONAL TEST PASSED. >>>"
fi

exit 0
