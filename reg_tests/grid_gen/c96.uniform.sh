#!/bin/bash

set -x

export TMPDIR=${WORK_DIR}/c96.uniform.work
export out_dir=${WORK_DIR}/c96.uniform

export res=96
export gtype=uniform

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

echo "Starting at: " `date`

$home_dir/ush/fv3gfs_driver_grid.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C96 UNIFORM TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output to baseline set of data.
#-----------------------------------------------------------------------------

cd $out_dir

test_failed=0
for files in *.nc ./fix_sfc/*.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96.uniform/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 UNIFORM TEST PASSED. >>>"
else
  echo "<<< C96 UNIFORM TEST PASSED. >>>"
fi

exit 0
