#!/bin/bash

#-----------------------------------------------------------------------------
# Run ice_blend consistency test.  This script called from it machine-specific
# driver script.
#-----------------------------------------------------------------------------

export IMS_FILE=$HOMEreg/input_data/imssnow96.grib2.gdas.2018120618
export FIVE_MIN_ICE_FILE=$HOMEreg/input_data/seaice.5min.grib2.gdas.2018120618

${HOMEgfs}/ush/emcsfc_ice_blend.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< ICE_BLEND TEST FAILED. <<<"
  echo "<<< ICE_BLEND TEST FAILED. <<<"  > ./summary.log
  exit $iret
fi

cmp ${DATA}/seaice.5min.blend $HOMEreg/baseline_data/seaice.5min.blend
iret=$?
test_failed=0
if [ $iret -ne 0 ]; then
  test_failed=1
fi

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "*********************************"
  echo "<<< ICE BLEND TEST FAILED. >>>"
  echo "*********************************"
  echo "<<< ICE BLEND TEST FAILED. >>>" > ./summary.log
else
  echo
  echo "*********************************"
  echo "<<< ICE BLEND TEST PASSED. >>>"
  echo "*********************************"
  echo "<<< ICE BLEND TEST PASSED. >>>" > ./summary.log
fi

exit 0
