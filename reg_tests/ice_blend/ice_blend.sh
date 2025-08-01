#!/bin/bash

#-----------------------------------------------------------------------------
# Run ice_blend consistency test.  This script called from it machine-specific
# driver script.
#-----------------------------------------------------------------------------

set -x

rm -fr $DATA
mkdir -p $DATA
cd $DATA

IMS_FILE=${HOMEreg}/input_data/imssnow96.grib2.gdas.2018120618
cp ${IMS_FILE} ./ims.grib2

# Follow OPS procedure by converting the IMS data to the 5 minute
# global grid.
$WGRIB2 ims.grib2 -match "ICEC" -grib ims.icec.grib2
grid173="0 0 0 0 0 0 0 0 4320 2160 0 0 89958000 42000 48 -89958000 359958000 83000 83000 0"
$COPYGB2 -x -i3 -g "$grid173" ims.icec.grib2 ims.icec.5min.grib2

FIVE_MIN_ICE_FILE=${HOMEreg}/input_data/seaice.5min.grib2.gdas.2018120618
FIVE_MIN_ICE_MASK_FILE=${HOMEgfs}/fix/am/emcsfc_gland5min.grib2
BLENDED_ICE_FILE="seaice.5min.blend"

# These are input files.
export FORT17="$FIVE_MIN_ICE_MASK_FILE"
export FORT11="ims.icec.5min.grib2"
export FORT15="$FIVE_MIN_ICE_FILE"

# This is the output blended file
export FORT51=${BLENDED_ICE_FILE}

# Run the emcsfc_ice_blend executable.
${HOMEgfs}/exec/emcsfc_ice_blend >> OUTPUT 2> errfile
iret=$?

if [ $iret -ne 0 ]; then
  set +x
  echo "<<< ICE_BLEND TEST FAILED. <<<"
  echo "<<< ICE_BLEND TEST FAILED. <<<"  > ${HOMEgfs}/reg_tests/ice_blend/summary.log
  exit $iret
else
# Follow the OPS procedure by converting the output to grib1 and replace the
# bitmap with a flag value of 1.57. That flag is expected by the global_cycle
# program.
  $WGRIB2 -set_int 3 51 42000 ${BLENDED_ICE_FILE} -grib ${BLENDED_ICE_FILE}.corner
  $CNVGRIB -g21 ${BLENDED_ICE_FILE}.corner ${BLENDED_ICE_FILE}.bitmap
  $COPYGB -M "#1.57" -x ${BLENDED_ICE_FILE}.bitmap $BLENDED_ICE_FILE
fi

cmp ${BLENDED_ICE_FILE} ${HOMEreg}/baseline_data/t1534/seaice.5min.blend
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
  echo "<<< ICE BLEND TEST FAILED. >>>" > ${HOMEgfs}/reg_tests/ice_blend/summary.log
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    ${HOMEgfs}/reg_tests/update_baseline.sh $HOMEreg "t1534" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< ICE BLEND TEST PASSED. >>>"
  echo "*********************************"
  echo "<<< ICE BLEND TEST PASSED. >>>" > ${HOMEgfs}/reg_tests/ice_blend/summary.log
fi

exit 0
