#!/bin/bash

#--------------------------------------------------------------------------
# Create a snow file from afwa global data and ims data.  This script
# is run from its machine-specific driver.
#
# This test mimics current GFS OPS which uses the global afwa data.
#
# Note, this test uses the "snow2mdl.nml.tmpl" template to create
# the fort.41 namelist as is done by the global workflow.
#--------------------------------------------------------------------------

echo "BEGIN SNOW2MDL GLOBAL TEST."

set -x

HOMEush="${HOMEglobal}/ush"
HOMEparm="${HOMEglobal}/parm"
HOMEexec="${HOMEglobal}/exec"
HOMEfix="${HOMEglobal}/fix/am"

source "${HOMEush}/atparse.bash"  # include function atparse for parsing @[XYZ] templated files

SNOW2MDLNMLTMPL="${HOMEparm}/prep_sfc/snow2mdl.nml.tmpl"

CLIMO_QC="${HOMEfix}/emcsfc_snow_cover_climo.grib2"

MODEL_LATITUDE_FILE="${HOMEfix}/global_latitudes.t1534.3072.1536.grb"
MODEL_LONGITUDE_FILE="${HOMEfix}/global_longitudes.t1534.3072.1536.grb"
MODEL_SLMASK_FILE="${HOMEfix}/global_slmask.t1534.3072.1536.grb"
GFS_LONSPERLAT_FILE="${HOMEfix}/global_lonsperlat.t1534.3072.1536.txt"

MODEL_SNOW_FILE="snogrb_model"
OUTPUT_GRIB2=".false."

IMSYEAR=2025
IMSMONTH=03
IMSDAY=25
IMSHOUR=0

rm -fr $DATA
mkdir -p $DATA
cd $DATA

cp "$HOMEreg/input_data/global/gfs.t00z.imssnow96.grib2" "./imssnow96.grib2"
cp "$HOMEreg/input_data/global/gfs.t00z.snow.usaf.grib2" "./snow.usaf.grib2"

atparse < "${SNOW2MDLNMLTMPL}" >> "./fort.41"
echo "Rendered fort.41"
cat "./fort.41"

eval ${HOMEexec}/emcsfc_snow2mdl >> OUTPUT 2> errfile
iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< SNOW2MDL GLOBAL TEST FAILED. <<<"
  exit $iret
fi

test_failed=0

cmp ${DATA}/snogrb_model $HOMEreg/baseline_data/t1534.global/snogrb_model
iret=$?
if [ $iret -ne 0 ]; then
  test_failed=1
fi

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL GLOBAL TEST FAILED. >>>"
  echo "*********************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    cd $DATA
    $HOMEglobal/reg_tests/update_baseline.sh $HOMEreg "t1534.global" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL GLOBAL TEST PASSED. >>>"
  echo "*********************************"
fi

exit
