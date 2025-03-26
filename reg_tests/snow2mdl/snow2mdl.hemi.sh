#!/bin/bash

#--------------------------------------------------------------------------
# Mimic v16 and prior GFS OPS, which used hemispheric afwa/airforce data.  
# This script is run from its machine-specific driver.
#--------------------------------------------------------------------------

echo "BEGIN SNOW2MDL HEMI TEST."

set -x

rm -fr $DATA
mkdir -p $DATA
cd $DATA

cat << EOF > ./fort.41
 &source_data
  autosnow_file=""
  nesdis_snow_file="$HOMEreg/input_data/hemi/imssnow96.grb"
  nesdis_lsmask_file=""
  afwa_snow_global_file=""
  afwa_snow_nh_file="$HOMEreg/input_data/hemi/NPR.SNWN.SP.S1200.MESH16"
  afwa_snow_sh_file="$HOMEreg/input_data/hemi/NPR.SNWS.SP.S1200.MESH16"
  afwa_lsmask_nh_file=""
  afwa_lsmask_sh_file=""
 /
 &qc
  climo_qc_file="$HOMEgfs/fix/am/emcsfc_snow_cover_climo.grib2"
 /
 &model_specs
  model_lat_file="$HOMEgfs/fix/am/global_latitudes.t1534.3072.1536.grb"
  model_lon_file="$HOMEgfs/fix/am/global_longitudes.t1534.3072.1536.grb"
  model_lsmask_file="$HOMEgfs/fix/am/global_slmask.t1534.3072.1536.grb"
  gfs_lpl_file="$HOMEgfs/fix/am/global_lonsperlat.t1534.3072.1536.txt"
  /
 &output_data
  model_snow_file="./snogrb_model"
  output_grib2=.false.
 /
 &output_grib_time
  grib_year=2012
  grib_month=10
  grib_day=29
  grib_hour=0
 /
 &parameters
  lat_threshold=55.0
  min_snow_depth=0.05
  snow_cvr_threshold=50.0
 /
EOF

eval $HOMEgfs/exec/emcsfc_snow2mdl >> OUTPUT 2> errfile
iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< SNOW2MDL HEMI TEST FAILED. <<<"
  exit $iret
fi

test_failed=0

cmp ${DATA}/snogrb_model $HOMEreg/baseline_data/t1534.hemi/snogrb_model
iret=$?
if [ $iret -ne 0 ]; then
  test_failed=1
fi

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL HEMI TEST FAILED. >>>"
  echo "*********************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    cd $DATA
    $HOMEgfs/reg_tests/update_baseline.sh $HOMEreg "t1534.hemi" $commit_num
  fi
else
  echo
  echo "*********************************"
  echo "<<< SNOW2MDL HEMI TEST PASSED. >>>"
  echo "*********************************"
fi

exit
