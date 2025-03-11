#!/bin/bash

#------------------------------------------------------------------
# Run regrid_sfc utility to regrid a gsi (Gaussian) increment 
# to the cube sphere grid.
# Compare output to a baseline set of files using the 'nccmp' 
# utility.
#------------------------------------------------------------------

set -x

NCCMP=${NCCMP:-$(which nccmp)}

REGRID_EXEC=${NWPROD}/exec/regridStates.x
FIXorog=${NWPROD}/fix/orog/

COMIN_REGTEST=$HOMEreg/input_data_noahmp

# resolutions
CASE_IN=C192
CASE_OUT=C192
OCNRES_OUT=050

LONB_CASE_IN=$((4*${CASE_IN:1}))
LATB_CASE_IN=$((2*${CASE_IN:1}))

# create clean data dir 
if [[ -d $DATA ]]; then
   rm -rf ${DATA}
fi

mkdir -p $DATA

cd ${DATA}

# input, increment
/bin/cp -p ${COMIN_REGTEST}/sfcincr_gsi ${DATA}/sfcincr_gsi

# input, fixed files
ln -sf "${FIXorog}/${CASE_IN}/gaussian.${LONB_CASE_IN}.${LATB_CASE_IN}.nc" \
        "${DATA}/gaussian_scrip.nc"

# output, fixed files
ln -sf "${FIXorog}/${CASE_OUT}/${CASE_OUT}_mosaic.nc" \
        "${DATA}/${CASE_OUT}_mosaic.nc"

ntiles=6
for n in $(seq 1 $ntiles); do
    ln -sf ${FIXorog}/${CASE_OUT}/sfc/${CASE_OUT}.mx${OCNRES_OUT}.vegetation_type.tile${n}.nc  ${DATA}/vegetation_type.tile${n}.nc
    ln -sf ${FIXorog}/${CASE_OUT}/${CASE_OUT}_grid.tile${n}.nc ${DATA}/${CASE_OUT}_grid.tile${n}.nc
done

# namelist
cat << EOF > regrid.nml
 &config
  n_vars=4,
  variable_list="soilt1_inc", "soilt2_inc", "slc1_inc", "slc2_inc",
  missing_value=0.,
  extrap_levs=2,
 /
 &input
  gridtype="gau_inc",
  ires=${LONB_CASE_IN},
  jres=${LATB_CASE_IN},
  fname="sfcincr_gsi",
  dir="./",
  fname_coord="gaussian_scrip.nc",
  dir_coord="./"
 /

 &output
  gridtype="fv3_rst",
  ires=${CASE_OUT:1},
  jres=${CASE_OUT:1},
  fname="sfci",
  dir="./",
  fname_mask="vegetation_type"
  dir_mask="./"
  dir_coord="./",
 /
EOF

# run the executable
$APRUN_REGRID $REGRID_EXEC

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< REGRID SFC TEST FAILED. >>>"
  exit $iret
fi

test_failed=0

# check the ouput
for files in sfci*
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/gauss2fv3incr/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo
  echo "**********************************************"
  echo "<<< REGRID SFC TEST FAILED. >>>"
  echo "**********************************************"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    ${NWPROD}/reg_tests/update_baseline.sh $HOMEreg "gauss2fv3incr" $commit_num
  fi
else
  echo
  echo "*****************************************"
  echo "<<< REGRID SFC TEST PASSED >>>"
  echo "*****************************************"
fi

exit
