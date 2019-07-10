#!/bin/sh

set -x

OUTDIR=$OUTDIR/c96_regional
rm -fr $OUTDIR
mkdir -p $OUTDIR
cd $OUTDIR

cat << EOF > ./fort.41
&config
 mosaic_file_target_grid="${HOMEreg}/fix/C96.regional/C96_mosaic.nc"
 fix_dir_target_grid="${HOMEreg}/fix/C96.regional/fix_sfc"
 orog_dir_target_grid="${HOMEreg}/fix/C96.regional"
 orog_files_target_grid="C96_oro_data.tile7.nc"
 vcoord_file_target_grid="${HOMEufs}/fix/fix_am/global_hyblev.l64.txt"
 mosaic_file_input_grid="NULL"
 orog_dir_input_grid="NULL"
 orog_files_input_grid="NULL"
 data_dir_input_grid="${INPUT_DATA}"
 atm_files_input_grid="gfs.t12z.atmf000.nemsio"
 sfc_files_input_grid="gfs.t12z.sfcf000.nemsio"
 cycle_mon=7
 cycle_day=4
 cycle_hour=12
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 input_type="gaussian"
 regional=1
 halo_bndy=4
 halo_blend=0
 tracers="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
 tracers_input="spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"
/

EOF

date

$APRUN ${HOMEufs}/exec/chgres_cube.exe

iret=$?
if [ $iret -ne 0 ]; then
  echo "<<< C96 REGIONAL TEST FAILED. <<<"
  exit $iret
fi

date

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_regional/$files
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
