#!/bin/sh

set -x

OUTDIR=$OUTDIR/c96_fv3_restart
rm -fr $OUTDIR
mkdir -p $OUTDIR
cd $OUTDIR

cat << EOF > ./fort.41
&config
 mosaic_file_target_grid="${HOMEreg}/fix/C96/C96_mosaic.nc"
 fix_dir_target_grid="${HOMEreg}/fix/C96/fix_sfc"
 orog_dir_target_grid="${HOMEreg}/fix/C96"
 orog_files_target_grid="C96_oro_data.tile1.nc","C96_oro_data.tile2.nc","C96_oro_data.tile3.nc","C96_oro_data.tile4.nc","C96_oro_data.tile5.nc","C96_oro_data.tile6.nc"
 vcoord_file_target_grid="${HOMEufs}/fix/fix_am/global_hyblev.l64.txt"
 mosaic_file_input_grid="${HOMEreg}/fix/C384/C384_mosaic.nc"
 orog_dir_input_grid="${HOMEreg}/fix/C384"
 orog_files_input_grid="C384_oro_data.tile1.nc","C384_oro_data.tile2.nc","C384_oro_data.tile3.nc","C384_oro_data.tile4.nc","C384_oro_data.tile5.nc","C384_oro_data.tile6.nc"
 data_dir_input_grid="${INPUT_DATA}"
 atm_core_files_input_grid="20190706.120000.fv_core.res.tile1.nc","20190706.120000.fv_core.res.tile2.nc","20190706.120000.fv_core.res.tile3.nc","20190706.120000.fv_core.res.tile4.nc","20190706.120000.fv_core.res.tile5.nc","20190706.120000.fv_core.res.tile6.nc","20190706.120000.fv_core.res.nc"
 atm_tracer_files_input_grid="20190706.120000.fv_tracer.res.tile1.nc","20190706.120000.fv_tracer.res.tile2.nc","20190706.120000.fv_tracer.res.tile3.nc","20190706.120000.fv_tracer.res.tile4.nc","20190706.120000.fv_tracer.res.tile5.nc","20190706.120000.fv_tracer.res.tile6.nc"
 sfc_files_input_grid="20190706.120000.sfc_data.tile1.nc","20190706.120000.sfc_data.tile2.nc","20190706.120000.sfc_data.tile3.nc","20190706.120000.sfc_data.tile4.nc","20190706.120000.sfc_data.tile5.nc","20190706.120000.sfc_data.tile6.nc"
 cycle_mon=7
 cycle_day=6
 cycle_hour=12
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 input_type="restart"
 tracers="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
 tracers_input="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
/

EOF

date

$APRUN ${HOMEufs}/exec/chgres_cube.exe

iret=$?
if [ $iret -ne 0 ]; then
  echo "<<< C96 FV3 RESTART TEST FAILED. <<<"
  exit $iret
fi

date

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_fv3_restart/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 FV3 RESTART TEST FAILED. >>>"
else
  echo "<<< C96 FV3 RESTART TEST PASSED. >>>"
fi

exit 0
