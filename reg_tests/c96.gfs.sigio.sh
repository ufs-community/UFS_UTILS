#!/bin/sh

set -x

OUTDIR=$OUTDIR/c96_gfs_sigio
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
 mosaic_file_input_grid="NULL"
 orog_dir_input_grid="NULL"
 orog_files_input_grid="NULL"
 data_dir_input_grid="${INPUT_DATA}"
 atm_files_input_grid="gdas.t00z.sanl"
 nst_files_input_grid="NULL"
 sfc_files_input_grid="gdas.t00z.sfcanl"
 cycle_mon=7
 cycle_day=17
 cycle_hour=0
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.false.
 input_type="gfs_spectral"
 tracers_input="spfh","o3mr","clwmr"
 tracers="sphum","o3mr","liq_wat"
/

EOF

date

$APRUN ${HOMEufs}/exec/chgres_cube.exe

iret=$?
if [ $iret -ne 0 ]; then
  echo "<<< C96 GFS SIGIO TEST FAILED. <<<"
  exit $iret
fi

date

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c96_gfs_sigio/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C96 GFS SIGIO TEST FAILED. >>>"
else
  echo "<<< C96 GFS SIGIO TEST PASSED. >>>"
fi

exit 0
