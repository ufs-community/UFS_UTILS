#!/bin/sh

set -x

# Threads useful when ingesting spectral gfs sigio files.
# Otherwise set to 1.
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1024M

OUTDIR=$OUTDIR/c192_fv3_history
rm -fr $OUTDIR
mkdir -p $OUTDIR
cd $OUTDIR

cat << EOF > ./fort.41
&config
 mosaic_file_target_grid="${HOMEreg}/fix/C192/C192_mosaic.nc"
 fix_dir_target_grid="${HOMEreg}/fix/C192/fix_sfc"
 orog_dir_target_grid="${HOMEreg}/fix/C192"
 orog_files_target_grid="C192_oro_data.tile1.nc","C192_oro_data.tile2.nc","C192_oro_data.tile3.nc","C192_oro_data.tile4.nc","C192_oro_data.tile5.nc","C192_oro_data.tile6.nc"
 vcoord_file_target_grid="${HOMEufs}/fix/fix_am/global_hyblev.l64.txt"
 mosaic_file_input_grid="${HOMEreg}/fix/C96/C96_mosaic.nc"
 orog_dir_input_grid="${HOMEreg}/fix/C96"
 orog_files_input_grid="C96_oro_data.tile1.nc","C96_oro_data.tile2.nc","C96_oro_data.tile3.nc","C96_oro_data.tile4.nc","C96_oro_data.tile5.nc","C96_oro_data.tile6.nc"
 data_dir_input_grid="${INPUT_DATA}"
 atm_files_input_grid="dynf000.tile1.nc","dynf000.tile2.nc","dynf000.tile3.nc","dynf000.tile4.nc","dynf000.tile5.nc","dynf000.tile6.nc"
 sfc_files_input_grid="phyf000.tile1.nc","phyf000.tile2.nc","phyf000.tile3.nc","phyf000.tile4.nc","phyf000.tile5.nc","phyf000.tile6.nc"
 cycle_mon=10
 cycle_day=3
 cycle_hour=0
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 input_type="history"
 tracers="sphum","liq_wat","o3mr"
 tracers_input="spfh","clwmr","o3mr"
/

EOF

date

$APRUN ${HOMEufs}/exec/chgres_cube.exe

iret=$?
if [ $iret -ne 0 ]; then
  echo "<<< C192 FV3 HISTORY TEST FAILED. <<<"
  exit $iret
fi

date

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    nccmp -dmfqS $files $HOMEreg/baseline_data/c192_fv3_history/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C192 FV3 HISTORY TEST FAILED. >>>"
else
  echo "<<< C192 FV3 HISTORY TEST PASSED. >>>"
fi

exit 0
