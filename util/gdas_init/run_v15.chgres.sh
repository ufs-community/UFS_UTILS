#!/bin/bash

#---------------------------------------------------------------------------
# Run chgres using gfs v15 data as input.
#---------------------------------------------------------------------------

set -x

MEMBER=$1

FIX_FV3=$UFS_DIR/fix
FIX_ORO=${FIX_FV3}/orog
FIX_AM=${FIX_FV3}/am

date10=`$NDATE -6 $yy$mm$dd$hh`
yy_d=$(echo $date10 | cut -c1-4)
mm_d=$(echo $date10 | cut -c5-6)
dd_d=$(echo $date10 | cut -c7-8)
hh_d=$(echo $date10 | cut -c9-10)

YMDH=${yy}${mm}${dd}.${hh}0000

WORKDIR=${WORKDIR:-$OUTDIR/work.${MEMBER}}

if [ ${MEMBER} == 'gdas' ]; then
  CINP=${CINP:-"C768"}
  CTAR=${CRES_HIRES}
  INPUT_DATA_DIR="${EXTRACT_DIR}/gdas.${yy_d}${mm_d}${dd_d}/${hh_d}/RESTART"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/gdas.${yy}${mm}${dd}/${hh}"
else  
  CINP=${CINP:-"C384"}
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkfgdas.${yy_d}${mm_d}${dd_d}/${hh_d}/mem${MEMBER}/RESTART"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/enkfgdas.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
fi

rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

source $GDAS_INIT_DIR/set_fixed_files.sh

cat << EOF > fort.41

&config
 fix_dir_target_grid="${FIX_ORO}/${ORO_DIR}/sfc"
 mosaic_file_target_grid="${FIX_ORO}/${ORO_DIR}/${CTAR}_mosaic.nc"
 orog_dir_target_grid="${FIX_ORO}/${ORO_DIR}"
 orog_files_target_grid="${ORO_NAME}.tile1.nc","${ORO_NAME}.tile2.nc","${ORO_NAME}.tile3.nc","${ORO_NAME}.tile4.nc","${ORO_NAME}.tile5.nc","${ORO_NAME}.tile6.nc"
 mosaic_file_input_grid="${FIX_ORO_INPUT}/${CINP}/${CINP}_mosaic.nc"
 orog_dir_input_grid="${FIX_ORO_INPUT}/${CINP}"
 orog_files_input_grid="${CINP}_oro_data.tile1.nc","${CINP}_oro_data.tile2.nc","${CINP}_oro_data.tile3.nc","${CINP}_oro_data.tile4.nc","${CINP}_oro_data.tile5.nc","${CINP}_oro_data.tile6.nc"
 data_dir_input_grid="${INPUT_DATA_DIR}"
 atm_core_files_input_grid="${YMDH}.fv_core.res.tile1.nc","${YMDH}.fv_core.res.tile2.nc","${YMDH}.fv_core.res.tile3.nc","${YMDH}.fv_core.res.tile4.nc","${YMDH}.fv_core.res.tile5.nc","${YMDH}.fv_core.res.tile6.nc","${YMDH}.fv_core.res.nc"
 atm_tracer_files_input_grid="${YMDH}.fv_tracer.res.tile1.nc","${YMDH}.fv_tracer.res.tile2.nc","${YMDH}.fv_tracer.res.tile3.nc","${YMDH}.fv_tracer.res.tile4.nc","${YMDH}.fv_tracer.res.tile5.nc","${YMDH}.fv_tracer.res.tile6.nc"
 vcoord_file_target_grid="${FIX_AM}/global_hyblev.l${LEVS}.txt"
 sfc_files_input_grid="${YMDH}.sfc_data.tile1.nc","${YMDH}.sfc_data.tile2.nc","${YMDH}.sfc_data.tile3.nc","${YMDH}.sfc_data.tile4.nc","${YMDH}.sfc_data.tile5.nc","${YMDH}.sfc_data.tile6.nc"
 cycle_mon=$mm
 cycle_day=$dd
 cycle_hour=$hh
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 tracers="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
 tracers_input="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
/
EOF

$APRUN $EXEC_DIR/chgres_cube
rc=$?

if [ $rc != 0 ]; then
  exit $rc
fi

$GDAS_INIT_DIR/copy_coldstart_files.sh $MEMBER $OUTDIR $yy $mm $dd $hh $RADSTAT_DATA_DIR

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER $MEMBER

exit 0
