#!/bin/bash

#---------------------------------------------------------------------------
# Run chgres using v16 netcdf history data as input.  These history
# files are part of the OPS v16 gfs/gdas/enkf tarballs, and the
# v16 retro parallel gfs tarballs.  To run using the v16 retro
# gdas tarballs (which contain warm restart files), the 
# run_v16retro.chgres.sh is used.
#---------------------------------------------------------------------------

set -x

MEMBER=$1

FIX_FV3=$UFS_DIR/fix
FIX_ORO=${FIX_FV3}/orog
FIX_AM=${FIX_FV3}/am

WORKDIR=${WORKDIR:-$OUTDIR/work.${MEMBER}}

if [ ${MEMBER} == 'gdas' ] || [ ${MEMBER} == 'gfs' ] ; then
  CTAR=${CRES_HIRES}
#---------------------------------------------------------------------------
# Some gfs tarballs from the v16 retro parallels dont have 'atmos'
# in their path.  Account for this.
#---------------------------------------------------------------------------
  INPUT_DATA_DIR="${EXTRACT_DIR}/${MEMBER}.${yy}${mm}${dd}/${hh}/atmos"
  if [ ! -d ${INPUT_DATA_DIR} ]; then
    INPUT_DATA_DIR="${EXTRACT_DIR}/${MEMBER}.${yy}${mm}${dd}/${hh}"
  fi
  ATMFILE="${MEMBER}.t${hh}z.atmanl.nc"
  SFCFILE="${MEMBER}.t${hh}z.sfcanl.nc"
else  
  date10=`$NDATE -6 $yy$mm$dd$hh`
  yy_d=$(echo $date10 | cut -c1-4)
  mm_d=$(echo $date10 | cut -c5-6)
  dd_d=$(echo $date10 | cut -c7-8)
  hh_d=$(echo $date10 | cut -c9-10)
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkfgdas.${yy_d}${mm_d}${dd_d}/${hh_d}/atmos/mem${MEMBER}"
  ATMFILE="gdas.t${hh_d}z.atmf006.nc"
  SFCFILE="gdas.t${hh_d}z.sfcf006.nc"
fi

rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

source $UFS_DIR/util/gdas_init/set_fixed_files.sh

cat << EOF > fort.41

&config
 fix_dir_target_grid="${FIX_ORO}/${ORO_DIR}/fix_sfc"
 mosaic_file_target_grid="${FIX_ORO}/${ORO_DIR}/${CTAR}_mosaic.nc"
 orog_dir_target_grid="${FIX_ORO}/${ORO_DIR}"
 orog_files_target_grid="${ORO_NAME}.tile1.nc","${ORO_NAME}.tile2.nc","${ORO_NAME}.tile3.nc","${ORO_NAME}.tile4.nc","${ORO_NAME}.tile5.nc","${ORO_NAME}.tile6.nc"
 data_dir_input_grid="${INPUT_DATA_DIR}"
 atm_files_input_grid="${ATMFILE}"
 sfc_files_input_grid="${SFCFILE}"
 vcoord_file_target_grid="${FIX_AM}/global_hyblev.l${LEVS}.txt"
 cycle_mon=$mm
 cycle_day=$dd
 cycle_hour=$hh
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 input_type="gaussian_netcdf"
 tracers="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
 tracers_input="spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"
/
EOF

$APRUN $UFS_DIR/exec/chgres_cube
rc=$?

if [ $rc != 0 ]; then
  exit $rc
fi

$UFS_DIR/util/gdas_init/copy_coldstart_files.sh $MEMBER $OUTDIR $yy $mm $dd $hh $INPUT_DATA_DIR

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER $MEMBER

exit 0
