#!/bin/bash

#----------------------------------------------------------------
# Run chgres using gfs v14 data as input.
#----------------------------------------------------------------

set -x

MEMBER=$1

FIX_FV3=$UFS_DIR/fix
FIX_ORO=${FIX_FV3}/fix_fv3_gmted2010
FIX_AM=${FIX_FV3}/fix_am

date10=$yy$mm$dd$hh
yy=$(echo $date10 | cut -c1-4)
mm=$(echo $date10 | cut -c5-6)
dd=$(echo $date10 | cut -c7-8)
hh=$(echo $date10 | cut -c9-10)

YMDH=${yy}${mm}${dd}.${hh}0000

WORKDIR=$OUTDIR/work.$MEMBER

if [ ${MEMBER} == 'hires' ]; then
  CTAR=${CRES_HIRES}
  INPUT_DATA_DIR="${EXTRACT_DIR}/gdas.${yy}${mm}${dd}/${hh}"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/gdas.${yy}${mm}${dd}/${hh}"
  OUTDIR=$OUTDIR/gdas.${yy}${mm}${dd}/${hh}/atmos
  ATMFILE="gdas.t${hh}z.atmanl.nemsio"
  SFCFILE="gdas.t${hh}z.sfcanl.nemsio"
  NSTFILE="gdas.t${hh}z.nstanl.nemsio"
else  
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  OUTDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/atmos/mem${MEMBER}
  ATMFILE="gdas.t${hh}z.ratmanl.mem${MEMBER}.nemsio"
  SFCFILE="gdas.t${hh}z.sfcanl.mem${MEMBER}.nemsio"
  NSTFILE="gdas.t${hh}z.nstanl.mem${MEMBER}.nemsio"
fi

rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

rm -fr $OUTDIR
mkdir -p $OUTDIR
mkdir -p $OUTDIR/INPUT

cat << EOF > fort.41

&config
 fix_dir_target_grid="${FIX_ORO}/${CTAR}/fix_sfc"
 mosaic_file_target_grid="${FIX_ORO}/${CTAR}/${CTAR}_mosaic.nc"
 orog_dir_target_grid="${FIX_ORO}/${CTAR}"
 orog_files_target_grid="${CTAR}_oro_data.tile1.nc","${CTAR}_oro_data.tile2.nc","${CTAR}_oro_data.tile3.nc","${CTAR}_oro_data.tile4.nc","${CTAR}_oro_data.tile5.nc","${CTAR}_oro_data.tile6.nc"
 data_dir_input_grid="${INPUT_DATA_DIR}"
 atm_files_input_grid="$ATMFILE"
 sfc_files_input_grid="$SFCFILE"
 nst_files_input_grid="$NSTFILE"
 vcoord_file_target_grid="${FIX_AM}/global_hyblev.l${LEVS}.txt"
 cycle_mon=$mm
 cycle_day=$dd
 cycle_hour=$hh
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 input_type="gfs_gaussian_nemsio"
 tracers="sphum","liq_wat","o3mr"
 tracers_input="spfh","clwmr","o3mr"
/
EOF

$APRUN $UFS_DIR/exec/chgres_cube
rc=$?

if [ $rc != 0 ]; then
  exit $rc
fi

mv gfs_ctrl.nc ${OUTDIR}/INPUT

for tile in 'tile1' 'tile2' 'tile3' 'tile4' 'tile5' 'tile6'
do
  mv out.atm.${tile}.nc  ${OUTDIR}/INPUT/gfs_data.${tile}.nc
  mv out.sfc.${tile}.nc  ${OUTDIR}/INPUT/sfc_data.${tile}.nc 
done

if [ ${MEMBER} == 'hires' ]; then
  cp ${RADSTAT_DATA_DIR}/*radstat* $OUTDIR
  cp ${RADSTAT_DATA_DIR}/*abias* $OUTDIR
  touch $OUTDIR/gdas.t${hh}z.loginc.txt
else
  touch $OUTDIR/enkfgdas.t${hh}z.loginc.txt
fi

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER $MEMBER

exit 0
