#!/bin/bash

#----------------------------------------------------------------
# Run chgres using pre-v14 gfs data (sigio/sfcio format
# from the spectral gfs).
#----------------------------------------------------------------

set -x

MEMBER=$1

FIX_FV3=$UFS_DIR/fix
FIX_ORO=${FIX_FV3}/orog
FIX_AM=${FIX_FV3}/am

WORKDIR=${WORKDIR:-$OUTDIR/work.${MEMBER}}

if [ "${MEMBER}" = "gdas" ] || [ "${MEMBER}" = "gfs" ]; then
  CTAR=${CRES_HIRES}
  INPUT_DATA_DIR="${EXTRACT_DIR}/${MEMBER}.${yy}${mm}${dd}/${hh}"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/${MEMBER}.${yy}${mm}${dd}/${hh}"
  OUTDIR=$OUTDIR/${MEMBER}.${yy}${mm}${dd}/${hh}/atmos
  if [ "${MEMBER}" = "gdas" ]; then
    ATMFILE="gdas1.t${hh}z.sanl"
    SFCFILE="gdas1.t${hh}z.sfcanl"
  else
    ATMFILE="gfs.t${hh}z.sanl"
    SFCFILE="gfs.t${hh}z.sfcanl"
  fi
else  
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  OUTDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/atmos/mem${MEMBER}
  ATMFILE="siganl_${yy}${mm}${dd}${hh}_mem${MEMBER}"
  SFCFILE="sfcanl_${yy}${mm}${dd}${hh}_mem${MEMBER}"
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
 vcoord_file_target_grid="${FIX_AM}/global_hyblev.l${LEVS}.txt"
 cycle_mon=$mm
 cycle_day=$dd
 cycle_hour=$hh
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.false.
 input_type="gfs_sigio"
 tracers_input="spfh","o3mr","clwmr"
 tracers="sphum","o3mr","liq_wat"
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

if [ "${MEMBER}" = "gdas" ]; then
  cp ${RADSTAT_DATA_DIR}/*radstat* $OUTDIR
  cp ${RADSTAT_DATA_DIR}/*abias* $OUTDIR
  touch $OUTDIR/gdas.t${hh}z.loginc.txt
elif [ "${MEMBER}" = "gfs" ]; then
  touch $OUTDIR/gfs.t${hh}z.loginc.txt
else
  touch $OUTDIR/enkfgdas.t${hh}z.loginc.txt
fi

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER $MEMBER

exit 0
