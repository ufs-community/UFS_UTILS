#!/bin/bash

#----------------------------------------------------------------
# Run chgres using v15 nemsio data as input.  This is used
# for initializing GFS free forecasts.
#----------------------------------------------------------------

set -x

FIX_FV3=$UFS_DIR/fix
FIX_ORO=${FIX_FV3}/fix_fv3_gmted2010
FIX_AM=${FIX_FV3}/fix_am

WORKDIR=${WORKDIR:-$OUTDIR/work.gfs}

CTAR=${CRES_HIRES}
INPUT_DATA_DIR="${EXTRACT_DIR}/gfs.${yy}${mm}${dd}/${hh}"
OUTDIR=$OUTDIR/gfs.${yy}${mm}${dd}/${hh}/atmos
ATMFILE="gfs.t${hh}z.atmanl.nemsio"
SFCFILE="gfs.t${hh}z.sfcanl.nemsio"

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
 convert_nst=.true.
 input_type="gaussian_nemsio"
 tracers="sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"
 tracers_input="spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"
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

touch $OUTDIR/gfs.t${hh}z.loginc.txt

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER gfs

exit 0
