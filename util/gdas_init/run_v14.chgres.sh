#!/bin/bash

#----------------------------------------------------------------
# Run chgres using gfs v14 data as input.
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
  ATMFILE="${MEMBER}.t${hh}z.atmanl.nemsio"
  SFCFILE="${MEMBER}.t${hh}z.sfcanl.nemsio"
  NSTFILE="${MEMBER}.t${hh}z.nstanl.nemsio"
else  
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  RADSTAT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  OUTDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/mem${MEMBER}/atmos
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

source $UFS_DIR/util/gdas_init/set_fixed_files.sh

cat << EOF > fort.41

&config
 fix_dir_target_grid="${FIX_ORO}/${ORO_DIR}/fix_sfc"
 mosaic_file_target_grid="${FIX_ORO}/${ORO_DIR}/${CTAR}_mosaic.nc"
 orog_dir_target_grid="${FIX_ORO}/${ORO_DIR}"
 orog_files_target_grid="${ORO_NAME}.tile1.nc","${ORO_NAME}.tile2.nc","${ORO_NAME}.tile3.nc","${ORO_NAME}.tile4.nc","${ORO_NAME}.tile5.nc","${ORO_NAME}.tile6.nc"
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
