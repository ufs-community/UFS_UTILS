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
  ATMFILE="${MEMBER}.t${hh}z.atmanl.nemsio"
  SFCFILE="${MEMBER}.t${hh}z.sfcanl.nemsio"
  NSTFILE="${MEMBER}.t${hh}z.nstanl.nemsio"
else  
  CTAR=${CRES_ENKF}
  INPUT_DATA_DIR="${EXTRACT_DIR}/enkf.${yy}${mm}${dd}/${hh}/mem${MEMBER}"
  ATMFILE="gdas.t${hh}z.ratmanl.mem${MEMBER}.nemsio"
  SFCFILE="gdas.t${hh}z.sfcanl.mem${MEMBER}.nemsio"
  NSTFILE="gdas.t${hh}z.nstanl.mem${MEMBER}.nemsio"
fi

rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

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

$UFS_DIR/util/gdas_init/copy_coldstart_files.sh $MEMBER $OUTDIR $yy $mm $dd $hh $INPUT_DATA_DIR

rm -fr $WORKDIR

set +x
echo CHGRES COMPLETED FOR MEMBER $MEMBER

exit 0
