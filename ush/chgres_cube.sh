#!/bin/bash

#----------------------------------------------------------------------------
# Script name:  chgres_cube.sh
#
# Abstract:  Run the chgres program to initialize an FV3 run.
#----------------------------------------------------------------------------

CDATE=${CDATE:?}

iy=$(echo $CDATE|cut -c1-4)
im=$(echo $CDATE|cut -c5-6)
id=$(echo $CDATE|cut -c7-8)
ih=$(echo $CDATE|cut -c9-10)

CRES=${CRES:-96}

REGIONAL=${REGIONAL:-0}
HALO_BNDY=${HALO_BNDY:-0}
HALO_BLEND=${HALO_BLEND:-0}

HOMEufs=${HOMEufs:-/nwprod2/ufs_util}
EXECufs=${EXECufs:-$HOMEufs/exec}
FIXfv3=${FIXfv3:-$HOMEufs/fix/fix_fv3_gmted2010/C${CRES}}
FIXsfc=${FIXsfc:-$FIXfv3/fix_sfc}
FIXam=${FIXam:-$HOMEufs/fix/fix_am}

# Input grid stuff

MOSAIC_FILE_INPUT_GRID=${MOSAIC_FILE_INPUT_GRID:-NULL}
OROG_DIR_INPUT_GRID=${OROG_DIR_INPUT_GRID:-NULL}
OROG_FILES_INPUT_GRID=${OROG_FILES_INPUT_GRID:-NULL}

# Input data stuff

INPUT_TYPE=${INPUT_TYPE:-"gaussian"}
CONVERT_ATM=${CONVERT_ATM:-.true.}
CONVERT_SFC=${CONVERT_SFC:-.true.}
CONVERT_NST=${CONVERT_NST:-.true.}
TRACERS_INPUT=${TRACERS_INPUT:-'"spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"'}
COMIN=${COMIN:-$PWD}
ATM_FILES_INPUT=${ATM_FILES_INPUT:-gfs.t${ih}z.atmf000.nemsio}
SFC_FILES_INPUT=${SFC_FILES_INPUT:-gfs.t${ih}z.sfcf000.nemsio}
NST_FILES_INPUT=${NST_FILES_INPUT:-NULL}

# Target grid stuff

TRACERS_OUTPUT=${TRACERS_OUTPUT:-'"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'}

VCOORD_FILE=${VCOORD_FILE:-${FIXam}/global_hyblev.l65.txt}

MOSAIC_FILE_TARGET_GRID=${MOSAIC_FILE_TARGET_GRID:-${FIXfv3}/C${CRES}_mosaic.nc}

OROG_FILES_TARGET_GRID=${OROG_FILES_TARGET_GRID:-NULL}
if [ $OROG_FILES_TARGET_GRID == NULL ]; then
  OROG_FILES_TARGET_GRID='"C'${CRES}'_oro_data.tile1.nc"'
  tile=2
  while [ $tile -le 6 ]
  do
    OROG_FILES_TARGET_GRID=${OROG_FILES_TARGET_GRID}',"C'${CRES}'_oro_data.tile'${tile}'.nc"'
    ((tile=tile+1)) 
  done
fi

APRUN=${APRUN:-time}
CHGRESEXEC=${CHGRESEXEC:-${EXECufs}/chgres_cube.exe}

export OMP_NUM_THREADS=${OMP_NUM_THREADS_CY:-1}

DATA=${DATA:-$PWD/chgres}
mkdir -p $DATA
cd $DATA || exit 99

rm -f ./fort.41

cat << EOF > ./fort.41
 &config
  mosaic_file_target_grid="${MOSAIC_FILE_TARGET_GRID}"
  fix_dir_target_grid="${FIXsfc}"
  orog_dir_target_grid="${FIXfv3}"
  orog_files_target_grid=$OROG_FILES_TARGET_GRID
  vcoord_file_target_grid="${VCOORD_FILE}"
  mosaic_file_input_grid="${MOSAIC_FILE_INPUT_GRID}"
  orog_dir_input_grid="${OROG_DIR_INPUT_GRID}"
  orog_files_input_grid="${OROG_FILES_INPUT_GRID}"
  data_dir_input_grid="${COMIN}"
  atm_files_input_grid="${ATM_FILES_INPUT}"
  sfc_files_input_grid="${SFC_FILES_INPUT}"
  nst_files_input_grid="${NST_FILES_INPUT}"
  cycle_mon=$im
  cycle_day=$id
  cycle_hour=$ih
  convert_atm=$CONVERT_ATM
  convert_sfc=$CONVERT_SFC
  convert_nst=$CONVERT_NST
  input_type="${INPUT_TYPE}"
  tracers=$TRACERS_OUTPUT
  tracers_input=$TRACERS_INPUT
  regional=$REGIONAL
  halo_bndy=$HALO_BNDY
  halo_blend=$HALO_BLEND
 /
EOF

$APRUN $CHGRESEXEC

iret=$?
if [ $iret -ne 0 ]; then
  echo "FATAL ERROR RUNNING CHGRES"
  exit $iret
fi

exit
