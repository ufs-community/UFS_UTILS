#!/bin/bash

#----------------------------------------------------------------------------
# Script name:  chgres_cube.sh
#
# Abstract:  Run the chgres program to initialize an FV3 run.
#
# See comments for variable definitions and setup information.
#----------------------------------------------------------------------------

set -x

#----------------------------------------------------------------------------
# Resolution of target grid.
#----------------------------------------------------------------------------

CRES=${CRES:-96}

#----------------------------------------------------------------------------
# Set up environment paths.
#
# HOMEufs - Location of root ufs_utils directory.
# EXECufs - Location of ufs_utils executable directory.
# FIXufs  - Location of ufs_utils root fixed data directory.
# FIXfv3  - Location of target grid orography and 'grid' files.
# FIXsfc  - Location of target grid surface climatological files.
# FIXam   - Location of vertical coordinate definition file for target grid.
#----------------------------------------------------------------------------

ufs_ver=${ufs_ver:-v1.0.0}
envir=${envir:-prod}
NWROOT=${NWROOT:-/nw${envir}}
HOMEufs=${HOMEufs:-${NWROOT}/ufs_util.${ufs_ver}}
EXECufs=${EXECufs:-$HOMEufs/exec}
FIXufs=${FIXufs:-$HOMEufs/fix}
FIXfv3=${FIXfv3:-$FIXufs/fix_fv3_gmted2010/C${CRES}}
FIXsfc=${FIXsfc:-$FIXfv3/fix_sfc}
FIXam=${FIXam:-$FIXufs/fix_am}

#----------------------------------------------------------------------------
# CDATE - YYMMDDHH of your run.
#----------------------------------------------------------------------------

CDATE=${CDATE:?}
iy=$(echo $CDATE|cut -c1-4)
im=$(echo $CDATE|cut -c5-6)
id=$(echo $CDATE|cut -c7-8)
ih=$(echo $CDATE|cut -c9-10)

#----------------------------------------------------------------------------
# Variables for stand-alone regional grids.
#
# REGIONAL - Set to 1 to create remove halo and create lateral boundary
#            file.  Set to 2 for lateral boundary file only.  Set to
#            0 for non-regional grids.
# HALO_BNDY  - Number of rows/cols for lateral boundaries.
# HALO_BLEND - Number of rows/cols for blending zone.
#----------------------------------------------------------------------------

REGIONAL=${REGIONAL:-0}
HALO_BNDY=${HALO_BNDY:-0}
HALO_BLEND=${HALO_BLEND:-0}

#----------------------------------------------------------------------------
# INPUT_TYPE - Input data type:  
#        'restart' for tiled fv3 warm restart files.  
#        'history' for tiled fv3 history files.
#        'gaussian_nemsio' for fv3 gaussian nemsio files.
#        'gaussian_netcdf' for fv3 gaussian netcdf files.
#        'grib2' for fv3gfs grib2 files.
#        'gfs_gaussain_nemsio' for spectral gfs nemsio files.
#        'gfs_sigio' for spectral gfs sigio/sfcio files.
#
# MOSAIC_FILE_INPUT_GRID - Path/Name of mosaic file for input grid.  Only
#                          used for 'history' and 'restart' INPUT_TYPE.
#                          Set to NULL otherwise.
#
# OROG_DIR_INPUT_GRID - Location of orography and grid files for input grid.
#                       Only used for 'history' and 'restart' INPUT_TYPE.
#                       Set to NULL otherwise.
#
# OROG_FILES_INPUT_GRID - List of orography files for input grid.  Only
#                         used for 'history' and 'restart' INPUT_TYPE.
#                         Set to NULL otherwise.
#----------------------------------------------------------------------------

INPUT_TYPE=${INPUT_TYPE:-"gaussian_nemsio"}
MOSAIC_FILE_INPUT_GRID=${MOSAIC_FILE_INPUT_GRID:-NULL}
OROG_DIR_INPUT_GRID=${OROG_DIR_INPUT_GRID:-NULL}
OROG_FILES_INPUT_GRID=${OROG_FILES_INPUT_GRID:-NULL}

#----------------------------------------------------------------------------
# COMIN       - Location of input data
# CONVERT_ATM - Convert atmospheric fields when true
# CONVERT_SFC - Convert surface fields when true
# CONVERT_NST - Convert nst fields when true
#----------------------------------------------------------------------------

CONVERT_ATM=${CONVERT_ATM:-.true.}
CONVERT_SFC=${CONVERT_SFC:-.true.}
CONVERT_NST=${CONVERT_NST:-.true.}

COMIN=${COMIN:-$PWD}

#----------------------------------------------------------------------------
# ATM_FILES_INPUT - Input atmospheric data file(s).  Not used for 'restart'
#                   or 'grib2' INPUT_TYPE.
#
# ATM_CORE_FILES - Input atmospheric core files.  Used for 'restart' 
#                  INPUT_TYPE only.  The first six entries are the tiled
#                  files.  The seventh is the file containing the 
#                  vertical coord definition.
#
# ATM_TRACER_FILES_INPUT - Input atmospheric tracer files for each tile.
#                          Used for 'restart' INPUT_TYPE only.
#
# SFC_FILES_INPUT - Input surface data file(s).  Not used for 'grib2'
#                   INPUT_TYPE.
#
# NST_FILES_INPUT - Input nst data file.  'gfs_gaussian_nemsio' INPUT_TYPE only.
#
# GRIB2_FILE_INPUT - Input gfs grib2 data file.  Only used for 'grib2'
#                    INPUT_TYPE.
#
# TRACERS_INPUT - List of input atmospheric tracer records to be processed.
#                 Not used for 'grib2' INPUT_TYPE.
#----------------------------------------------------------------------------

ATM_FILES_INPUT=${ATM_FILES_INPUT:-NULL}
ATM_CORE_FILES_INPUT=${ATM_CORE_FILES_INPUT:-NULL}
ATM_TRACER_FILES_INPUT=${ATM_TRACER_FILES_INPUT:-NULL}
SFC_FILES_INPUT=${SFC_FILES_INPUT:-NULL}
NST_FILES_INPUT=${NST_FILES_INPUT:-NULL}
GRIB2_FILE_INPUT=${GRIB2_FILE_INPUT:-NULL}
TRACERS_INPUT=${TRACERS_INPUT:-'"spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"'}

#----------------------------------------------------------------------------
#
# VARMAP_FILE - Variable mapping table.  Only used for 'grib2' INPUT_TYPE.
#
# TRACERS_TARGET - List of target tracer records. Must corresponde with
#                  with TRACERS_INPUT.  Not used for 'grib2' INPUT_TYPE.
#
# VCOORD_FILE - File containing vertical coordinate definition for target
#               grid.
#
# MOSAIC FILE_TARGET_GRID - Mosaic file for target grid (include path).
#                           The associated 'grid' files assumed to be in
#                           FIXfv3.
#
# OROG_FILES_TARGET_GRID - Orography file(s) for target grid.  Assumed to
#                          be located in FIXfv3.
#----------------------------------------------------------------------------

VARMAP_FILE=${VARMAP_FILE:-NULL}

TRACERS_TARGET=${TRACERS_TARGET:-'"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'}

VCOORD_FILE=${VCOORD_FILE:-${FIXam}/global_hyblev.l65.txt}

MOSAIC_FILE_TARGET_GRID=${MOSAIC_FILE_TARGET_GRID:-${FIXfv3}/C${CRES}_mosaic.nc}

OROG_FILES_TARGET_GRID=${OROG_FILES_TARGET_GRID:-NULL}
if [ $OROG_FILES_TARGET_GRID == NULL ]; then
  OROG_FILES_TARGET_GRID='C'${CRES}'_oro_data.tile1.nc","C'${CRES}'_oro_data.tile2.nc"'
  OROG_FILES_TARGET_GRID=${OROG_FILES_TARGET_GRID}',"C'${CRES}'_oro_data.tile3.nc","C'${CRES}'_oro_data.tile4.nc"'
  OROG_FILES_TARGET_GRID=${OROG_FILES_TARGET_GRID}',"C'${CRES}'_oro_data.tile5.nc","C'${CRES}'_oro_data.tile6.nc'
fi

#----------------------------------------------------------------------------
# APRUN - machine specific command to run program.
# CHGRESEXEC - program executable.
# OMP_NUM_THREADS - threads most useful for 'gfs_sigio' INPUT_TYPE.
# DATA - working directory.
# PGMOUT - standard output file
# PGMERR - standard error file
# REDOUT - standard output redirect
# REDERR - standard error redirect
#----------------------------------------------------------------------------

APRUN=${APRUN:-time}
CHGRESEXEC=${CHGRESEXEC:-${EXECufs}/chgres_cube}

export OMP_NUM_THREADS=${OMP_NUM_THREADS_CH:-1}

PGMOUT=${PGMOUT:-${pgmout:-'&1'}}
PGMERR=${PGMERR:-${pgmerr:-'&2'}}
REDOUT=${REDOUT:-'1>'}
REDERR=${REDERR:-'2>'}

DATA=${DATA:-$PWD/chgres}
mkdir -p $DATA
cd $DATA || exit 99

rm -f ./fort.41

cat << EOF > ./fort.41
 &config
  mosaic_file_target_grid="${MOSAIC_FILE_TARGET_GRID}"
  fix_dir_target_grid="${FIXsfc}"
  orog_dir_target_grid="${FIXfv3}"
  orog_files_target_grid="${OROG_FILES_TARGET_GRID}"
  vcoord_file_target_grid="${VCOORD_FILE}"
  mosaic_file_input_grid="${MOSAIC_FILE_INPUT_GRID}"
  orog_dir_input_grid="${OROG_DIR_INPUT_GRID}"
  orog_files_input_grid="${OROG_FILES_INPUT_GRID}"
  data_dir_input_grid="${COMIN}"
  atm_files_input_grid="${ATM_FILES_INPUT}"
  atm_core_files_input_grid="${ATM_CORE_FILES_INPUT}"
  atm_tracer_files_input_grid="${ATM_TRACER_FILES_INPUT}"
  sfc_files_input_grid="${SFC_FILES_INPUT}"
  nst_files_input_grid="${NST_FILES_INPUT}"
  grib2_file_input_grid="${GRIB2_FILE_INPUT}"
  varmap_file="${VARMAP_FILE}"
  cycle_mon=$im
  cycle_day=$id
  cycle_hour=$ih
  convert_atm=$CONVERT_ATM
  convert_sfc=$CONVERT_SFC
  convert_nst=$CONVERT_NST
  input_type="${INPUT_TYPE}"
  tracers=$TRACERS_TARGET
  tracers_input=$TRACERS_INPUT
  regional=$REGIONAL
  halo_bndy=$HALO_BNDY
  halo_blend=$HALO_BLEND
 /
EOF

$APRUN $CHGRESEXEC $REDOUT$PGMOUT $REDERR$PGMERR

iret=$?
if [ $iret -ne 0 ]; then
  echo "FATAL ERROR RUNNING CHGRES"
  exit $iret
fi

exit
