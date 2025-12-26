#!/bin/bash
################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         sfc_climo_gen.sh
# Script description:  Create surface climatology fixed fields
#
# Abstract: This script runs the sfc_climo_gen program to create surface
#           climatology fixed fields such as vegetation type, soil type,
#           substrate temperature, albedo, and other surface properties
#           needed for model initialization.
#
# Usage:  sfc_climo_gen.sh
#
#   Required Shell Variables:
#     BASE_DIR          Location of UFS_UTILS repository. REQUIRED.
#     input_sfc_climo_dir  Location of raw input surface climatology data.
#                          REQUIRED.
#     FIX_FV3           Location of grid and orog fixed files directory. REQUIRED.
#
#   Optional Shell Variables:
#     Grid configuration:
#     res               Resolution of cubed-sphere grid (e.g., 96, 384, 768).
#                       Defaults to 96.
#     GRIDTYPE          Grid type flag. Valid values: 'nest', 'regional', or NULL.
#                       Defaults to NULL (global uniform).
#     ocn               Ocean grid resolution (e.g., 025, 050, 100).
#                       When declared, uses orog files for coupled model.
#                       Defaults to undefined (atmosphere-only).
#     HALO              Number of halo rows/cols for regional grids.
#                       Defaults to 0.
#
#     Directory paths:
#     WORK_DIR          Temporary working directory.
#                       Defaults to /scratch3/NCEPDEV/stmp1/$LOGNAME/sfc_climo_gen.C${res}
#     SAVE_DIR          Directory where output is saved.
#                       Defaults to $WORK_DIR.
#     exec_dir          Location of sfc_climo_gen executable.
#                       Defaults to $BASE_DIR/exec
#
#     Input file configuration:
#     mosaic_file       Path/name of mosaic file.
#                       Defaults to $FIX_FV3/C${res}_mosaic.nc
#     VEG_TYPE_FILE     Path/name of input vegetation type data.
#                       Defaults to ${input_sfc_climo_dir}/vegetation_type.${veg_type_src}.nc
#     SOIL_TYPE_FILE    Path/name of input soil type data.
#                       Defaults to ${input_sfc_climo_dir}/soil_type.${soil_type_src}.nc
#     veg_type_src      Vegetation type source identifier.
#                       Defaults to "modis.igbp.0.05"
#     soil_type_src     Soil type source identifier.
#                       Defaults to "statsgo.0.05"
#
#     Processing options:
#     vegsoilt_frac     When .true., outputs dominant soil/vegetation category AND
#                       fractional values of each category. When .false., outputs
#                       dominant category only. Defaults to .false.
#
#     Execution:
#     APRUN_SFC         Command to run executable with MPI.
#                       Defaults to "aprun -j 1 -n 6 -N 6"
#                       NOTE: Must use task count that is MULTIPLE OF 6 (ESMF requirement).
#
#   Example:
#     export BASE_DIR=/path/to/UFS_UTILS
#     export input_sfc_climo_dir=/path/to/fix/sfc_climo
#     export FIX_FV3=/path/to/fix/orog/C96
#     export res=96
#     ./sfc_climo_gen.sh
#
# Remarks:
#   - Stand-alone regional grids may run with any number of MPI tasks
#   - All other configurations MUST run with task count that is MULTIPLE OF 6
#     (This is an ESMF library requirement for the sfc_climo_gen executable)
#   - Large grids may require tasks spread across multiple nodes
#   - For regional grids, script creates both halo and no-halo output versions
#   - Output files are named: C${res}.${field_name}.tileX.nc (global)
#     or C${res}.${field_name}.haloX.nc (regional)
#
# Attributes:
#   Language: POSIX shell
#
################################################################################

set -eux

res=${res:-96}
WORK_DIR=${WORK_DIR:-/scratch3/NCEPDEV/stmp1/$LOGNAME/sfc_climo_gen.C${res}}
SAVE_DIR=${SAVE_DIR:-$WORK_DIR}
BASE_DIR=${BASE_DIR:?"ERROR: BASE_DIR is required. Location of UFS_UTILS repository"}
exec_dir=${exec_dir:-$BASE_DIR/exec}
GRIDTYPE=${GRIDTYPE:-NULL}
FIX_FV3=${FIX_FV3:?"ERROR: FIX_FV3 is required. Location of grid/orog fixed files directory"}
input_sfc_climo_dir=${input_sfc_climo_dir:?"ERROR: input_sfc_climo_dir is required. Location of raw surface climatology data"}
mosaic_file=${mosaic_file:-$FIX_FV3/C${res}_mosaic.nc}
HALO=${HALO:-0}
vegsoilt_frac=${vegsoilt_frac:-.false.}
veg_type_src=${veg_type_src:-"modis.igbp.0.05"}
VEG_TYPE_FILE=${VEG_TYPE_FILE:-${input_sfc_climo_dir}/vegetation_type.${veg_type_src}.nc}
soil_type_src=${soil_type_src:-"statsgo.0.05"}
SOIL_TYPE_FILE=${SOIL_TYPE_FILE:-${input_sfc_climo_dir}/soil_type.${soil_type_src}.nc}

if [ ! -d $SAVE_DIR ]; then
  mkdir -p $SAVE_DIR
fi

rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

#----------------------------------------------------------------------------------
# The stand-alone regional and global nest are assumed to be tile 7.
#----------------------------------------------------------------------------------

if [[ $GRIDTYPE == "nest" ]] || [[ $GRIDTYPE == "regional" ]] ; then
  the_orog_files='"C'${res}'_oro_data.tile7.nc"'
else
  if declare -p ocn &>/dev/null;then	
     the_orog_files='"C'${res}.mx${ocn}'_oro_data.tile1.nc","C'${res}.mx${ocn}'_oro_data.tile2.nc","C'${res}.mx${ocn}'_oro_data.tile3.nc","C'${res}.mx${ocn}'_oro_data.tile4.nc","C'${res}.mx${ocn}'_oro_data.tile5.nc","C'${res}.mx${ocn}'_oro_data.tile6.nc"'
  else
     the_orog_files='"C'${res}'_oro_data.tile1.nc","C'${res}'_oro_data.tile2.nc","C'${res}'_oro_data.tile3.nc","C'${res}'_oro_data.tile4.nc","C'${res}'_oro_data.tile5.nc","C'${res}'_oro_data.tile6.nc"'
  fi
fi

cat << EOF > ./fort.41
&config
input_facsf_file="${input_sfc_climo_dir}/facsf.1.0.nc"
input_substrate_temperature_file="${input_sfc_climo_dir}/substrate_temperature.gfs.0.5.nc"
input_maximum_snow_albedo_file="${input_sfc_climo_dir}/maximum_snow_albedo.0.05.nc"
input_snowfree_albedo_file="${input_sfc_climo_dir}/snowfree_albedo.4comp.0.05.nc"
input_slope_type_file="${input_sfc_climo_dir}/slope_type.1.0.nc"
input_soil_type_file="${SOIL_TYPE_FILE}"
input_soil_color_file="${input_sfc_climo_dir}/soil_color.clm.0.05.nc"
input_vegetation_type_file="${VEG_TYPE_FILE}"
input_vegetation_greenness_file="${input_sfc_climo_dir}/vegetation_greenness.0.144.nc"
mosaic_file_mdl="$mosaic_file"
orog_dir_mdl="$FIX_FV3"
orog_files_mdl=$the_orog_files
halo=$HALO
maximum_snow_albedo_method="bilinear"
snowfree_albedo_method="bilinear"
vegetation_greenness_method="bilinear"
fract_vegsoil_type=${vegsoilt_frac}
/
EOF


APRUN_SFC=${APRUN_SFC:-"aprun -j 1 -n 6 -N 6"}
$APRUN_SFC $exec_dir/sfc_climo_gen

rc=$?

if [[ $rc == 0 ]]; then
  if [[ $GRIDTYPE != "regional" ]]; then
    for files in *.nc
    do
      if [[ -f $files ]]; then
	if declare -p ocn &>/dev/null; then
        	mv $files ${SAVE_DIR}/C${res}.mx${ocn}.${files}
	else
		mv $files ${SAVE_DIR}/C${res}.${files}
	fi
      fi
    done
  else
    for files in *.halo.nc
    do
      if [[ -f $files ]]; then
        file2=${files%.halo.nc}
        mv $files ${SAVE_DIR}/C${res}.${file2}.halo${HALO}.nc
      fi
    done
    for files in *.nc
    do
      if [[ -f $files ]]; then
        file2=${files%.nc}
        mv $files ${SAVE_DIR}/C${res}.${file2}.halo0.nc
      fi
    done
  fi  # is regional?
else
  exit $rc
fi

exit 0
