#!/bin/bash

#-------------------------------------------------------------------------
# Run sfc_climo_gen program to create surface fixed fields,
# such as vegetation type.
#
# Stand-alone regional grids may be run with any number of
# tasks.  All other configurations must be run with a
# MULTIPLE OF SIX MPI TASKS. 
#
# Some variable definitions:
#
# BASE_DIR                      Location of your repository.
# input_sfc_climo_dir           Location of raw input surface climo data
# exec_dir                      Location of program executable
# FIX_DIR                       Location of 'grid' and 'orog' files
# GRIDTYPE                      Flag to invoke logic for global nests
#                               and regional grids.  Valid values are
#                               'nest' and 'regional'.
# HALO                          Number of halo row/cols to remove
#                               for regional grid.
# mosaic_file                   Path/name of mosaic file.
# res                           Resolution of cubed-sphere grid
# SAVE_DIR                      Directory where output is saved
# WORK_DIR                      Temporary working directory
# SOIL_TYPE_FILE                Path/name of input soil type data.
# VEG_TYPE_FILE                 Path/name of input vegetation type data.
# vegsoilt_frac                 When true, outputs dominant soil and
#                               vegetation type category and the 
#                               fractional value of each category.
#                               When false, outputs dominant category.
#-------------------------------------------------------------------------

set -eux

res=${res:-96}
WORK_DIR=${WORK_DIR:-/scratch3/NCEPDEV/stmp1/$LOGNAME/sfc_climo_gen.C${res}}
SAVE_DIR=${SAVE_DIR:-$WORK_DIR}
BASE_DIR=${BASE_DIR:?}
exec_dir=${exec_dir:-$BASE_DIR/exec}
GRIDTYPE=${GRIDTYPE:-NULL}
FIX_FV3=${FIX_FV3:-/scratch4/NCEPDEV/global/save/glopara/git/fv3gfs/fix/fix_fv3_gmted2010/C${res}}
input_sfc_climo_dir=${input_sfc_climo_dir:?}
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
#input_substrate_temperature_file="${input_sfc_climo_dir}/substrate_temperature.2.6x1.5.nc"
#input_substrate_temperature_file="${input_sfc_climo_dir}/substrate_temperature.gfs.0.5.nc"
# update the above depending on the fix files
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
#the_orog_files='"C'${res}'_oro_data.tile1.nc","C'${res}'_oro_data.tile2.nc","C'${res}'_oro_data.tile3.nc","C'${res}'_oro_data.tile4.nc","C'${res}'_oro_data.tile5.nc","C'${res}'_oro_data.tile6.nc"'

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
