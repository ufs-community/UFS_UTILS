#!/bin/bash

#-----------------------------------------------------------------------
#
# This script is run by the machine specific driver script.
#
# Set the following variables:
#
# res      - Grid resolution. Example: 384 or 384.mx025.
#
# FIX_FV3  - Location of the pre-existing 'grid' and 'orography'
#            files. Defaults to ${BASE_DIR}/fix/orog/C${res}, where
#            BASE_DIR is the location of the checked out repository.
#
#            The required files are:
#
#            'mosaic' file - C${res}_mosaic.nc (Note: 'res' without
#                                               the 'mx' extension.)
#
#            'grid' files - C${res}_grid.tile7.halo${HALO}.nc (regional grids).
#                           C${res}_grid.tile[1-6].nc (global grids).
#                           Note: 'res' without the 'mx' extension.
#          
#            'orog' files - C${res}_oro_data.tile7.halo${HALO}.nc (regional grids).
#                           C${res}_oro_data.tile[1-6].nc (global grids).
#
# GRIDTYPE - set to 'regional' for regional grids. Otherwise,
#            comment out.
#
# FIX_REG  - For regional grids. Hold links to the 'grid' and 'orog' files
#            with names expected by the program.
#
# HALO     - The number of halo rows/cols. Only for regional grids.
#            Otherwise, comment out.
#
# WORK_DIR - Working directory.
#
# SAVE_DIR - Directory where the surface files will be saved.
#
# veg_type_src - Input vegetation type data. Choices are:
#                  For viirs-based vegetation type data, set to:
#                  - "viirs.v3.igbp.30s" for global 30s data
#                  For the modis-based vegetation data, set to:
#                  - "modis.igbp.0.05" for global 0.05-deg data
#                  - "modis.igbp.0.03" for global 0.03-deg data
#                  - "modis.igbp.conus.30s" for CONUS 30s data
#                  - "modis.igbp.nh.30s" for NH 30s data
#                  - "modis.igbp.30s" for global 30s data
#
# soil_type_src - Input soil type data. Choices are:
#                   For Beijing Norm. Univ. soil type data
#                   - "bnu.v3.30s" for global 30s data
#                   - "statsgo.0.03" for global 0.03-deg data
#                   - "statsgo.conus.30s" for CONUS 30s data
#                   - "statsgo.nh.30s" for NH 30s data
#                   - "statsgo.30s" for global 30s data
#
# vegsoilt_frac - When .true., output the fraction of each
#                 vegetation and soil type and the dominant
#                 category. When .false., output dominant
#                 category only.
#-----------------------------------------------------------------------

set -x

export res=768.mx025

#HALO=4
#export GRIDTYPE=regional
#FIX_REG=/lfs/h2/emc/stmp/$LOGNAME/fix.reg


export veg_type_src="viirs.v3.igbp.30s"

export soil_type_src="bnu.v3.30s"

export WORK_DIR=/scratch1/NCEPDEV/stmp2/$LOGNAME/work.sfc
export SAVE_DIR=/scratch1/NCEPDEV/stmp2/$LOGNAME/sfc.C${res}

export FIX_FV3=${BASE_DIR}/fix/orog/C${res}

# Requires much more resources when true. On hera, use 6 nodes,
# 12 tasks per node. On WCOSS2, use 5 nodes, 12 tasks per node.

export vegsoilt_frac=.false.

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# Should not have to touch anything below here.
#------------------------------------------------------------------------
#------------------------------------------------------------------------

if [[ "$GRIDTYPE" = "regional" ]]; then
  mkdir -p $FIX_REG
  ln -fs $FIX_FV3/C${res}_grid.tile7.halo${HALO}.nc $FIX_REG/C${res}_grid.tile7.halo${HALO}.nc
  ln -fs $FIX_FV3/C${res}_oro_data.tile7.halo${HALO}.nc $FIX_REG/C${res}_oro_data.tile7.nc
  ln -fs $FIX_FV3/C${res}_mosaic.nc $FIX_REG/C${res}_mosaic.nc
  export mosaic_file=$FIX_REG/C${res}_mosaic.nc
  export FIX_FV3=$FIX_REG
  HALO=$(( $HALO + 1 ))
  export HALO
else
  res2=${res//".mx"*}
  export mosaic_file=$FIX_FV3/C${res2}_mosaic.nc
fi

export input_sfc_climo_dir=${BASE_DIR}/fix/sfc_climo

ulimit -a 
ulimit -s unlimited

source ${BASE_DIR}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${BASE_DIR}/modulefiles
module load build.$target.intel
module list

rm -fr $WORK_DIR $SAVE_DIR

export APRUN

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
