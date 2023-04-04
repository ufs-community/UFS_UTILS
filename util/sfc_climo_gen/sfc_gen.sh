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
# HALO     - The number of halo rows/cols. Only for regional grids.
#            Otherwise, comment out.
#
# WORK_DIR - Working directory.
#
# SAVE_DIR - Directory where the surface files will be saved.
#
# veg_type_src - Input vegetation type data. Choices are:
#                  For viirs-based vegetation type data, set to:
#                  - "viirs.igbp.0.1" for global 0.10-deg data
#                  - "viirs.igbp.0.05" for global 0.05-deg data
#                  - "viirs.igbp.0.03" for global 0.03-deg data
#                  - "viirs.igbp.conus.30s" for CONUS 30s data
#                  - "viirs.igbp.nh.30s" for NH 30s data
#                  - "viirs.igbp.30s" for global 30s data
#                  For the modis-based vegetation data, set to:
#                  - "modis.igbp.0.05" for global 0.05-deg data
#                  - "modis.igbp.0.03" for global 0.03-deg data
#                  - "modis.igbp.conus.30s" for CONUS 30s data
#                  - "modis.igbp.nh.30s" for NH 30s data
#                  - "modis.igbp.30s" for global 30s data
#
# soil_type_src - Input soil type data. Choices are:
#                   For STATSGO soil type data
#                   - "statsgo.0.05" for global 0.05-deg data
#                   - "statsgo.0.03" for global 0.03-deg data
#                   - "statsgo.conus.30s" for CONUS 30s data
#                   - "statsgo.nh.30s" for NH 30s data
#                   - "statsgo.30s" for global 30s data
#                   For Beijing Norm. Univ. soil type data
#                   - "bnu.30s" for global 30s data
#-----------------------------------------------------------------------

#export res=384
export res=384.mx025

##HALO=3
##export GRIDTYPE=regional

export veg_type_src="modis.igbp.0.05"

export soil_type_src="statsgo.0.05"

export WORK_DIR=/lfs/h2/emc/stmp/$LOGNAME/work.sfc
export SAVE_DIR=/lfs/h2/emc/stmp/$LOGNAME/sfc.C${res}

export FIX_FV3=${BASE_DIR}/fix/orog/C${res}

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# Should not have to touch anything below here.
#------------------------------------------------------------------------
#------------------------------------------------------------------------

if [[ $GRIDTYPE = "regional" ]]; then
  HALO=$(( $HALO + 1 ))
  export HALO
  ln -fs $FIX_FV3/C${res}_grid.tile7.halo${HALO}.nc $FIX_FV3/C${res}_grid.tile7.nc
  ln -fs $FIX_FV3/C${res}_oro_data.tile7.halo${HALO}.nc $FIX_FV3/C${res}_oro_data.tile7.nc
fi

res2=${res//".mx"*}
export mosaic_file=$FIX_FV3/C${res2}_mosaic.nc

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
