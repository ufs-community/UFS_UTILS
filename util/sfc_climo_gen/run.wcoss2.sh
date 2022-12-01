#!/bin/bash

#------------------------------------------------------------
# Run the sfc_climo_gen program stand-alone on WCOSS2 using
# pre-exiting 'grid' and 'orography files.
#
# To run, type: 'qsub $script'
#------------------------------------------------------------

#PBS -o log
#PBS -e log
#PBS -q debug
#PBS -A GFS-DEV
#PBS -N grid_fv3
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=24:mem=250GB

set -x

# Adjust according to the PBS -l statement.
export APRUN_SFC="mpiexec -n 24 -ppn 24 -cpu-bind core"

export BASE_DIR=$PBS_O_WORKDIR/../..

source ${BASE_DIR}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${BASE_DIR}/modulefiles
module load build.$target.intel
module list

#-------------------------------------
# Set model resolution.
#-------------------------------------

export res=384

#-------------------------------------
# Where the model "grid", "mosaic" and "oro" files reside.
#-------------------------------------

export FIX_FV3=${BASE_DIR}/fix/orog/C${res}

#-------------------------------------
# Uncomment for regional grids.
#-------------------------------------

##HALO=3
##export GRIDTYPE=regional

#-------------------------------------------------------------
# Choose which soil type and vegetation type data to use.
#
# For viirs-based vegetation type data, set to:
#   1) "viirs.igbp.0.1" for global 0.10-deg data
#   2) "viirs.igbp.0.05" for global 0.05-deg data
#   3) "viirs.igbp.0.03" for global 0.03-deg data
#   4) "viirs.igbp.conus.30s" for CONUS 30s data
#   4) "viirs.igbp.nh.30s" for NH 30s data
#   4) "viirs.igbp.30s" for global 30s data
#
# For the modis-based vegetation data, set to:
#   1) "modis.igbp.0.05" for global 0.05-deg data
#   2) "modis.igbp.0.03" for global 0.03-deg data
#   3) "modis.igbp.conus.30s" for CONUS 30s data
#   4) "modis.igbp.nh.30s" for NH 30s data
#   5) "modis.igbp.30s" for global 30s data
#
# Soil type data
#   1) "statsgo.0.05" for global 0.05-deg data
#   2) "statsgo.0.03" for global 0.03-deg data
#   3) "statsgo.conus.30s" for CONUS 30s data
#   4) "statsgo.nh.30s" for NH 30s data
#   5) "statsgo.30s" for global 30s data
#-------------------------------------------------------------

export veg_type_src="modis.igbp.0.05"

export soil_type_src="statsgo.0.05"

export vegsoilt_frac='.true.'  # When true, outputs percent of each
                               # soil and veg type category and a 
                               # dominate category. When false, only
                               # outputs the dominate category. A
                               # Fortran logical, so include the dots.

#-------------------------------------
# Set working directory and directory where output files will be saved.
#-------------------------------------

export WORK_DIR=/lfs/h2/emc/stmp/$LOGNAME/work.sfc
export SAVE_DIR=/lfs/h2/emc/stmp/$LOGNAME/sfc.C${res}

#-------------------------------------
# Should not have to touch anything below here.
#-------------------------------------

if [[ $GRIDTYPE = "regional" ]]; then
  HALO=$(( $HALO + 1 ))
  export HALO
  ln -fs $FIX_FV3/C${res}_grid.tile7.halo${HALO}.nc $FIX_FV3/C${res}_grid.tile7.nc
  ln -fs $FIX_FV3/C${res}_oro_data.tile7.halo${HALO}.nc $FIX_FV3/C${res}_oro_data.tile7.nc
fi

export input_sfc_climo_dir=${BASE_DIR}/fix/sfc_climo

ulimit -a
ulimit -s unlimited

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
