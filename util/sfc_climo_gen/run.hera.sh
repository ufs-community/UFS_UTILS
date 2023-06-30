#!/bin/bash

#------------------------------------------------------------
# Run the sfc_climo_gen program stand-alone on Hera using
# pre-exiting 'grid' and 'orography files.
#
# To run, type: 'sbatch $script'
#------------------------------------------------------------

#SBATCH -J sfc_climo_gen
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --partition=bigmem
#SBATCH -q debug
#SBATCH -t 00:30:00

set -x

export APRUN_SFC="srun"

export BASE_DIR=$SLURM_SUBMIT_DIR/../..

source ${BASE_DIR}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${BASE_DIR}/modulefiles
module load build.$target.intel
module list

#-------------------------------------
# Set model resolution.
#-------------------------------------

export res=3342
#export res=384.mx025

#-------------------------------------
# Where the model "grid", "mosaic" and "oro" files reside.
#-------------------------------------

#export FIX_FV3=${BASE_DIR}/fix/orog/C${res}

export FIX_FV3=/scratch1/NCEPDEV/global/glopara/fix/orog/20220805/C${res}
#-------------------------------------
# Uncomment for regional grids.
#-------------------------------------

#HALO=3
#export GRIDTYPE=regional

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
# For STATSGO soil type data
#   1) "statsgo.0.05" for global 0.05-deg data
#   2) "statsgo.0.03" for global 0.03-deg data
#   3) "statsgo.conus.30s" for CONUS 30s data
#   4) "statsgo.nh.30s" for NH 30s data
#   5) "statsgo.30s" for global 30s data
#
# For Beijing Norm. Univ. soil type data
#   1) "bnu.30s" for global 30s data
#-------------------------------------------------------------

#export veg_type_src="viirs.igbp.30s"

#export soil_type_src="statsgo.0.05"

export veg_type_src="viirs.igbp.0.05"
export soil_type_src="statsgo.0.05"


#
#soil type has been hard coded tyo use BNU V2
#


#-------------------------------------
# Set working directory and directory where output files will be saved.
#-------------------------------------

export WORK_DIR=/scratch2/NCEPDEV/stmp1/$LOGNAME/work.sfc
export SAVE_DIR=/scratch2/NCEPDEV/stmp1/$LOGNAME/repo_test/sfc.C${res}

#-------------------------------------
# Should not have to touch anything below here.
#-------------------------------------

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

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
