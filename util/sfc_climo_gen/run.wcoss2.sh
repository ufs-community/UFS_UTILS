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
#PBS -l select=1:ncpus=24:mem=75GB

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

export FIX_FV3=${BASE_DIR}/fix/fix_fv3_gmted2010/C${res}

#-------------------------------------
# Uncomment for regional grids.
#-------------------------------------

##HALO=3
##export GRIDTYPE=regional

#-------------------------------------
# Choose which virrs data to use.
#-------------------------------------

export veg_type_src="viirs.igbp.0.05"    # Use global 0.05-degree viirs data
#export veg_type_src="viirs.igbp.0.1"     # Use global 0.1-degree viirs data
#export veg_type_src="viirs.igbp.0.03"   # Use global 0.03-degree viirs data
#export veg_type_src="viirs.igbp.conus.0.01"   # Use CONUS 0.01-degree virrs data. Do not
                                               # use for global grids.

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

export input_sfc_climo_dir=${BASE_DIR}/fix/fix_sfc_climo

ulimit -a
ulimit -s unlimited

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
