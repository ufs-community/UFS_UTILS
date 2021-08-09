#!/bin/bash

#------------------------------------------------------------
# Run the sfc_climo_gen program stand-alone on Dell using
# pre-exiting 'grid' and 'orography files.
#------------------------------------------------------------

#BSUB -oo log
#BSUB -eo log
#BSUB -q debug
#BSUB -P GFS-DEV
#BSUB -J grid_fv3
#BSUB -W 0:10
#BSUB -x                 # run not shared
#BSUB -n 24              # total tasks
#BSUB -R span[ptile=24]   # tasks per node
#BSUB -R affinity[core(1):distribute=balance]

set -x

export BASE_DIR=$LS_SUBCWD/../..

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

export WORK_DIR=/gpfs/dell1/stmp/$LOGNAME/work.sfc
export SAVE_DIR=/gpfs/dell1/stmp/$LOGNAME/sfc.C${res}

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
export APRUN_SFC="mpirun -l"

ulimit -a
ulimit -s unlimited

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
