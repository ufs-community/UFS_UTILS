#!/bin/bash

# Run sfc_climo_gen stand-alone on Orion.

#SBATCH -J fv3_grid_driver
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH -q debug
#SBATCH -t 00:10:00

set -x

export BASE_DIR=$SLURM_SUBMIT_DIR/../../..

source ${BASE_DIR}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${BASE_DIR}/modulefiles
module load build.$target.intel
module list

export res=768

# where the 'grid', 'mosaic' and 'oro' files reside.
export FIX_FV3=${BASE_DIR}/fix/fix_fv3_gmted2010/C${res}
####export FIX_FV3=/work/noaa/stmp/ggayno/yihua

# uncomment for regional grid
##HALO=3
##export GRIDTYPE=regional

# use global 0.05-degree viirs data
export VEG_FILE=/work/noaa/da/ggayno/save/ufs_utils.git/fv3.vegt.new.tundra.netcdf/fix_sfc_climo/vegetation_type.viirs.igbp.0.05.nc
# use global 0.03-degree viirs data
##export VEG_FILE=/work/noaa/da/ggayno/save/ufs_utils.git/fv3.vegt.new.tundra.netcdf/fix_sfc_climo/vegetation_type.viirs.igbp.0.03.nc
# use conus 0.01-degree viirs data (do not use for global fv3 grids).
##export VEG_FILE=/work/noaa/da/ggayno/save/ufs_utils.git/fv3.vegt.new.tundra.netcdf/fix_sfc_climo/vegetation_type.viirs.igbp.conus.0.01.nc

# Set working directory and directory where output files will reside.

export WORK_DIR=/work/noaa/stmp/$LOGNAME/work.sfc
export SAVE_DIR=/work/noaa/stmp/$LOGNAME/sfc.C${res}

# Should not have to touch anything below here.

if [[ $GRIDTYPE = "regional" ]]; then
  HALO=$(( $HALO + 1 ))
  export HALO
  ln -fs $FIX_FV3/C${res}_grid.tile7.halo${HALO}.nc $FIX_FV3/C${res}_grid.tile7.nc
  ln -fs $FIX_FV3/C${res}_oro_data.tile7.halo${HALO}.nc $FIX_FV3/C${res}_oro_data.tile7.nc
fi

export input_sfc_climo_dir=${BASE_DIR}/fix/fix_sfc_climo
export APRUN_SFC="srun"

ulimit -a
ulimit -s 199000000

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
