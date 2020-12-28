#!/bin/bash

! Run the sfc_climo_gen program stand-alone on Dell.

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

export BASE_DIR=$LS_SUBCWD/../../..

source ${BASE_DIR}/sorc/machine-setup.sh > /dev/null 2>&1
module use ${BASE_DIR}/modulefiles
module load build.$target.intel
module list

export res=768

# where the 'grid', 'mosaic' and 'oro' files reside.
export FIX_FV3=${BASE_DIR}/fix/fix_fv3_gmted2010/C${res}
##export FIX_FV3=/gpfs/dell1/stmp/George.Gayno/my_grids/C3379

# uncomment for regional grid
##HALO=3
##export GRIDTYPE=regional

if [[ $GRIDTYPE = "regional" ]]; then
  HALO=$(( $HALO + 1 ))
  export HALO
  ln -fs $FIX_FV3/C${res}_grid.tile7.halo${HALO}.nc $FIX_FV3/C${res}_grid.tile7.nc
  ln -fs $FIX_FV3/C${res}_oro_data.tile7.halo${HALO}.nc $FIX_FV3/C${res}_oro_data.tile7.nc
fi

export WORK_DIR=/gpfs/dell1/stmp/$LOGNAME/work.sfc
export SAVE_DIR=/gpfs/dell1/stmp/$LOGNAME/sfc.C${res}
export input_sfc_climo_dir=${BASE_DIR}/fix/fix_sfc_climo
export APRUN_SFC="mpirun -l"

ulimit -a
ulimit -s unlimited

rm -fr $WORK_DIR $SAVE_DIR

${BASE_DIR}/ush/sfc_climo_gen.sh

exit
