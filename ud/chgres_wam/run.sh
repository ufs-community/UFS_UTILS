#!/bin/bash

#
# to submit this job to batch queue do  `sbatch <run.sh`
#
#-----------------------------------------------------------
# Invoke as: sbatch $script
#-----------------------------------------------------------

#SBATCH --ntasks-per-node=6 --nodes=2
#SBATCH -t 0:05:00
#SBATCH -A fv3-cpu
#SBATCH -q debug
#SBATCH -J fv3
#SBATCH -o ./log
#SBATCH -e ./log

set -x

PACKDIR=/scratch1/NCEPDEV/global/Henry.Juang/testGC/
EXPN=v15
LEVS=150

source $PACKDIR/UFS_UTILS/sorc/machine-setup.sh > /dev/null 2>&1
source $PACKDIR/UFS_UTILS/modulefiles/build.$target

# Threads useful when ingesting spectral gfs sigio files.
# Otherwise set to 1.
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1024M

WORKDIR=/scratch2/NCEPDEV/stmp1/$LOGNAME/chgres.wam.$EXPN.$LEVS
rm -fr $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR

ln -fs ${SLURM_SUBMIT_DIR}/config.$EXPN.l$LEVS.nml ./fort.41

date

srun $PACKDIR/UFS_UTILS/exec/chgres_cube

date

for tile in tile1 tile2 tile3 tile4 tile5 tile6
do

  mv out.atm.$tile.nc gfs_data.$tile.nc
  mv out.sfc.$tile.nc sfc_data.$tile.nc

done

exit 0
