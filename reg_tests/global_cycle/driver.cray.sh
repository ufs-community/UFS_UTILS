#!/bin/bash

#BSUB -oo regression.log
#BSUB -eo regression.log
#BSUB -q debug
#BSUB -P GDAS-T2O
#BSUB -J cycle_regt
#BSUB -M 2400
#BSUB -W 00:05
#BSUB -extsched 'CRAYLINUX[]'

module load PrgEnv-intel cfp-intel-sandybridge/1.1.0
module list

export NODES=1

export DATA=/gpfs/hps3/stmp/$LOGNAME/reg_tests.cycle

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=4

export KMP_AFFINITY=disabled

export APRUNCY="aprun -n 6 -N 6 -j 1 -d $OMP_NUM_THREADS_CY -cc depth"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
