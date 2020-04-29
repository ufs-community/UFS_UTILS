#!/bin/bash

#BSUB -W 00:05
#BSUB -n 6
#BSUB -R span[ptile=6]
#BSUB -x
#BSUB -o regression.log
#BSUB -e regression.log
#BSUB -R "affinity[core(1)]"
#BSUB -q debug
#BSUB -M 2400
#BSUB -J glc_regt
#BSUB -P GFS-DEV

set -x

module purge
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load lsf/10.1
module load impi/18.0.1
module use /usrx/local/nceplibs/dev/NCEPLIBS/modulefiles
module load netcdf_parallel/4.7.4

set -x


export DATA=/gpfs/dell1/stmp/$LOGNAME/reg_tests.cycle


export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/global_cycle
export OMP_NUM_THREADS_CY=2

export APRUNCY="mpirun -l"

export NWPROD=$PWD/../..

export COMOUT=$DATA

export NCCMP=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/util/nccmp/nccmp-nc4.7.4/src/nccmp

reg_dir=$PWD

./C768.fv3gfs.sh

cp $DATA/summary.log  $reg_dir

exit
