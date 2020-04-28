#!/bin/bash

module purge
module load EnvVars/1.0.2
module load lsf/10.1
module load ips/18.0.1.163
module load impi/18.0.1
module load NetCDF/4.5.0
module load HDF5-serial/1.10.1
module list

set -x

LOG_FILE=regression.log
SUM_FILE=summary.log
QUEUE="debug"
PROJECT_CODE="GFS-DEV"

export home_dir=$PWD/../..
export WORK_DIR=/gpfs/dell1/stmp/$LOGNAME/reg_tests.grid
rm -f $WORK_DIR

export APRUN=time
export APRUN_SFC="mpirun -l"
export OMP_STACKSIZE=2048m
export machine=WCOSS_DELL_P3

export NCCMP=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp
export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/grid_gen

ulimit -a
ulimit -s unlimited

export OMP_NUM_THREADS=24
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.uniform -W 0:15 -x -n 24 \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/c96.uniform.sh"

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -W 0:10 -x -n 24 -w 'ended(c96.uniform)' \
        -R "span[ptile=24]" -R "affinity[core(1):distribute=balance]" "$PWD/c96.regional.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c96.regional)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"
