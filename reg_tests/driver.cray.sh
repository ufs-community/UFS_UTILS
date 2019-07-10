#!/bin/sh

set -x

module purge
module load PrgEnv-intel/5.2.56
module rm intel
module load intel/16.3.210
module load cray-mpich/7.2.0
module load craype-haswell
module load cray-netcdf
module load xt-lsfhpc/9.1.3
module list

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/chgres_cube
export HOMEufs=$PWD/..
export OUTDIR=/gpfs/hps3/stmp/George.Gayno/chgres_reg_tests

QUEUE="debug"
PROJECT_CODE="FV3GFS-T2O"

LOG_FILE=regression.log
SUM_FILE=summary.log

export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp

export OMP_STACKSIZE=1024M
export KMP_AFFINITY=disabled

export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/fv3.restart
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.restart -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.fv3.restart.sh"

export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/fv3.history
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c192.fv3.history -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        -w 'ended(c96.fv3.restart)' "export NODES=1; $PWD/c192.fv3.history.sh"

export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.nemsio -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        -w 'ended(c192.fv3.history)' "export NODES=1; $PWD/c96.fv3.nemsio.sh"

export OMP_NUM_THREADS=4
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/gfs.sigio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.sigio -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        -w 'ended(c96.fv3.nemsio)' "export NODES=1; $PWD/c96.gfs.sigio.sh"

export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/gfs.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.nemsio -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        -w 'ended(c96.gfs.sigio)' "export NODES=1; $PWD/c96.gfs.nemsio.sh"

export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        -w 'ended(c96.gfs.nemsio)' "export NODES=1; $PWD/c96.regional.sh"

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 -w 'ended(c96.regional)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
