#!/bin/sh

set -x

module purge
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load impi/18.0.1
module load lsf/10.1
module use /usrx/local/dev/modulefiles
module load NetCDF/4.5.0
module list

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/chgres_cube
export HOMEufs=$PWD/..
export OUTDIR=/gpfs/dell1/stmp/$LOGNAME/chgres_reg_tests

QUEUE="debug"
PROJECT_CODE="FV3GFS-T2O"

LOG_FILE=regression.log
SUM_FILE=summary.log
rm -f $LOG_FILE $SUM_FILE

export NCCMP=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp

export OMP_STACKSIZE=1024M
export APRUN=mpirun

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.restart
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.restart -W 0:15 -x -n 6 \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.restart.sh"

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.history
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c192.fv3.history -W 0:15 -x -n 6 -w 'ended(c96.fv3.restart)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c192.fv3.history.sh"

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.nemsio -W 0:15 -x -n 6 -w 'ended(c192.fv3.history)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.nemsio.sh"

export OMP_NUM_THREADS=4
export INPUT_DATA=${HOMEreg}/input_data/gfs.sigio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.sigio -W  0:15 -x -n 6 -w 'ended(c96.fv3.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.sigio.sh"

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/gfs.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.nemsio -W  0:15 -x -n 6 -w 'ended(c96.gfs.sigio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.nemsio.sh"

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -W  0:15 -x -n 6 -w 'ended(c96.gfs.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.regional.sh"

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c96.regional)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
