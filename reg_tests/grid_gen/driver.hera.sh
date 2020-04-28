#!/bin/bash

. /apps/lmod/lmod/init/sh
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load hdf5/1.10.5
module load netcdf/4.7.0
module list

set -x

LOG_FILE=regression.log
SUM_FILE=summary.log
QUEUE="batch"
PROJECT_CODE="fv3-cpu"

export home_dir=$PWD/../..
export WORK_DIR=/scratch2/NCEPDEV/stmp1/$LOGNAME/reg_tests.grid
rm -f $WORK_DIR

export APRUN=time
export APRUN_SFC=srun
export OMP_STACKSIZE=2048m
export machine=HERA

export NCCMP=/apps/nccmp/1.8.5/intel/18.0.3.051/bin/nccmp
export HOMEreg=/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/grid_gen/baseline_data

ulimit -a
ulimit -s unlimited

export OMP_NUM_THREADS=24
TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
      -o $LOG_FILE -e $LOG_FILE ./c96.uniform.sh)

TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST2 << EOF
#!/bin/sh
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF
