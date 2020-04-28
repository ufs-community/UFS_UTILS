#!/bin/bash

. /apps/lmod/lmod/init/sh
module purge
module load intel/18.0.5.274
module load impi/2018.4.274
module load szip
module load hdf5
module load netcdf/4.2.1.1
module list

set -x

LOG_FILE=regression.log
SUM_FILE=summary.log
QUEUE="windfall"
PROJECT_CODE="emcda"

export home_dir=$PWD/../..
export WORK_DIR=/mnt/lfs3/projects/emcda/$LOGNAME/stmp/reg_tests.grid
rm -f $WORK_DIR

export APRUN=time
export APRUN_SFC=srun
export OMP_STACKSIZE=2048m
export machine=JET

export NCCMP=/apps/nccmp/1.8.2.1/intel/18.0.3.222/bin/nccmp
export HOMEreg=/mnt/lfs3/projects/emcda/George.Gayno/reg_tests/grid_gen/baseline_data

ulimit -a
ulimit -s unlimited

export OMP_NUM_THREADS=24
TEST1=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.uniform \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE ./c96.uniform.sh)

TEST2=$(sbatch --parsable --ntasks-per-node=24 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      --partition=xjet -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J grid_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST2 << EOF
#!/bin/sh
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF
