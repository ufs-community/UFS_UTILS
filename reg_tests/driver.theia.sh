#!/bin/sh

set -x

source /apps/lmod/lmod/init/sh
module purge
module load intel/18.1.163
module load impi/5.1.1.109
module load netcdf/4.3.0
module list

export HOMEufs=$PWD/..
export OUTDIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/chgres_reg_tests
export HOMEreg=/scratch4/NCEPDEV/da/noscrub/George.Gayno/reg_tests/chgres_cube

PROJECT_CODE="fv3-cpu"
QUEUE="debug"
export NCCMP=/apps/nccmp/1.8.2-gcc/bin/nccmp

LOG_FILE=regression.log
SUM_FILE=summary.log

export OMP_STACKSIZE=1024M

export APRUN=srun

rm -fr $OUTDIR

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.restart
TEST1=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.restart.sh)

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.history
TEST2=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c192.fv3.history.sh)

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
TEST3=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST2 ./c96.fv3.nemsio.sh)

export OMP_NUM_THREADS=6
export INPUT_DATA=${HOMEreg}/input_data/gfs.sigio
TEST4=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.sigio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST3 ./c96.gfs.sigio.sh)

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/gfs.nemsio
TEST5=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST4 ./c96.gfs.nemsio.sh)

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
TEST6=$(sbatch --parsable --ntasks=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST5 ./c96.regional.sh)

sbatch --ntasks=1 --mem=100M -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST6 << EOF
#!/bin/sh
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF

exit 0
