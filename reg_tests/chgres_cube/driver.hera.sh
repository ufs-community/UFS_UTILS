#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on Hera.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube. 
# Set the PROJECT_CODE and QUEUE as appropriate.  To see which projects you 
# are authorized to use, type "account_params".
#
# Invoke the script with no arguments.  A series of daily-chained
# consistency tests will be submitted.  To check the queue, type:
# "squeue -u USERNAME".
#
# The run output will be stored in OUTDIR.  Log output from the suite
# will be in LOG_FILE.  Once the suite has completed, a summary is
# placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

export OUTDIR="${WORK_DIR:-/scratch2/NCEPDEV/stmp1/$LOGNAME}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $LOG_FILE* $SUM_FILE

export OMP_STACKSIZE=1024M

export APRUN=srun
export NCCMP=${NCCMP:-nccmp}
rm -fr $OUTDIR

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log01
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log02
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST2=$(sbatch --parsable --ntasks-per-node=6 --nodes=2 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      -o $LOG_FILE -e $LOG_FILE ./c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log03
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST3=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.nemsio \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log04
export OMP_NUM_THREADS=6   # should match cpus-per-task
TEST4=$(sbatch --parsable --ntasks-per-node=3 --cpus-per-task=6 --nodes=2 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.sigio \
      -o $LOG_FILE -e $LOG_FILE ./c96.gfs.sigio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log05
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST5=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.nemsio \
      -o $LOG_FILE -e $LOG_FILE ./c96.gfs.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log06
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST6=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      -o $LOG_FILE -e $LOG_FILE ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log07
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST7=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.netcdf.sh)

#-----------------------------------------------------------------------------
# Initialize global C192 using GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log08
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST8=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c192.gfs.grib2 \
      -o $LOG_FILE -e $LOG_FILE ./c192.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log09
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST9=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 25km.conus.gfs.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./25km.conus.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GFS PHYSICS.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log10
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST10=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.gfssdf.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./3km.conus.hrrr.gfssdf.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GSD PHYSICS AND SFC VARS FROM FILE.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log11
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST11=$(sbatch --parsable --ntasks-per-node=6 --nodes=2 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.newsfc.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./3km.conus.hrrr.newsfc.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING NAM GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log12
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST12=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.nam.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./13km.conus.nam.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log13
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST13=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.rap.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./13km.conus.rap.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM NA USING NCEI GFS GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log14
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST14=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.na.gfs.ncei.grib2.conus \
      -o $LOG_FILE -e $LOG_FILE ./13km.na.gfs.ncei.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log15
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST15=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf2wam \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.netcdf2wam.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------
LOG_FILE=consistency.log
sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q $QUEUE -d\
      afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6:$TEST7:$TEST8:$TEST9:$TEST10:$TEST11:$TEST12:$TEST13:$TEST14:$TEST15 << EOF
#!/bin/bash
grep -a '<<<' $LOG_FILE*  > $SUM_FILE
EOF

exit 0
