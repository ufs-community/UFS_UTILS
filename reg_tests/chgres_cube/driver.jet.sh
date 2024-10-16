#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on JET.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube. 
# Set the PROJECT_CODE and QUEUE as appropriate.  To see which projects you 
# are authorized to use, type "account_params".
#
# Invoke the script with no arguments.  A series of daisy-
# chained jobs will be submitted.  To check the queue, type:
# "squeue -u USERNAME".
#
# The run output will be stored in OUTDIR.  Standard output from
# each test will be placed in its own log file. Once the suite
# has completed, a summary of results is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export OUTDIR="${WORK_DIR:-/lfs5/HFIP/emcda/$LOGNAME/stmp}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
QUEUE="${QUEUE:-batch}"

export machine="jet"
export HDF5_DISABLE_VERSION_CHECK=2

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEufs=$PWD/../..

export HOMEreg=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/chgres_cube

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
export OMP_NUM_THREADS=1
TEST1=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log02
export OMP_NUM_THREADS=1
TEST2=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log03
export OMP_NUM_THREADS=1
TEST3=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c96.fv3.netcdf.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using GFS GRIB2 data.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log04
export OMP_NUM_THREADS=1
TEST4=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c192.gfs.grib2 \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c192.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log05
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST5=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 25km.conus.gfs.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./25km.conus.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GFS PHYSICS.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log06
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST6=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=2 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.gfssdf.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./3km.conus.hrrr.gfssdf.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GSD PHYSICS AND SFC VARS FROM FILE.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log07
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST7=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=3 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.newsfc.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./3km.conus.hrrr.newsfc.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING NAM GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log08
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST8=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.nam.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./13km.conus.nam.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log09
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST9=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.rap.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./13km.conus.rap.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM NA USING NCEI GFS GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log10
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST10=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=2 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 13km.na.gfs.ncei.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./13km.na.gfs.ncei.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log11
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST11=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=2 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf2wam \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c96.fv3.netcdf2wam.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS PGRIB2+BGRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log12
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST12=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 25km.conus.gfs.pbgrib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./25km.conus.gfs.pbgrib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using GEFS GRIB2 data.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log13
export OMP_NUM_THREADS=1
TEST13=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c96.gefs.grib2 \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./c96.gefs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP-SMOKE GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log14
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST14=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.rap-smoke.grib2.conus \
      --exclusive -o $LOG_FILE -e $LOG_FILE ./13km.conus.rap-smoke.grib2.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log
sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d\
       afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6:$TEST7:$TEST8:$TEST9:$TEST10:$TEST11:$TEST12:$TEST13:$TEST14 << EOF
#!/bin/bash
grep -a '<<<' $LOG_FILE*  > $SUM_FILE
EOF

exit 0
