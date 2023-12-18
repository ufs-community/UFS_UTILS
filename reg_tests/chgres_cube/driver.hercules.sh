#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on Hercules.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube.
# Set the PROJECT_CODE and QUEUE as appropriate.  To see which projects you 
# are authorized to use, type "saccount_params".
#
# Invoke the script with no arguments.  A series of daily-chained
# consistency tests will be submitted.  To check the queue, type:
# "squeue -u $LOGNAME".
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

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

ulimit -s unlimited

export OUTDIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

PROJECT_CODE="${PROJECT_CODE:-nesdis-rdo2}"
QUEUE="${QUEUE:-batch}"

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

export HOMEreg=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/chgres_cube

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $SUM_FILE ${LOG_FILE}*

export OMP_STACKSIZE=1024M

export APRUN=srun

export machine=hercules

export NCCMP=${NCCMP:-nccmp}

rm -fr $OUTDIR

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

LOG_FILE1=${LOG_FILE}01
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      -o $LOG_FILE1 -e $LOG_FILE1 ./c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

LOG_FILE2=${LOG_FILE}02
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST2=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      --open-mode=append -o $LOG_FILE2 -e $LOG_FILE2 ./c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE3=${LOG_FILE}03
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST3=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.nemsio \
      --open-mode=append -o $LOG_FILE3 -e $LOG_FILE3 ./c96.fv3.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

LOG_FILE4=${LOG_FILE}04
export OMP_NUM_THREADS=6  # needs to match cpus-per-task
TEST4=$(sbatch --parsable --ntasks-per-node=3 --cpus-per-task=6 --nodes=2 --mem=75G -t 0:25:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.sigio \
      --open-mode=append -o $LOG_FILE4 -e $LOG_FILE4 ./c96.gfs.sigio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE5=${LOG_FILE}05
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST5=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.nemsio \
      --open-mode=append -o $LOG_FILE5 -e $LOG_FILE5 ./c96.gfs.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE6=${LOG_FILE}06
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST6=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      --open-mode=append -o $LOG_FILE6 -e $LOG_FILE6 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Initialize global C192 using GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE7=${LOG_FILE}07
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST7=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c192.gfs.grib2 \
      --open-mode=append -o $LOG_FILE7 -e $LOG_FILE7 ./c192.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize global C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE8=${LOG_FILE}08
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST8=$(sbatch --parsable --ntasks-per-node=6 --nodes=2 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf \
      --open-mode=append -o $LOG_FILE8 -e $LOG_FILE8 ./c96.fv3.netcdf.sh)

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE9=${LOG_FILE}09
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST9=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=100G -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf2wam \
      --open-mode=append -o $LOG_FILE9 -e $LOG_FILE9 ./c96.fv3.netcdf2wam.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE10=${LOG_FILE}10
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST10=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 25km.conus.gfs.grib2 \
      --open-mode=append -o $LOG_FILE10 -e $LOG_FILE10 ./25km.conus.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GFS PHYSICS.
#-----------------------------------------------------------------------------

LOG_FILE11=${LOG_FILE}11
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST11=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.gfssdf.grib2 \
      --open-mode=append -o $LOG_FILE11 -e $LOG_FILE11 ./3km.conus.hrrr.gfssdf.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GSD PHYSICS AND SFC VARS FROM FILE.
#-----------------------------------------------------------------------------

LOG_FILE12=${LOG_FILE}12
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST12=$(sbatch --parsable --ntasks-per-node=6 --nodes=2 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 3km.conus.hrrr.newsfc.grib2 \
      --open-mode=append -o $LOG_FILE12 -e $LOG_FILE12 ./3km.conus.hrrr.newsfc.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING NAM GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE13=${LOG_FILE}13
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST13=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.nam.grib2 \
      --open-mode=append -o $LOG_FILE13 -e $LOG_FILE13 ./13km.conus.nam.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE14=${LOG_FILE}14
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST14=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 13km.conus.rap.grib2 \
      --open-mode=append -o $LOG_FILE14 -e $LOG_FILE14 ./13km.conus.rap.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM NA USING NCEI GFS GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE15=${LOG_FILE}15
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST15=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 13km.na.gfs.ncei.grib2 \
      --open-mode=append -o $LOG_FILE15 -e $LOG_FILE15 ./13km.na.gfs.ncei.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS PGRIB2+BGRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE16=${LOG_FILE}16
export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST16=$(sbatch --parsable --ntasks-per-node=12 --nodes=1 --mem=75G -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J 25km.conus.gfs.pbgrib2 \
      --open-mode=append -o $LOG_FILE16 -e $LOG_FILE16 ./25km.conus.gfs.pbgrib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using GEFS GRIB2 file.
#-----------------------------------------------------------------------------

LOG_FILE17=${LOG_FILE}17
export OMP_NUM_THREADS=1  # needs to match cpus-per-task
TEST17=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 --mem=75G -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c96.gefs.grib2 \
      --open-mode=append -o $LOG_FILE17 -e $LOG_FILE17 ./c96.gefs.grib2.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE \
       -d afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6:$TEST7:$TEST8:$TEST9:$TEST10:$TEST11:$TEST12:$TEST13:$TEST14:$TEST15:$TEST16:$TEST17 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}*  > $SUM_FILE
EOF

exit 0
