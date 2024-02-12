#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on WCOSS2.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube. 
#
# Set the PROJECT_CODE and QUEUE as appropriate.
#
# Invoke the script with no arguments.  To check the queue, type:
# "qstat -u USERNAME".
#
# The run output will be stored in OUTDIR.  Log output will be placed
# in LOG_FILE??.  Once the suite has completed, a summary is placed
# in SUM_FILE.
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

export OUTDIR="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
QUEUE="${QUEUE:-dev}"

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

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/chgres_cube

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $LOG_FILE* $SUM_FILE

export OMP_STACKSIZE=1024M

export NCCMP=/lfs/h2/emc/global/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp
#export NCCMP=${NCCMP:-nccmp}
rm -fr $OUTDIR

this_dir=$PWD

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log01
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST1=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.fv3.restart -l select=1:ncpus=6:ompthreads=1:mem=10GB $PWD/c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log02
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST2=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c192.fv3.history -l select=1:ncpus=6:ompthreads=1:mem=10GB $PWD/c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log03
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST3=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.fv3.nemsio -l select=1:ncpus=6:ompthreads=1:mem=45GB $PWD/c96.fv3.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log04
export OMP_PLACES=cores
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core --depth 4"
TEST4=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:10:00 \
        -N c96.gfs.sigio -l select=1:ncpus=24:ompthreads=4:mem=45GB $PWD/c96.gfs.sigio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log05
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST5=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.gfs.nemsio -l select=1:ncpus=6:ompthreads=1:mem=35GB $PWD/c96.gfs.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize global C96 using GEFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log06
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST6=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.gefs.grib2 -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/c96.gefs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log07
export APRUN="mpiexec -n 12 -ppn 12 --cpu-bind core"
TEST7=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.fv3.netcdf -l select=1:ncpus=12:ompthreads=1:mem=80GB $PWD/c96.fv3.netcdf.sh)

#-----------------------------------------------------------------------------
# Initialize global C192 using GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log08
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST8=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c192.gfs.grib2 -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/c192.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log09
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST9=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 25km.conus.gfs.grib2.conus -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/25km.conus.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GFS PHYSICS.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log10
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST10=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 3km.conus.hrrr.gfssdf.grib2.conus -l select=1:ncpus=6:ompthreads=1:mem=75GB $PWD/3km.conus.hrrr.gfssdf.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GSD PHYSICS AND SFC VARS FROM FILE.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log11
export APRUN="mpiexec -n 12 -ppn 12 --cpu-bind core"
TEST11=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 3km.conus.hrrr.newsfc.grib2.conus -l select=1:ncpus=12:ompthreads=1:mem=75GB $PWD/3km.conus.hrrr.newsfc.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING NAM GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log12
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST12=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 13km.conus.nam.grib2.conus -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/13km.conus.nam.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log13
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST13=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 13km.conus.rap.grib2.conus -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/13km.conus.rap.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM NA USING NCEI GFS GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log14
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST14=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 13km.na.gfs.ncei.grib2.conus -l select=1:ncpus=6:ompthreads=1:mem=25GB $PWD/13km.na.gfs.ncei.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log15
export APRUN="mpiexec -n 12 -ppn 12 --cpu-bind core"
TEST15=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c96.fv3.netcdf2wam -l select=1:ncpus=12:ompthreads=1:mem=75GB $PWD/c96.fv3.netcdf2wam.sh)

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING  GFS PGRIB2+BGRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log16
export APRUN="mpiexec -n 6 -ppn 6 --cpu-bind core"
TEST16=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N 25km.conus.gfs.pbgrib2.conus -l select=1:ncpus=6:ompthreads=1:mem=15GB $PWD/25km.conus.gfs.pbgrib2.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

LOG_FILE=consistency.log
qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N chgres_summary -l select=1:ncpus=1:mem=100MB \
        -W depend=afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5:$TEST6:$TEST7:$TEST8:$TEST9:$TEST10:$TEST11:$TEST12:$TEST13:$TEST14:$TEST15:$TEST16 << EOF
#!/bin/bash
cd ${this_dir}
grep -a '<<<' ${LOG_FILE}?? | grep -v echo > $SUM_FILE
EOF

exit 0
