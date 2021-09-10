#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on WCOSS-Cray.
#
# Set WORK_DIR to a general working location outside the UFS_UTILS directory.
# The exact working directory (OUTDIR) will be WORK_DIR/reg_tests/chgres-cube.
# Set the PROJECT_CODE and QUEUE as appropriate. 
#
# Invoke the script with no arguments.   A series of daily-
# chained jobs will be submitted.  To check the queue, type: "bjobs".
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

export OUTDIR="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

QUEUE="${QUEUE:-dev}"
PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/chgres_cube

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f ${LOG_FILE}* $SUM_FILE

export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp

export OMP_STACKSIZE=1024M

export KMP_AFFINITY=disabled

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS GRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE1=${LOG_FILE}01
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE1 -o $LOG_FILE1 -q $QUEUE -P $PROJECT_CODE -J chgres01 -M 1000 -W 0:05 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/25km.conus.gfs.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GFS PHYSICS.
#-----------------------------------------------------------------------------

LOG_FILE2=${LOG_FILE}02
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE2 -o $LOG_FILE2 -q $QUEUE -P $PROJECT_CODE -J chgres02 -M 1000 -W 0:07 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/3km.conus.hrrr.gfssdf.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 3-KM USING HRRR GRIB2 file WITH GSD PHYSICS AND SFC VARS FROM FILE.
#-----------------------------------------------------------------------------

LOG_FILE3=${LOG_FILE}03
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 12 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE3 -o $LOG_FILE3 -q $QUEUE -P $PROJECT_CODE -J chgres03 -M 1000 -W 0:07 -extsched 'CRAYLINUX[]' \
        "export NODES=2; $PWD/3km.conus.hrrr.newsfc.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING NAM GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE4=${LOG_FILE}04
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE4 -o $LOG_FILE4 -q $QUEUE -P $PROJECT_CODE -J chgres04 -M 1000 -W 0:07 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/13km.conus.nam.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM USING RAP GRIB2 file WITH GSD PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE5=${LOG_FILE}05
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE5 -o $LOG_FILE5 -q $QUEUE -P $PROJECT_CODE -J chgres05 -M 1000 -W 0:07 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/13km.conus.rap.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 13-KM NA USING NCEI GFS GRIB2 file WITH GFS PHYSICS .
#-----------------------------------------------------------------------------

LOG_FILE6=${LOG_FILE}06
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE6 -o $LOG_FILE6 -q $QUEUE -P $PROJECT_CODE -J chgres06 -M 1000 -W 0:07 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/13km.na.gfs.ncei.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

LOG_FILE7=${LOG_FILE}07
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE7 -o $LOG_FILE7 -q $QUEUE -P $PROJECT_CODE -J chgres07 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.fv3.restart.sh"

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

LOG_FILE8=${LOG_FILE}08
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE8 -o $LOG_FILE8 -q $QUEUE -P $PROJECT_CODE -J chgres08 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c192.fv3.history.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE9=${LOG_FILE}09
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE9 -o $LOG_FILE9 -q $QUEUE -P $PROJECT_CODE -J chgres09 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.fv3.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

LOG_FILE10=${LOG_FILE}10
export OMP_NUM_THREADS=4
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE10 -o $LOG_FILE10 -q $QUEUE -P $PROJECT_CODE -J chgres10 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.gfs.sigio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE11=${LOG_FILE}11
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE11 -o $LOG_FILE11 -q $QUEUE -P $PROJECT_CODE -J chgres11 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.gfs.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

LOG_FILE12=${LOG_FILE}12
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE12 -o $LOG_FILE12 -q $QUEUE -P $PROJECT_CODE -J chgres12 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c96.regional.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE13=${LOG_FILE}13
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 12 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE13 -o $LOG_FILE13 -q $QUEUE -P $PROJECT_CODE -J chgres13 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=2; $PWD/c96.fv3.netcdf.sh"

#-----------------------------------------------------------------------------
# Initialize global C192 using GFS GRIB2 data.
#-----------------------------------------------------------------------------

LOG_FILE14=${LOG_FILE}14
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE14 -o $LOG_FILE14 -q $QUEUE -P $PROJECT_CODE -J chgres14 -M 1000 -W 0:05 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/c192.gfs.grib2.sh"

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

LOG_FILE15=${LOG_FILE}15
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 12 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE15 -o $LOG_FILE15 -q $QUEUE -P $PROJECT_CODE -J chgres15 -M 1000 -W 0:15 -extsched 'CRAYLINUX[]' \
        "export NODES=2; $PWD/c96.fv3.netcdf2wam.sh"

#-----------------------------------------------------------------------------
# Initialize CONUS 25-KM USING GFS PGRIB2+BGRIB2 files.
#-----------------------------------------------------------------------------

LOG_FILE16=${LOG_FILE}16
export OMP_NUM_THREADS=1
export APRUN="aprun -j 1 -n 6 -N 6 -d ${OMP_NUM_THREADS} -cc depth"
bsub -e $LOG_FILE16 -o $LOG_FILE16 -q $QUEUE -P $PROJECT_CODE -J chgres16 -M 1000 -W 0:05 -extsched 'CRAYLINUX[]' \
        "export NODES=1; $PWD/25km.conus.gfs.pbgrib2.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 -w 'ended(chgres*)' "grep -a '<<<' "*.log*" >> $SUM_FILE"

exit
