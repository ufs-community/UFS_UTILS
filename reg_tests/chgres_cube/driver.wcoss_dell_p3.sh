#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on WCOSS-Dell.
#
# Set OUTDIR to your working directory.  Set the PROJECT_CODE and QUEUE
# as appropriate.
#
# Invoke the script as follows with no arguments.  A series of daily-
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

export OUTDIR="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

QUEUE="${QUEUE:-debug}"
PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/chgres_cube

LOG_FILE=consistency.log

SUM_FILE=summary.log

rm -f $LOG_FILE $SUM_FILE

export OMP_STACKSIZE=1024M

export APRUN=mpirun

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.restart -W 0:15 -x -n 6 \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.restart.sh"

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c192.fv3.history -W 0:15 -x -n 6 -w 'ended(c96.fv3.restart)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c192.fv3.history.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.nemsio -W 0:15 -x -n 6 -w 'ended(c192.fv3.history)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=4
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.sigio -W  0:15 -x -n 6 -w 'ended(c96.fv3.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.sigio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.nemsio -W  0:15 -x -n 6 -w 'ended(c96.gfs.sigio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -W  0:15 -x -n 6 -w 'ended(c96.gfs.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.regional.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.netcdf -W 0:15 -x -n 12 -w 'ended(c96.regional)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.netcdf.sh"

#-----------------------------------------------------------------------------
# Initialize global C192 using GFS GRIB2 file.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c192.gfs.grib2 -W 0:05 -x -n 6 -w 'ended(c96.fv3.netcdf)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c192.gfs.grib2.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c192.gfs.grib2)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
