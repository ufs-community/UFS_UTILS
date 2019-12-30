#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube regression tests on WCOSS-Dell.
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

module purge
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load impi/18.0.1
module load lsf/10.1
module use /usrx/local/dev/modulefiles
module load NetCDF/4.5.0
module list

export OUTDIR=/gpfs/dell1/stmp/$LOGNAME/chgres_reg_tests
QUEUE="debug"
PROJECT_CODE="GFS-DEV"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/chgres_cube

LOG_FILE=regression.log

SUM_FILE=summary.log

rm -f $LOG_FILE $SUM_FILE

export NCCMP=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp

export OMP_STACKSIZE=1024M

export APRUN=mpirun

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.restart
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.restart -W 0:15 -x -n 6 \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.restart.sh"

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.history
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c192.fv3.history -W 0:15 -x -n 6 -w 'ended(c96.fv3.restart)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c192.fv3.history.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.fv3.nemsio -W 0:15 -x -n 6 -w 'ended(c192.fv3.history)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.fv3.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=4
export INPUT_DATA=${HOMEreg}/input_data/gfs.sigio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.sigio -W  0:15 -x -n 6 -w 'ended(c96.fv3.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.sigio.sh"

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/gfs.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.gfs.nemsio -W  0:15 -x -n 6 -w 'ended(c96.gfs.sigio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.gfs.nemsio.sh"

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export INPUT_DATA=${HOMEreg}/input_data/fv3.nemsio
bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.regional -W  0:15 -x -n 6 -w 'ended(c96.gfs.nemsio)' \
        -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" "$PWD/c96.regional.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 \
     -w 'ended(c96.regional)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
