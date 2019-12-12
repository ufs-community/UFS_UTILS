#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube regression tests on JET.
#
# Set OUTDIR to your working directory.  Set the PROJECT_CODE and QUEUE
# as appropriate.  To see which projects you are authorized to use,
# type "account_params".
#
# Invoke the script with no arguments.  A series of daily-
# chained jobs will be submitted.  To check the queue, type:
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

source /apps/lmod/lmod/init/sh
module purge
module load intel/18.0.5.274
module load impi/2018.4.274
module load szip
module load hdf5
module load netcdf/4.2.1.1
module list

export OUTDIR=/mnt/lfs3/projects/emcda/$LOGNAME/stmp/chgres_reg_tests
PROJECT_CODE="hfv3gfs"
QUEUE="debug"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/mnt/lfs3/projects/emcda/George.Gayno/reg_tests/chgres_cube

export NCCMP=/apps/nccmp/1.8.2.1/intel/18.0.3.222/bin/nccmp

LOG_FILE=regression.log
SUM_FILE=summary.log
rm -f $LOG_FILE $SUM_FILE

export OMP_STACKSIZE=1024M

export APRUN=srun

rm -fr $OUTDIR

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

TEST1=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

TEST2=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

TEST3=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST2 ./c96.fv3.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

TEST4=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=3 --cpus-per-task=6 -t 0:15:00 \
      -A $PROJECT_CODE -q $QUEUE -J c96.gfs.sigio -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST3 ./c96.gfs.sigio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

TEST5=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST4 ./c96.gfs.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

TEST6=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST5 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST6 << EOF
#!/bin/sh
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF

exit 0
