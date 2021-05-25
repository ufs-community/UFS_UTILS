#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run the chgres_cube consistency tests on JET.
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

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export OUTDIR="${WORK_DIR:-/lfs4/HFIP/emcda/$LOGNAME/stmp}"
export OUTDIR="${OUTDIR}/reg-tests/chgres-cube"

PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
QUEUE="${QUEUE:-debug}"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.  HOMEufs is the root
# directory of your UFS_UTILS clone.  HOMEreg contains the input data
# and baseline data for each test.
#-----------------------------------------------------------------------------

export HOMEufs=$PWD/../..

export HOMEreg=/lfs4/HFIP/emcda/George.Gayno/reg_tests/chgres_cube

LOG_FILE=consistency.log
SUM_FILE=summary.log
rm -f $LOG_FILE $SUM_FILE

export OMP_STACKSIZE=1024M

export APRUN=srun

rm -fr $OUTDIR

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 warm restart files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST1=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.restart \
      -o $LOG_FILE -e $LOG_FILE ./c96.fv3.restart.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using FV3 tiled history files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST2=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c192.fv3.history \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST1 ./c192.fv3.history.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST3=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST2 ./c96.fv3.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS sigio/sfcio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=6   # should match cpus-per-task
TEST4=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=3 --cpus-per-task=6 -t 0:15:00 \
      -A $PROJECT_CODE -q $QUEUE -J c96.gfs.sigio -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST3 ./c96.gfs.sigio.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using spectral GFS gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST5=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.gfs.nemsio \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST4 ./c96.gfs.nemsio.sh)

#-----------------------------------------------------------------------------
# Initialize regional C96 using FV3 gaussian nemsio files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST6=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.regional \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST5 ./c96.regional.sh)

#-----------------------------------------------------------------------------
# Initialize C96 using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST7=$(sbatch --parsable --partition=xjet --nodes=2 --ntasks-per-node=6 -t 0:10:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST6 ./c96.fv3.netcdf.sh)

#-----------------------------------------------------------------------------
# Initialize C192 using GFS GRIB2 data.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1
TEST8=$(sbatch --parsable --partition=xjet --nodes=1 --ntasks-per-node=6 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c192.gfs.grib2 \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST7 ./c192.gfs.grib2.sh)

#-----------------------------------------------------------------------------
# Initialize C96 WAM IC using FV3 gaussian netcdf files.
#-----------------------------------------------------------------------------

export OMP_NUM_THREADS=1   # should match cpus-per-task
TEST9=$(sbatch --parsable --partition=xjet --ntasks-per-node=6 --nodes=2 -t 0:07:00 -A $PROJECT_CODE -q $QUEUE -J c96.fv3.netcdf2wam \
      -o $LOG_FILE -e $LOG_FILE -d afterok:$TEST8 ./c96.fv3.netcdf2wam.sh)

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

sbatch --partition=xjet --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d afterok:$TEST9 << EOF
#!/bin/bash
grep -a '<<<' $LOG_FILE  > $SUM_FILE
EOF

exit 0
