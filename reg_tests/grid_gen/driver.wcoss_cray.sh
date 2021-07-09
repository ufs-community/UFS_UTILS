#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run grid generation consistency tests on WCOSS-Cray.
#
# Set WORK_DIR to your working directory. Set the PROJECT_CODE and QUEUE
# as appropriate.
#
# Invoke the script with no arguments.  A series of daily-
# chained jobs will be submitted.  To check the queue, type: "bjobs".
#
# Log output from the suite will be in LOG_FILE.  Once the suite
# has completed, a summary is placed in SUM_FILE.
#
# A test fails when its output does not match the baseline files as
# determined by the "nccmp" utility.  The baseline files are stored in
# HOMEreg
#
#-----------------------------------------------------------------------------

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

set -x

QUEUE="${QUEUE:-debug}"
PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
export WORK_DIR="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"
export WORK_DIR="${WORK_DIR}/reg-tests/grid-gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below here.
#-----------------------------------------------------------------------------

export home_dir=$PWD/../..
LOG_FILE=consistency.log
SUM_FILE=summary.log
export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
export APRUN_SFC="aprun -j 1 -n 24 -N 24"
export OMP_STACKSIZE=2048m
export OMP_NUM_THREADS=6
export machine=WCOSS_C
export KMP_AFFINITY=disabled
export NCCMP=/gpfs/hps3/emc/global/noscrub/George.Gayno/util/netcdf/nccmp
export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/grid_gen/baseline_data
export NCDUMP=/gpfs/hps/usrx/local/prod/NetCDF/4.2/intel/sandybridge/bin/ncdump

rm -fr $WORK_DIR

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------------
# C96 uniform grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.uniform -W 0:15 -M 2400 \
        -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/c96.uniform.sh"

#-----------------------------------------------------------------------------
# C96 uniform grid using viirs vegetation type data.
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J c96.viirs.vegt -W 0:15 -M 2400 \
        -w 'ended(c96.uniform)' -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/c96.viirs.vegt.sh"

#-----------------------------------------------------------------------------
# gfdl regional grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J gfdl.regional -W 0:10 -M 2400 \
        -w 'ended(c96.viirs.vegt)' -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/gfdl.regional.sh"

#-----------------------------------------------------------------------------
# ESG regional grid
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J esg.regional -W 0:10 -M 2400 \
        -w 'ended(gfdl.regional)' -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/esg.regional.sh"

#-----------------------------------------------------------------------------
# Regional GSL gravity wave drag.
#-----------------------------------------------------------------------------

bsub -e $LOG_FILE -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J reg.gsl.gwd -W 0:08 -M 2400 \
        -w 'ended(esg.regional)' -extsched 'CRAYLINUX[]' "export NODES=1; $PWD/regional.gsl.gwd.sh"

#-----------------------------------------------------------------------------
# Create summary log.
#-----------------------------------------------------------------------------

bsub -o $LOG_FILE -q $QUEUE -P $PROJECT_CODE -J summary -R "rusage[mem=100]" -W 0:01 -w 'ended(reg.gsl.gwd)' "grep -a '<<<' $LOG_FILE >> $SUM_FILE"

exit
