#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency test on WCOSS2.
#
# Set $DATA to your working directory.  Set the project code (PBS -A)
# and queue (PBS -q) as appropriate.
#
# Invoke the script as follows:  qsub $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is 
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#PBS -l walltime=00:05:00
#PBS -o consistency.log
#PBS -e consistency.log
#PBS -N s2m_regt
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l select=1:ncpus=1:mem=2500MB

cd $PBS_O_WORKDIR

source ../../modulefiles/modulefile.global_emcsfc_snow2mdl.wcoss2_cray
module load grib_util/1.2.2
module load wgrib2/2.0.8
module list

set -x

export DATA="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/snow2mdl"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export HOMEreg=/lfs/h2/emc/eib/noscrub/George.Gayno/ufs_utils.git/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..

rm -fr $DATA

./snow2mdl.sh

exit 0
