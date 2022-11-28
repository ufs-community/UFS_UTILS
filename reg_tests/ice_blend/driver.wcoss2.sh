#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on WCOSS2.
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

#PBS -l walltime=00:02:00
#PBS -o consistency.log
#PBS -e consistency.log
#PBS -N iceb_regt
#PBS -q dev
#PBS -A GFS-DEV
#PBS -l select=1:ncpus=1:mem=2500MB

cd $PBS_O_WORKDIR

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module load grib_util/1.2.3
module load wgrib2/2.0.8
module list

set -x

export DATA="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/ice_blend
export HOMEgfs=$PBS_O_WORKDIR/../..

rm -fr $DATA

./ice_blend.sh

exit 0
