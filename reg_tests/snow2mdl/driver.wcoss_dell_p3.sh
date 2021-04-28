#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency test on WCOSS-Dell.
#
# Set $DATA to your working directory.  Set the project code (BSUB -P)
# and queue (BSUB -q) as appropriate.
#
# Invoke the script as follows:  cat $script | bsub
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is 
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#BSUB -W 0:02
#BSUB -o consistency.log
#BSUB -e consistency.log
#BSUB -J s2m_regt
#BSUB -q debug
#BSUB -R "affinity[core(1)]"
#BSUB -P GFS-DEV

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

set -x

export DATA="${WORK_DIR:-/gpfs/dell1/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/snow2mdl"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib
export WGRIB2=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib2

rm -fr $DATA

./snow2mdl.sh

exit 0
