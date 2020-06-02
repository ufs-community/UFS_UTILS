#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl regression test on WCOSS-Dell.
#
# Set $DATA to your working directory.  Set the project code (BSUB -P)
# and queue (BSUB -q) as appropriate.
#
# Invoke the script as follows:  cat $script | bsub
#
# Log output is placed in regression.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is 
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#BSUB -W 0:02
#BSUB -o regression.log
#BSUB -e regression.log
#BSUB -J s2m_regt
#BSUB -q debug
#BSUB -R "affinity[core(1)]"
#BSUB -P GFS-DEV

source ../../sorc/machine-setup.sh > /dev/null 2>&1
source ../../modulefiles/build.$target

set -x

export DATA=/gpfs/dell1/stmp/$LOGNAME/reg_tests.snow2mdl

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
