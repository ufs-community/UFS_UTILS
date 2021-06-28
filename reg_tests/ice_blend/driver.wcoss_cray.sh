#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on WCOSS-Cray.
#
# Set $DATA to your working directory.  Set the project code (BSUB -P)
# and queue (BSUB -q) as appropriate.
#
# Invoke the script as follows:  cat $script | bsub
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log.
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#BSUB -W 0:02
#BSUB -o consistency.log
#BSUB -e consistency.log
#BSUB -J iceb_regt
#BSUB -q debug
#BSUB -R "rusage[mem=2000]"
#BSUB -P GFS-DEV

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export DATA="${WORK_DIR:-/gpfs/hps3/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export WGRIB=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.5/exec/wgrib
export WGRIB2=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.5/exec/wgrib2
export COPYGB2=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.5/exec/copygb2
export COPYGB=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.5/exec/copygb
export CNVGRIB=/gpfs/hps/nco/ops/nwprod/grib_util.v1.0.5/exec/cnvgrib

export HOMEreg=/gpfs/hps3/emc/global/noscrub/George.Gayno/ufs_utils.git/reg_tests/ice_blend
export HOMEgfs=$PWD/../..

rm -fr $DATA

./ice_blend.sh

exit 0
