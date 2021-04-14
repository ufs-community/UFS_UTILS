#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend regression test on Orion.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in regression.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J ice_blend
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o regression.log
#SBATCH -e regression.log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

ulimit -s unlimited

export DATA="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export WGRIB=/apps/contrib/NCEPLIBS/orion/utils/grib_util.v1.2.0/exec/wgrib
export WGRIB2=/apps/contrib/NCEPLIBS/orion/utils/grib_util.v1.2.0/exec/wgrib2
export COPYGB=/apps/contrib/NCEPLIBS/lib/NCEPLIBS-grib_util/v1.1.1/exec/copygb
export COPYGB2=/apps/contrib/NCEPLIBS/orion/utils/grib_util.v1.2.0/exec/copygb2
export CNVGRIB=/apps/contrib/NCEPLIBS/orion/utils/grib_util.v1.2.0/exec/cnvgrib

export HOMEreg=/work/noaa/da/ggayno/save/ufs_utils.git/reg_tests/ice_blend
export HOMEgfs=$PWD/../..

rm -fr $DATA

./ice_blend.sh

exit 0
