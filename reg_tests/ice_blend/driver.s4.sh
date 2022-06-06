#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on S4.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J ice_blend
#SBATCH -A s4
#SBATCH --open-mode=truncate
#SBATCH -o consistency.log
#SBATCH -e consistency.log
#SBATCH --ntasks=1
#SBATCH -q s4
#SBATCH -t 00:03:00

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

export DATA="${WORK_DIR:-/scratch/short/users/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export WGRIB=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/wgrib
export WGRIB2=/data/prod/hpc-stack/intel-2022.1/wgrib2/2.0.8/bin/wgrib2
export COPYGB=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/copygb
export COPYGB2=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/copygb2
export CNVGRIB=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/cnvgrib

export HOMEreg=/data/users/dhuber/save/nems/role.ufsutils/ufs_utils/reg_tests/ice_blend
export HOMEgfs=$PWD/../..

rm -fr $DATA

./ice_blend.sh

exit 0
