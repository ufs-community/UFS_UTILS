#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency test on S4.
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

#SBATCH -J snow
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o consistency.log
#SBATCH -e consistency.log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

export DATA="${WORK_DIR:-/scratch/short/users/$LOGNAME}"
export DATA="${DATA}/reg-tests/snow2mdl"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

rm -fr $DATA

export HOMEreg=/data/users/dhuber/save/nems/role.ufsutils/ufs_utils/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/data/prod/hpc-stack/intel-2022.1/grib_util/1.2.2/bin/wgrib
export WGRIB2=/data/prod/hpc-stack/intel-2022.1/impi-2022.1/wgrib2/2.0.8/bin/wgrib2

./snow2mdl.sh

exit 0
