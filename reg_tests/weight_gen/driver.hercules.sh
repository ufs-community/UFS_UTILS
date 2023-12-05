#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run weight_gen consistency test on Hercules.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'nccmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J weight_gen
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

export DATA="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/weight_gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/work/noaa/nems/role-nems/ufs_utils/hercules/reg_tests/weight_gen
export HOMEufs=$PWD/../..

./weight_gen.sh

exit 0
