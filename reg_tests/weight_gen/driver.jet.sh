#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run weight_gen consistency test on Jet.
#
# Set $DATA to your working directory.  Set the project code 
# (SBATCH --account) as appropriate.
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

#SBATCH --nodes=1
#SBATCH --partition=sjet
#SBATCH --time 0:01
#SBATCH --account=emcda
#SBATCH --job-name=weight_gen
#SBATCH -o consistency.log
#SBATCH -e consistency.log

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

export DATA="${WORK_DIR:-/lfs4/HFIP/emcda/$LOGNAME/stmp}"
export DATA="${DATA}/reg-tests/weight_gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/lfs4/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/weight_gen
export HOMEufs=$PWD/../..

./weight_gen.sh

exit 0
