#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run weight_gen consistency test on WCOSS2.
#
# Set $DATA to your working directory.  Set the project code (PBS -A)
# and queue (PBS -q) as appropriate.
#
# Invoke the script as follows:  qsub $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'nccmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#PBS -l walltime=00:02:00
#PBS -o consistency.log
#PBS -e consistency.log
#PBS -N wgt_regt
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l select=1:ncpus=1:mem=250MB

set -x

cd $PBS_O_WORKDIR

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

export DATA="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/weight_gen"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/weight_gen
export HOMEufs=$PBS_O_WORKDIR/../..

export NCCMP=/lfs/h2/emc/global/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp

./weight_gen.sh

exit 0
