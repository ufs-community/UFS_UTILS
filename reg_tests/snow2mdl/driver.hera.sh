#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency test on Hera.
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

export DATA="${WORK_DIR:-/scratch2/NCEPDEV/stmp1/$LOGNAME}"
export DATA="${DATA}/reg-tests/snow2mdl"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

rm -fr $DATA

export HOMEreg=/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/scratch2/NCEPDEV/nwprod/NCEPLIBS/utils/grib_util.v1.1.1/exec/wgrib
export WGRIB2=/scratch2/NCEPDEV/nwprod/NCEPLIBS/utils/grib_util.v1.1.1/exec/wgrib2

./snow2mdl.sh

exit 0
