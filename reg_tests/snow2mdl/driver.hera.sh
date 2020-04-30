#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl regression test on Hera.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in regression.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline files
# as determined by the 'cmp' command.  The baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J snow
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o regression.log
#SBATCH -e regression.log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

set -x

module load intel

export DATA="/scratch2/NCEPDEV/stmp1/$LOGNAME/reg_tests.snow2mdl"

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
