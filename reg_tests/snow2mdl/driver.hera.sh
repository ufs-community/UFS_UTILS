#!/bin/bash

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

export DATA="/scratch2/NCEPDEV/stmp1/George.Gayno/snow2mdl"


rm -fr $DATA

export HOMEreg=/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/scratch2/NCEPDEV/nwprod/NCEPLIBS/utils/grib_util.v1.1.1/exec/wgrib
export WGRIB2=/scratch2/NCEPDEV/nwprod/NCEPLIBS/utils/grib_util.v1.1.1/exec/wgrib2

./snow2mdl.sh

exit 0
