#!/bin/bash

#BSUB -W 0:02
#BSUB -o regression.log
#BSUB -e regression.log
#BSUB -J s2m_regt
#BSUB -q debug
#BSUB -R "affinity[core(1)]"
#BSUB -P GFS-DEV

set -x


export DATA=/gpfs/dell1/stmp/George.Gayno/reg_tests.snow2mdl


rm -fr $DATA

export HOMEreg=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib
export WGRIB2=/gpfs/dell1/nco/ops/nwprod/grib_util.v1.0.6/exec/wgrib2

./snow2mdl.sh

exit 0
