#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=sjet
#SBATCH --time 0:01
#SBATCH --account=emcda
#SBATCH --job-name=snow2mdl
#SBATCH -o regression.log
#SBATCH -e regression.log

set -x

module unload intel
module load intel/18.0.5.274
module list

export DATA="/mnt/lfs3/projects/emcda/George.Gayno/stmp/snow2mdl"


rm -fr $DATA

export HOMEreg=/lfs3/HFIP/emcda/George.Gayno/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/apps/wgrib/1.8.1.0b/bin/wgrib
export WGRIB2=/apps/wgrib2/0.1.9.6a/bin/wgrib2

./snow2mdl.sh

exit 0
