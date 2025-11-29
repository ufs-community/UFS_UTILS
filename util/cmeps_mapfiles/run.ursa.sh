#!/bin/bash

# Edit account (-A) setting as required !

#SBATCH -J datmmesh_gen
#SBATCH -A nems
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -q debug
#SBATCH -t 00:15:00

set -x

UFS_DIR=$PWD/../..
source "${UFS_DIR}"/sorc/machine-setup.sh > /dev/null 2>&1
module use "${UFS_DIR}"/modulefiles
module load "build.${target}.intelllvm"
module list

export FIX_DIR=/scratch3/NCEPDEV/global/role.glopara/fix
export orog_ver=20240917
export ice_ver=20240416
export wav_ver=20250508

export OUTPUT_DIR=/scratch4/NCEPDEV/stmp/$USER/cmeps_mapfiles
mkdir -p "$OUTPUT_DIR"

export ATMRES=1760x880
#export ATMRES=C96
#export OCNRES=100
export OCNRES=008

"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

exit
