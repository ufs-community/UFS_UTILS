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

#-------------------------------------------------------------------------------
#
# Run the mapfile generation program on Ursa.
#
# Set the ATM resolution (either DATM mesh or ATM CSG), the
# target ocean resolution and a target wave resolution (optional)
#
# To run this script, do: 'sbatch $script'
#
#-------------------------------------------------------------------------------

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
export datm_ver=20220805

export OUTPUT_DIR=/scratch4/NCEPDEV/stmp/$USER/cmeps_mapfiles
mkdir -p "$OUTPUT_DIR"

# currently supported DATM RTS
# default cfsr
export ATMRES=1760x880
export OCNRES=100
"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

# default gefs
export ATMRES=1536x768
export OCNRES=100
"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

# 3072x1536_cfsr, gfs
export ATMRES=3072x1536
export OCNRES=100
"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

# mx025_cfsr
export ATMRES=1760x880
export OCNRES=025
"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

# mx025_gefs
export ATMRES=1536x768
export OCNRES=025
"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

#RTOFSv3.0
#export ATMRES=3072x1536
#export OCNRES=008
#"${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh

exit
