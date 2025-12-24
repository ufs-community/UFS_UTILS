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
# By default, the utility will generate mapfiles for multiple configurations of
# DATM and/or CSG. Setting a WAVRES will also create mapfiles for to/from WW3.
#
# For the unstructured WW3 meshes and the CSG, no nstod_bilnr mapping from WW3
# should be used at runtime, since no destination mask is available in the CSG
# atmmesh. This is an ESMF limitation because the mesh is created internally by ESMF
# from the supergrid and mosaic files, and no capability exists to provide a mask
# file for the reduced grid.
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

# Set WAVRES, optionally
# export WAVRES=270k

# Loop over DATM resolutions and ocean resolutions
for ATMRES in 1760x880 1536x768 3072x1536; do
    for OCNRES in 100 050 025 008; do
        export ATMRES
        export OCNRES
        "${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh
    done
done

# Set CSG resolutions, optionally
# Loop over FV3 cube-sphere resolutions and ocean resolutions
#for ATMRES in C96 C192 C384 C1152; do
#    for OCNRES in 100 050 025; do
#        export ATMRES
#        export OCNRES
#        "${UFS_DIR}"/util/cmeps_mapfiles/make_mapfiles.sh
#    done
#done

exit
