#!/bin/bash

# Edit account (-A) setting as required !

#SBATCH -J datmmesh_gen
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

#-------------------------------------------------------------------------------
#
# Run the datmmesh generation program on Ursa.
#
# Set NX and NY to your desired resolution. Valid choices are:
#   NX=3072, NY=1536 => 3072x1536 mesh
#   NX=1760, NY=880 => 1760x880 mesh
#   NX=1536, NY=768 => 1536x768 mesh
#
# Set N2S=.true. for N->S orientation, 1536x768 mesh only.
#
# To run this script, do: 'sbatch $script'
#
#-------------------------------------------------------------------------------

set -x

UFS_DIR=$PWD/../..
source $UFS_DIR/sorc/machine-setup.sh > /dev/null 2>&1
module use $UFS_DIR/modulefiles
module load build.$target.intelllvm
module list

export OUTPUT_DIR=/scratch4/NCEPDEV/stmp/$USER/datmmesh_gen
mkdir -p $OUTPUT_DIR

export NX=3072
export NY=1536
export N2S=.false.
${UFS_DIR}/util/datmmesh_gen/datmmesh.sh

export NX=1760
export NY=880
export N2S=.false.
${UFS_DIR}/util/datmmesh_gen/datmmesh.sh

export NX=1536
export NY=768
export N2S=.true.
${UFS_DIR}/util/datmmesh_gen/datmmesh.sh

exit
