#!/bin/bash

#SBATCH -J weight_gen
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

#-------------------------------------------------------------------------------
#
# Run the weight_gen program on Hera.
#
# Set WORK_DIR to your working directory.
#
# Set CRES to your desired resolution. Valid choices are:
#  - C48
#  - C96
#  - C128
#  - C192
#  - C384
#  - C768
#  - C1152
#  - C3072
#
# To run this script, do: 'sbatch $script'
#
#-------------------------------------------------------------------------------

set -x

UFS_DIR=$PWD/../..
source $UFS_DIR/sorc/machine-setup.sh > /dev/null 2>&1
module use $UFS_DIR/modulefiles
module load build.$target.intel
module list

export CRES="C48"

export WORK_DIR=/scratch2/NCEPDEV/stmp1/$USER/weight_gen

${UFS_DIR}/util/weight_gen/weight_gen.sh

exit
