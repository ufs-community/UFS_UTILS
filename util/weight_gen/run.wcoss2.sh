#!/bin/bash

#PBS -l walltime=00:02:00
#PBS -o log
#PBS -e log
#PBS -N weight_gen
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l select=1:ncpus=1:mem=100MB

#-------------------------------------------------------------------------------
#
# Run the weight_gen program on WCOSS2.
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
# To run this script, do: 'qsub $script'
#
#-------------------------------------------------------------------------------

cd $PBS_O_WORKDIR

set -x

UFS_DIR=$PBS_O_WORKDIR/../..
source $UFS_DIR/sorc/machine-setup.sh > /dev/null 2>&1
module use $UFS_DIR/modulefiles
module load build.$target.intel
module list

export CRES="C48"

export WORK_DIR=/lfs/h2/emc/stmp/$USER/weight_gen

${UFS_DIR}/util/weight_gen/weight_gen.sh

exit
