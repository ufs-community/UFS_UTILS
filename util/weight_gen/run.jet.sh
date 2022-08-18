#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=xjet
#SBATCH --time 0:02
#SBATCH --account=emcda
#SBATCH --job-name=weight_gen
#SBATCH -o log
#SBATCH -e log

#-------------------------------------------------------------------------------
#
# Run the weight_gen program on Jet.
#
# Set WORK_DIR to your working directory.
#
# Set CRES to your desired resolution. Valid choices are:
#  - C48  => 192x94 and 192x96 gaussian
#  - C96  => 384x192 and 384x194 gaussian
#  - C128 => 512x256 and 512x258 gaussian
#  - C192 => 768x384 and 768x386 gaussian
#  - C384 => 1536x768 and 1536x770 gaussian
#  - C768 => 3072x1536 and 3072x1538 gaussian
#  - C1152 => 4608x2304 and 4608x2406 gaussian
#  - C3072 => 12288x6144 and 12288x6146 gaussian
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

export WORK_DIR=/lfs4/HFIP/emcda/$USER/stmp/weight_gen

${UFS_DIR}/util/weight_gen/weight_gen.sh

exit
