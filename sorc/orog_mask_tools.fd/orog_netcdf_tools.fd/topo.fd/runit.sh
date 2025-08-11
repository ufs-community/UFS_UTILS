#!/bin/sh

# Run on Ursa.

#SBATCH --ntasks=1
#SBATCH --mem=50g
#SBATCH -t 0:03:00
#SBATCH -A fv3-cpu
#SBATCH -q debug
#SBATCH -J fv3
#SBATCH -o ./log
#SBATCH -e ./log

set -x

source ../../../machine-setup.sh > /dev/null 2>&1
module use ../../../../modulefiles
module load build.$target.intelllvm
module list

../../../../exec/topo.exe
