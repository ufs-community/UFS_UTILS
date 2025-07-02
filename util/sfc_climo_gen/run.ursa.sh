#!/bin/bash

#------------------------------------------------------------
# Run the sfc_climo_gen program stand-alone on Ursa using
# pre-exiting 'grid' and 'orography' files. See the
# sfc_gen.sh script for details.
#
# Set the configuration variables in sfc_gen.sh. Then
# run this script as follows: 'sbatch $script'
#------------------------------------------------------------

#SBATCH -J sfc_climo_gen
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o log
#SBATCH -e log
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --mem=300g
#SBATCH -q debug
#SBATCH -t 00:15:00

set -x

export APRUN_SFC="srun"

export BASE_DIR=$SLURM_SUBMIT_DIR/../..

$SLURM_SUBMIT_DIR/sfc_gen.sh

exit
