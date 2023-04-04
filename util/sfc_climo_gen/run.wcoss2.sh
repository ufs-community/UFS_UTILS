#!/bin/bash

#------------------------------------------------------------
# Run the sfc_climo_gen program stand-alone on WCOSS2 using
# pre-exiting 'grid' and 'orography files. See the
# sfc_gen.sh script for details.
#
# Set the configuration variables in sfc_gen.sh. Then,
# run this script as follows: 'qsub $script'
#------------------------------------------------------------

#PBS -o log
#PBS -e log
#PBS -q debug
#PBS -A GFS-DEV
#PBS -N grid_fv3
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=24:mem=75GB

set -x

# Adjust according to the PBS -l statement.
export APRUN_SFC="mpiexec -n 24 -ppn 24 -cpu-bind core"

export BASE_DIR=$PBS_O_WORKDIR/../..

$PBS_O_WORKDIR/sfc_gen.sh

exit
