#!/bin/bash

#------------------------------------------------------------
# Run the weight_gen program to create  ESMF 'scrip'
# files for a gaussian grid.
#
# Run this script using the driver script for your machine.
#------------------------------------------------------------

set -x

CRES=${CRES:-"C48"}

WORK_DIR=${WORK_DIR:-"/lfs/h2/emc/stmp/$USER/weight_gen"}

UFS_DIR=${UFS_DIR:-"$PWD/../.."}

rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

${UFS_DIR}/exec/weight_gen "$CRES"

exit
