#!/bin/bash

set -x

CRES=${CRES:-"C48"}

WORK_DIR=${WORK_DIR:-"/lfs/h2/emc/stmp/$USER/weight_gen"}

UFS_DIR=${UFS_DIR:-"$PWD/../.."}

rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

${UFS_DIR}/exec/scrip.x "$CRES"

exit
