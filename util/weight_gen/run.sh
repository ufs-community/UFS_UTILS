#!/bin/bash

set -x

CRES=C48
#CRES=C96
#CRES=C128
#CRES=C192
#CRES=C384
#CRES=C768
#CRES=C1152
#CRES=C3072

WORK_DIR=/work/noaa/stmp/$USER/weight_gen

UFS_DIR=$PWD/../../exec

rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

${UFS_DIR}/scrip.x "C96"

exit
