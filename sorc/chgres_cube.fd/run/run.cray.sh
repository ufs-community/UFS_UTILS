#!/bin/sh

#BSUB -oo log
#BSUB -eo log
#BSUB -q debug
#BSUB -J chgres_fv3
#BSUB -P FV3GFS-T2O
#BSUB -W 0:10
#BSUB -M 1000
#BSUB -extsched 'CRAYLINUX[]'

set -x

export NODES=1
# threads useful when using gfs sigio files as input
export OMP_NUM_THREADS=1
#export OMP_NUM_THREADS=4
export OMP_STACKSIZE=1024M

WORK_DIR=/gpfs/hps3/stmp/$LOGNAME/chgres_fv3
rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

cp $LS_SUBCWD/config.C48.cray.nml ./fort.41

EXEC_DIR=$LS_SUBCWD/../../../exec

export KMP_AFFINITY=disabled
aprun -j 1 -n 6 -N 6 -d${OMP_NUM_THREADS} -cc depth $EXEC_DIR/chgres_cube.exe
#aprun -j 1 -n 18 -N 18 -d${OMP_NUM_THREADS} -cc depth $EXEC_DIR/chgres_cube.exe

exit
