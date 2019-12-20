#!/bin/bash

set -x

source /apps/lmod/lmod/init/sh
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.7.0
module load hpss
module load prod_util
module list

PROJECT_CODE=fv3-cpu
QUEUE=batch

source config

if [ $EXTRACT_DATA == yes ]; then
  TEST1=$(sbatch --parsable --partition=service --ntasks=1 -t 4:00:00 -A $PROJECT_CODE -q $QUEUE -J get_data \
      -o log.data -e log.data ./get_data.sh)
  DEPEND="-d afterok:$TEST1"
else
  DEPEND=' '
fi

if [ $RUN_CHGRES == yes ]; then
  sbatch --parsable --ntasks-per-node=6 --nodes=3 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER} \
      -o log.${MEMBER} -e log.${MEMBER} ${DEPEND} run_chgres.sh ${MEMBER}

  MEMBER=1
  while [ $MEMBER -le 80 ]; do
    if [ $MEMBER -lt 10 ]; then
      MEMBER_CH="00${MEMBER}"
    else
      MEMBER_CH="0${MEMBER}"
    fi
    sbatch --parsable --ntasks-per-node=12 --nodes=1 -t 0:15:00 -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER_CH} \
      -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} run_chgres.sh ${MEMBER_CH}
    MEMBER=$(( $MEMBER + 1 ))
  done
fi
