#!/bin/bash

set -x

source /apps/lmod/lmod/init/sh
module purge
module use -a /scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.7.0
module load hpss
module load prod_util
module load nco/4.7.0
module list

PROJECT_CODE=fv3-cpu
QUEUE=batch

source config

if [ $EXTRACT_DATA == yes ]; then

  rm -fr $EXTRACT_DIR
  mkdir -p $EXTRACT_DIR

  MEM=6000M
  WALLT="2:00:00"

  DATAH=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_hires \
      -o log.data.hires -e log.data.hires ./get_data.sh hires)
  DATA1=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp1 \
      -o log.data.grp1 -e log.data.grp1 ./get_data.sh grp1)
  DATA2=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp2 \
      -o log.data.grp2 -e log.data.grp2 ./get_data.sh grp2)
  DATA3=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp3 \
      -o log.data.grp3 -e log.data.grp3 ./get_data.sh grp3)
  DATA4=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp4 \
      -o log.data.grp4 -e log.data.grp4 ./get_data.sh grp4)
  DATA5=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp5 \
      -o log.data.grp5 -e log.data.grp5 ./get_data.sh grp5)
  DATA6=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp6 \
      -o log.data.grp6 -e log.data.grp6 ./get_data.sh grp6)
  DATA7=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp7 \
      -o log.data.grp7 -e log.data.grp7 ./get_data.sh grp7)
  DATA8=$(sbatch --parsable --partition=service --ntasks=1 --mem=$MEM -t $WALLT -A $PROJECT_CODE -q $QUEUE -J get_grp8 \
      -o log.data.grp8 -e log.data.grp8 ./get_data.sh grp8)

  DEPEND="-d afterok:$DATAH:$DATA1:$DATA2:$DATA3:$DATA4:$DATA5:$DATA6:$DATA7:$DATA8"

else

  DEPEND=' '

fi

if [ $RUN_CHGRES == yes ]; then
  MEMBER=hires
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
