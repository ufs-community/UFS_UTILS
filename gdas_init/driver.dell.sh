#!/bin/bash

set -x

module purge
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load impi/18.0.1
module load lsf/10.1
module use /usrx/local/dev/modulefiles
module load NetCDF/4.5.0
module list

PROJECT_CODE=GFS-DEV
QUEUE=dev

source config

if [ $EXTRACT_DATA == yes ]; then

  rm -fr $EXTRACT_DIR
  mkdir -p $EXTRACT_DIR

  MEM=6000M
  WALLT="2:00"

  case $gfs_ver in
    v14)
      bsub -e log.data.hires -o log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.hres.data -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v14.data.sh hires"
      bsub -e log.data.enkf -o log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.enkf.data -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v14.data.sh enkf"
      DEPEND="ended(get.hres.data,get.enkf.data)"
      ;;
    v15)
      bsub -e log.data.hires -o log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.hres.data -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh hires"
      bsub -e log.data.grp1 -o log.data.grp1 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp1 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp1"
      bsub -e log.data.grp2 -o log.data.grp2 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp2 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp2"
      bsub -e log.data.grp3 -o log.data.grp3 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp3 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp3"
      bsub -e log.data.grp4 -o log.data.grp4 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp4 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp4"
      bsub -e log.data.grp5 -o log.data.grp5 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp5 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp5"
      bsub -e log.data.grp6 -o log.data.grp6 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp6 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp6"
      bsub -e log.data.grp7 -o log.data.grp7 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp7 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp7"
      bsub -e log.data.grp8 -o log.data.grp8 -q $QUEUE -P $PROJECT_CODE -J get.enkf.grp8 -W $WALLT -n 1 \
        -R span[ptile=1] -M $MEM "./get_v15.data.sh grp8"
      DEPEND="ended(get.hres.data,get.enkf.grp1,get.enkf.grp2,get.enkf.grp3,get.enkf.grp4,get.enkf.grp5,get.enkf.grp6,get.enkf.grp7,get.enkf.grp8)"
      ;;
 esac

else

  DEPEND=' '

fi

if [ $RUN_CHGRES == yes ]; then
  MEMBER=hires
  NODES=3
  WALLT="0:15"
  if [ $CRES_HIRES == 'C768' ] ; then
    NODES=5
  elif [ $CRES_HIRES == 'C1152' ] ; then
    NODES=8
    WALLT="0:20"
  fi
  case $gfs_ver in
    v14)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -W $WALLT \
        -x -n 6 -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" -w $DEPEND \
        "./run_v14.chgres.sh ${MEMBER}"
      ;;
    v15)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -W $WALLT \
        -x -n 6 -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" -w $DEPEND \
        "./run_v15.chgres.sh ${MEMBER}"
      ;;
  esac
  MEMBER=1
  while [ $MEMBER -le 80 ]; do
    if [ $MEMBER -lt 10 ]; then
      MEMBER_CH="00${MEMBER}"
    else
      MEMBER_CH="0${MEMBER}"
    fi
    case $gfs_ver in
      v14)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W 0:15 \
          -x -n 12 -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" -w $DEPEND \
          "./run_v14.chgres.sh ${MEMBER_CH}"
        ;;
      v15)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W 0:15 \
          -x -n 12 -R "span[ptile=6]" -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" -w $DEPEND \
          "./run_v15.chgres.sh ${MEMBER_CH}"
      ;;
    esac
    MEMBER=$(( $MEMBER + 1 ))
  done
fi
