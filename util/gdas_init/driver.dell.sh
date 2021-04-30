#!/bin/bash

#----------------------------------------------------------------------
# Driver script for running on Dell.
#
# Edit the 'config' file before running.
#----------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

PROJECT_CODE=GFS-DEV

source config

if [ $EXTRACT_DATA == yes ]; then

  rm -fr $EXTRACT_DIR
  mkdir -p $EXTRACT_DIR

  QUEUE=dev_transfer

  MEM=6000M
  WALLT="2:00"

  case $gfs_ver in
    v12 | v13 )
      bsub -o log.data.hires -e log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.data.hires -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_pre-v14.data.sh hires"
      bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_pre-v14.data.sh enkf"
      DEPEND="-w ended(get.data.*)"
      ;;
    v14)
      bsub -o log.data.hires -e log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.data.hires -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v14.data.sh hires"
      bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v14.data.sh enkf"
      DEPEND="-w ended(get.data.*)"
      ;;
    v15)
      bsub -o log.data.hires -e log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.data.hires -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh hires"
      bsub -o log.data.grp1 -e log.data.grp1 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf1 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp1"
      bsub -o log.data.grp2 -e log.data.grp2 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf2 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp2"
      bsub -o log.data.grp3 -e log.data.grp3 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf3 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp3"
      bsub -o log.data.grp4 -e log.data.grp4 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf4 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp4"
      bsub -o log.data.grp5 -e log.data.grp5 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf5 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp5"
      bsub -o log.data.grp6 -e log.data.grp6 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf6 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp6"
      bsub -o log.data.grp7 -e log.data.grp7 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf7 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp7"
      bsub -o log.data.grp8 -e log.data.grp8 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf8 -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh grp8"
      DEPEND="-w ended(get.data.*)"
      ;;
 esac

else

  DEPEND=' '

fi

if [ $RUN_CHGRES == yes ]; then
  QUEUE=dev
  MEMBER=hires
  WALLT="0:15"
  export OMP_NUM_THREADS=1
  NODES="-n 18 -R "span[ptile=9]""
  export APRUN="mpirun"
  if [ $CRES_HIRES == 'C768' ] ; then
    NODES="-n 24 -R "span[ptile=6]""
  elif [ $CRES_HIRES == 'C1152' ] ; then
    NODES="-n 36 -R "span[ptile=6]""
    WALLT="0:20"
  fi
  case $gfs_ver in
    v12 | v13)
      export OMP_STACKSIZE=1024M
      export OMP_NUM_THREADS=2
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -W $WALLT \
        -x $NODES -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" $DEPEND \
        "./run_pre-v14.chgres.sh ${MEMBER}"
      ;;
    v14)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -W $WALLT \
        -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
        "./run_v14.chgres.sh ${MEMBER}"
      ;;
    v15)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -W $WALLT \
        -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
        "./run_v15.chgres.sh ${MEMBER}"
      ;;
  esac

  NODES="-n 18 -R "span[ptile=9]""
  WALLT="0:15"
  MEMBER=1
  while [ $MEMBER -le 80 ]; do
    if [ $MEMBER -lt 10 ]; then
      MEMBER_CH="00${MEMBER}"
    else
      MEMBER_CH="0${MEMBER}"
    fi
    case $gfs_ver in
      v12 | v13)
        export OMP_STACKSIZE=1024M
        export OMP_NUM_THREADS=2
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W $WALLT \
          -x $NODES -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" $DEPEND \
          "./run_pre-v14.chgres.sh ${MEMBER_CH}"
        ;;
      v14)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v14.chgres.sh ${MEMBER_CH}"
        ;;
      v15)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v15.chgres.sh ${MEMBER_CH}"
      ;;
    esac
    MEMBER=$(( $MEMBER + 1 ))
  done
fi
