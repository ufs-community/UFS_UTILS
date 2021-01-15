#!/bin/bash

#----------------------------------------------------------------------
# Driver script for running on Cray.
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

  MEM=6000
  WALLT="2:00"

  case $gfs_ver in
    v12 | v13)
      bsub -o log.data.$CDUMP -e log.data.$CDUMP -q $QUEUE -P $PROJECT_CODE -J get.data.$CDUMP -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_pre-v14.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
          -R "rusage[mem=$MEM]" "./get_pre-v14.data.sh enkf"
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v14)
      bsub -o log.data.$CDUMP -e log.data.$CDUMP -q $QUEUE -P $PROJECT_CODE -J get.data.$CDUMP -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v14.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
          -R "rusage[mem=$MEM]" "./get_v14.data.sh enkf"
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v15)
      bsub -o log.data.hires -e log.data.hires -q $QUEUE -P $PROJECT_CODE -J get.data.hires -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh hires"
      bsub -o log.data.grp1 -e log.data.grp1 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf1 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp1"
      bsub -o log.data.grp2 -e log.data.grp2 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf2 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp2"
      bsub -o log.data.grp3 -e log.data.grp3 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf3 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp3"
      bsub -o log.data.grp4 -e log.data.grp4 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf4 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp4"
      bsub -o log.data.grp5 -e log.data.grp5 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf5 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp5"
      bsub -o log.data.grp6 -e log.data.grp6 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf6 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp6"
      bsub -o log.data.grp7 -e log.data.grp7 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf7 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp7"
      bsub -o log.data.grp8 -e log.data.grp8 -q $QUEUE -P $PROJECT_CODE -J get.data.enkf8 -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh grp8"
      DEPEND="-w ended(get.data.*)"
      ;;
    v16)
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v16.data2.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        for group in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8
        do
          bsub -o log.data.${group} -e log.data.${group} -q $QUEUE -P $PROJECT_CODE -J get.data.${group} -W $WALLT \
            -R "rusage[mem=$MEM]" "./get_v16.data2.sh ${group}"
        done
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
 esac

else

  DEPEND=' '

fi

if [ $RUN_CHGRES == yes ]; then
  MEM=2000
  QUEUE=dev
  MEMBER=$CDUMP
  WALLT="0:15"
  NUM_NODES=2
  case $gfs_ver in
    v12 | v13)
      export OMP_NUM_THREADS=2
      export OMP_STACKSIZE=1024M
      ;;
    *)
      export OMP_NUM_THREADS=1
      ;;
  esac
  export APRUN="aprun -j 1 -n 24 -N 12 -d ${OMP_NUM_THREADS} -cc depth"
  if [ $CRES_HIRES == 'C768' ] ; then
    WALLT="0:20"
    NUM_NODES=3
    export APRUN="aprun -j 1 -n 36 -N 12 -d ${OMP_NUM_THREADS} -cc depth"
  elif [ $CRES_HIRES == 'C1152' ] ; then
    WALLT="0:20"
    NUM_NODES=4
    export APRUN="aprun -j 1 -n 48 -N 12 -d ${OMP_NUM_THREADS} -cc depth"
  fi
  case $gfs_ver in
    v12 | v13)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_pre-v14.chgres.sh ${MEMBER}"
      ;;
    v14)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v14.chgres.sh ${MEMBER}"
      ;;
    v15)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v15.chgres.sh ${MEMBER}"
      ;;
    v16)
      bsub -e log.${MEMBER} -o log.${MEMBER} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16.chgres2.sh ${MEMBER}"
      ;;
  esac

  if [ "$CDUMP" = 'gdas' ]; then

    WALLT="0:15"
    NUM_NODES=1
    export APRUN="aprun -j 1 -n 12 -N 12 -d ${OMP_NUM_THREADS} -cc depth"
  
    MEMBER=1
    while [ $MEMBER -le 80 ]; do
      if [ $MEMBER -lt 10 ]; then
        MEMBER_CH="00${MEMBER}"
      else
        MEMBER_CH="0${MEMBER}"
      fi
      case $gfs_ver in
      v12 | v13)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -M $MEM -W $WALLT \
          -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_pre-v14.chgres.sh ${MEMBER_CH}"
        ;;
      v14)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -M $MEM -W $WALLT \
          -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v14.chgres.sh ${MEMBER_CH}"
        ;;
      v15)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -M $MEM -W $WALLT \
          -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v15.chgres.sh ${MEMBER_CH}"
      ;;
      v16)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -M $MEM -W $WALLT \
          -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16.chgres2.sh ${MEMBER_CH}"
      ;;
      esac
      MEMBER=$(( $MEMBER + 1 ))
    done

  fi

fi
