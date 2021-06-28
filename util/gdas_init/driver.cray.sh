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
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v15.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        for group in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8
        do
          bsub -o log.data.${group} -e log.data.${group} -q $QUEUE -P $PROJECT_CODE -J get.data.${group} -W $WALLT \
            -R "rusage[mem=$MEM]" "./get_v15.data.sh ${group}"
        done
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v16retro)
      bsub -o log.data.$CDUMP -e log.data.$CDUMP -q $QUEUE -P $PROJECT_CODE -J get.data.$CDUMP -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v16retro.data.sh ${CDUMP}"
      DEPEND="-w ended(get.data.${CDUMP})"
      ;;
    v16)
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "rusage[mem=$MEM]" "./get_v16.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        for group in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8
        do
          bsub -o log.data.${group} -e log.data.${group} -q $QUEUE -P $PROJECT_CODE -J get.data.${group} -W $WALLT \
            -R "rusage[mem=$MEM]" "./get_v16.data.sh ${group}"
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
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_pre-v14.chgres.sh ${CDUMP}"
      ;;
    v14)
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v14.chgres.sh ${CDUMP}"
      ;;
    v15)
      if [ "$CDUMP" = "gdas" ]; then
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
           -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v15.chgres.sh ${CDUMP}"
      else
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
           -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v15.chgres.gfs.sh"
      fi
      ;;
    v16retro)
      if [ "$CDUMP" = "gdas" ] ; then
        bsub -e log.gdas -o log.gdas -q $QUEUE -P $PROJECT_CODE -J chgres_gdas -M $MEM -W $WALLT \
           -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16retro.chgres.sh hires"
      else
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
           -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16.chgres.sh ${CDUMP}"
      fi 
      ;;
    v16)
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16.chgres.sh ${CDUMP}"
      ;;
  esac

  if [ "$CDUMP" = 'gdas' ]; then

    WALLT="0:15"
    NUM_NODES=1
    export APRUN="aprun -j 1 -n 12 -N 12 -d ${OMP_NUM_THREADS} -cc depth"
  
    if [ "$gfs_ver" = "v16retro" ]; then

      bsub -e log.enkf -o log.enkf -q $QUEUE -P $PROJECT_CODE -J chgres_enkf -M $MEM -W $WALLT \
         -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16retro.chgres.sh enkf"

    else

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
          -extsched 'CRAYLINUX[]' $DEPEND "export NODES=$NUM_NODES; ./run_v16.chgres.sh ${MEMBER_CH}"
      ;;
      esac
      MEMBER=$(( $MEMBER + 1 ))
    done

    fi # is this v16 retro?

  fi

fi
