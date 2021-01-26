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

#----------------------------------------------------------------------
# Extract data.
#----------------------------------------------------------------------

if [ "$EXTRACT_DATA" = "yes" ]; then

  rm -fr $EXTRACT_DIR
  mkdir -p $EXTRACT_DIR

  QUEUE=dev_transfer

  MEM=6000M
  WALLT="2:00"

  case $gfs_ver in
    v12 | v13 )
      bsub -o log.data.$CDUMP -e log.data.$CDUMP -q $QUEUE -P $PROJECT_CODE -J get.data.$CDUMP -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_pre-v14.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
          -R "affinity[core(1)]" -M $MEM "./get_pre-v14.data.sh enkf"
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v14)
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v14.data.sh ${CDUMP}"
      
      if [ "$CDUMP" = "gdas" ] ; then
        bsub -o log.data.enkf -e log.data.enkf -q $QUEUE -P $PROJECT_CODE -J get.data.enkf -W $WALLT \
          -R "affinity[core(1)]" -M $MEM "./get_v14.data.sh enkf"
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v15)
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        for group in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8
        do
          bsub -o log.data.enkf.${group} -e log.data.enkf.${group} -q $QUEUE -P $PROJECT_CODE -J get.data.enkf.${group} -W $WALLT \
            -R "affinity[core(1)]" -M $MEM "./get_v15.data.sh ${group}"
        done
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
    v16retro)
      bsub -o log.data.v16retro -e log.data.v16retro -q $QUEUE -P $PROJECT_CODE -J get.data.v16retro -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v16retro.data.sh ${CDUMP}"
      DEPEND="-w ended(get.data.v16retro)"
      ;;
    v16)
      bsub -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J get.data.${CDUMP} -W $WALLT \
        -R "affinity[core(1)]" -M $MEM "./get_v16.data.sh ${CDUMP}"
      if [ "$CDUMP" = "gdas" ] ; then
        for group in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8
        do
          bsub -o log.data.enkf.${group} -e log.data.enkf.${group} -q $QUEUE -P $PROJECT_CODE -J get.data.enkf.${group} -W $WALLT \
            -R "affinity[core(1)]" -M $MEM "./get_v16.data.sh ${group}"
        done
      fi
      DEPEND="-w ended(get.data.*)"
      ;;
 esac

else  # do not extract data.

  DEPEND=' '

fi # extract data?

#----------------------------------------------------------------------
# Run chgres.
#----------------------------------------------------------------------

if [ "$RUN_CHGRES" = "yes" ]; then

  QUEUE=dev2
  WALLT="0:15"
  export OMP_NUM_THREADS=1
  NODES="-n 18 -R "span[ptile=9]""
  export APRUN="mpirun"
  if [ "$CRES_HIRES" = "C768" ] ; then
    NODES="-n 24 -R "span[ptile=6]""
  elif [ "$CRES_HIRES" = "C1152" ] ; then
    NODES="-n 36 -R "span[ptile=6]""
    WALLT="0:20"
  fi

  case $gfs_ver in
    v12 | v13)
      export OMP_STACKSIZE=1024M
      export OMP_NUM_THREADS=2
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
        -x $NODES -R "affinity[core(${OMP_NUM_THREADS}):distribute=balance]" $DEPEND \
        "./run_pre-v14.chgres.sh ${CDUMP}"
      ;;
    v14)
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
        -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
        "./run_v14.chgres.sh ${CDUMP}"
      ;;
    v15)
      if [ "$CDUMP" = "gdas" ]; then
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v15.chgres.sh ${CDUMP}"
      else
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v15.chgres.gfs.sh"
      fi
      ;;
    v16retro)
      if [ "$CDUMP" = "gdas" ]; then
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v16retro.chgres.sh hires"
      else
        bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v16.chgres.sh ${CDUMP}"
      fi
     ;;
    v16)
      bsub -e log.${CDUMP} -o log.${CDUMP} -q $QUEUE -P $PROJECT_CODE -J chgres_${CDUMP} -W $WALLT \
        -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
        "./run_v16.chgres.sh ${CDUMP}"
     ;;
  esac

#----------------------------------------------------------------------
# If selected, run chgres for enkf members.
#----------------------------------------------------------------------

  if [ "$CDUMP" = "gdas" ]; then

    NODES="-n 18 -R "span[ptile=9]""
    WALLT="0:15"

    if [ "$gfs_ver" = "v16retro" ]; then

      bsub -e log.enkf -o log.enkf -q $QUEUE -P $PROJECT_CODE -J chgres_enkf -W $WALLT \
        -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
        "./run_v16retro.chgres.sh enkf"

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
      v16)
        bsub -e log.${MEMBER_CH} -o log.${MEMBER_CH} -q $QUEUE -P $PROJECT_CODE -J chgres_${MEMBER_CH} -W $WALLT \
          -x $NODES -R "affinity[core(1):distribute=balance]" $DEPEND \
          "./run_v16.chgres.sh ${MEMBER_CH}"
        ;;
      esac
      MEMBER=$(( $MEMBER + 1 ))
    done

    fi # is this v16 retro?

  fi # is this gdas? then process enkf.

fi  # run chgres?
