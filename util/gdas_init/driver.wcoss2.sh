#!/bin/bash

#---------------------------------------------------------------------
# Driver script for running on WCOSS2.
#
# Edit the 'config' file before running.
#---------------------------------------------------------------------

set -x

compiler=${compiler:-"intel"}
source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module list

# Needed for NDATE utility
module load prod_util/2.0.8

PROJECT_CODE=GFS-DEV

source config

this_dir=$PWD

if [ $EXTRACT_DATA == yes ]; then

  rm -fr $EXTRACT_DIR
  mkdir -p $EXTRACT_DIR

  QUEUE=dev_transfer
  MEM=2GB
  WALLT="02:00:00"

  case $gfs_ver in
    v12 | v13)
      DATAH=$(qsub -V -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
              -N get_${CDUMP} -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_pre-v14.data.sh ${CDUMP})
      DEPEND="-W depend=afterok:$DATAH"
      if [ "$CDUMP" = "gdas" ] ; then
        DATA1=$(qsub -V -o log.data.enkf -e log.data.enkf -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
              -N get_enkf -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_pre-v14.data.sh enkf)
        DEPEND="-W depend=afterok:$DATAH:$DATA1"
      fi
      ;;
    v14)
        DATAH=$(qsub -V -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_${CDUMP} -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v14.data.sh ${CDUMP})
        DEPEND="-W depend=afterok:$DATAH"
      if [ "$CDUMP" = "gdas" ] ; then
        DATA1=$(qsub -V -o log.data.enkf -e log.data.enkf -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_enkf -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v14.data.sh enkf)
        DEPEND="-W depend=afterok:$DATAH:$DATA1"
      fi
      ;;
    v15)
      if [ "$CDUMP" = "gfs" ] ; then
        DATAH=$(qsub -V -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_${CDUMP} -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh ${CDUMP})
        DEPEND="-W depend=afterok:$DATAH"
      else
        DATAH=$(qsub -V -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_${CDUMP} -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh ${CDUMP})
        DATA1=$(qsub -V -o log.data.grp1 -e log.data.grp1 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp1 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp1)
        DATA2=$(qsub -V -o log.data.grp2 -e log.data.grp2 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp2 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp2)
        DATA3=$(qsub -V -o log.data.grp3 -e log.data.grp3 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp3 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp3)
        DATA4=$(qsub -V -o log.data.grp4 -e log.data.grp4 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp4 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp4)
        DATA5=$(qsub -V -o log.data.grp5 -e log.data.grp5 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp5 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp5)
        DATA6=$(qsub -V -o log.data.grp6 -e log.data.grp6 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp6 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp6)
        DATA7=$(qsub -V -o log.data.grp7 -e log.data.grp7 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp7 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp7)
        DATA8=$(qsub -V -o log.data.grp8 -e log.data.grp8 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp8 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v15.data.sh grp8)
        DEPEND="-W depend=afterok:$DATAH:$DATA1:$DATA2:$DATA3:$DATA4:$DATA5:$DATA6:$DATA7:$DATA8"
      fi
      ;;
    v16retro)
      DATAH=$(qsub -V -o log.data.v16retro -e log.data.v16retro -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
              -N get_v16retro -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16retro.data.sh ${CDUMP})
      DEPEND="-W depend=afterok:$DATAH"
      ;;
    v16)
      DATAH=$(qsub -V -o log.data.${CDUMP} -e log.data.${CDUMP} -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
              -N get_${CDUMP} -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh ${CDUMP})
      DEPEND="-W depend=afterok:$DATAH"
      if [ "$CDUMP" = "gdas" ] ; then
        DATA1=$(qsub -V -o log.data.grp1 -e log.data.grp1 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp1 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp1)
        DATA2=$(qsub -V -o log.data.grp2 -e log.data.grp2 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp2 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp2)
        DATA3=$(qsub -V -o log.data.grp3 -e log.data.grp3 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp3 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp3)
        DATA4=$(qsub -V -o log.data.grp4 -e log.data.grp4 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp4 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp4)
        DATA5=$(qsub -V -o log.data.grp5 -e log.data.grp5 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp5 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp5)
        DATA6=$(qsub -V -o log.data.grp6 -e log.data.grp6 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp5 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp6)
        DATA7=$(qsub -V -o log.data.grp7 -e log.data.grp7 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp7 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp7)
        DATA8=$(qsub -V -o log.data.grp8 -e log.data.grp8 -q $QUEUE -A $PROJECT_CODE -l walltime=$WALLT \
                -N get_grp8 -l select=1:ncpus=1:mem=$MEM -- ${this_dir}/get_v16.data.sh grp8)
        DEPEND="-W depend=afterok:$DATAH:$DATA1:$DATA2:$DATA3:$DATA4:$DATA5:$DATA6:$DATA7:$DATA8"
      fi
      ;;
 esac

else  # do not extract data.

  DEPEND=' '

fi  # extract data?

if [ $RUN_CHGRES == yes ]; then

  QUEUE=dev
  NODES=1
  TASKS_PER_NODE=24
  WALLT="0:15:00"
  MEM=75GB
  if [ $CRES_HIRES == 'C768' ] ; then
    MEM=250GB
  elif [ $CRES_HIRES == 'C1152' ] ; then
    MEM=350GB
    NODES=1
    TASKS_PER_NODE=48
    WALLT="0:20:00"
  fi
  NCPUS=${TASKS_PER_NODE}
  (( TASKS = NODES * TASKS_PER_NODE ))
  export APRUN="mpiexec -n $TASKS -ppn $TASKS_PER_NODE --cpu-bind core"
  case $gfs_ver in
    v12 | v13)
      export OMP_NUM_THREADS=4
      export OMP_STACKSIZE=1024M
      export OMP_PLACES=cores
      export APRUN="$APRUN --depth ${OMP_NUM_THREADS}"
      (( NCPUS = NCPUS * OMP_NUM_THREADS ))
      qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=${OMP_NUM_THREADS}:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
      -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_pre-v14.chgres.sh ${CDUMP}
      ;;
    v14)
      qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
      -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_v14.chgres.sh ${CDUMP}
      ;;
    v15)
      if [ "$CDUMP" = "gdas" ]; then
        qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
        -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_v15.chgres.sh ${CDUMP}
      else
        qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
        -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} ${this_dir}/run_v15.chgres.gfs.sh
      fi
      ;;
    v16retro)
      if [ "$CDUMP" = "gdas" ] ; then
        qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
        -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_v16retro.chgres.sh hires
      else
        qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
        -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_v16.chgres.sh ${CDUMP}
      fi
      ;;
    v16)
      qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
      -N chgres_${CDUMP} -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} -- ${this_dir}/run_v16.chgres.sh ${CDUMP}
      ;;
  esac

  if [ "$CDUMP" = "gdas" ]; then

    WALLT="0:15:00"
    MEM=75GB
    NODES=1
    TASKS_PER_NODE=12
    NCPUS=${TASKS_PER_NODE}
    (( TASKS = NODES * TASKS_PER_NODE))
    export APRUN="mpiexec -n $TASKS -ppn $TASKS_PER_NODE --cpu-bind core"

    if [ "$gfs_ver" = "v16retro" ]; then

      qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
      -N chgres_enkf -o log.enkf -e log.enkf ${DEPEND} -- ${this_dir}/run_v16retro.chgres.sh enkf

    else

      case $gfs_ver in # use threads for v12/13 data.
        v12 | v13)
         export OMP_NUM_THREADS=2
         export OMP_STACKSIZE=1024M
         export OMP_PLACES=cores
         export APRUN="$APRUN --depth ${OMP_NUM_THREADS}"
         (( NCPUS = NCPUS * OMP_NUM_THREADS ))
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
          v12 | v13)
              qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=${OMP_NUM_THREADS}:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
              -N chgres_${MEMBER_CH} -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} -- ${this_dir}/run_pre-v14.chgres.sh ${MEMBER_CH}
            ;;
          v14)
              qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
             -N chgres_${MEMBER_CH} -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} -- ${this_dir}/run_v14.chgres.sh ${MEMBER_CH}
            ;;
          v15)
              qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
             -N chgres_${MEMBER_CH} -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} -- ${this_dir}/run_v15.chgres.sh ${MEMBER_CH}
            ;;
          v16)
              qsub -V -l select=${NODES}:ncpus=${NCPUS}:ompthreads=1:mem=${MEM} -l walltime=$WALLT -A $PROJECT_CODE -q $QUEUE \
             -N chgres_${MEMBER_CH} -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} -- ${this_dir}/run_v16.chgres.sh ${MEMBER_CH}
            ;;
        esac
        MEMBER=$(( $MEMBER + 1 ))
      done

    fi # v16 retro?

  fi  # which CDUMP?

fi  # run chgres?
