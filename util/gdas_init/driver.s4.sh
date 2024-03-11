#!/bin/bash

#---------------------------------------------------------------------
# Driver script for running on S4.
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
module load prod_util/2.1.1

PROJECT_CODE=star
QUEUE=s4

source config

export machine=s4

if [ $EXTRACT_DATA == yes ]; then

  echo "HPSS is not accessible from S4!  If you wish to run just the chgres portion, set EXTRACT_DATA=NO in the config file and try again."
  exit

else  # do not extract data.

  DEPEND=' '

fi  # extract data?

if [ $RUN_CHGRES == yes ]; then

  export APRUN=srun
  NODES=3
  WALLT="0:15:00"
  export OMP_NUM_THREADS=1
  if [ $CRES_HIRES == 'C768' ] ; then
    NODES=5
  elif [ $CRES_HIRES == 'C1152' ] ; then
    NODES=8
    WALLT="0:20:00"
  fi
  case $gfs_ver in
    v12 | v13)
      export OMP_NUM_THREADS=4
      export OMP_STACKSIZE=1024M
      sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} --cpus-per-task=$OMP_NUM_THREADS \
        -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
        -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_pre-v14.chgres.sh ${CDUMP}
      ;;
    v14)
      sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
      -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v14.chgres.sh ${CDUMP}
      ;;
    v15)
      if [ "$CDUMP" = "gdas" ]; then
        sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
        -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v15.chgres.sh ${CDUMP}
      else
        sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
        -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v15.chgres.gfs.sh
      fi
      ;;
    v16retro)
      if [ "$CDUMP" = "gdas" ] ; then
        sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
        -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v16retro.chgres.sh hires
      else
        sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
        -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v16.chgres.sh ${CDUMP}
      fi
      ;;
    v16)
      sbatch --parsable --ntasks-per-node=6 --nodes=${NODES} -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${CDUMP} \
      -o log.${CDUMP} -e log.${CDUMP} ${DEPEND} run_v16.chgres.sh ${CDUMP}
      ;;
  esac

  if [ "$CDUMP" = "gdas" ]; then

    WALLT="0:15:00"

    if [ "$gfs_ver" = "v16retro" ]; then

      sbatch --parsable --ntasks-per-node=12 --nodes=1 -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_enkf \
      -o log.enkf -e log.enkf ${DEPEND} run_v16retro.chgres.sh enkf

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
              export OMP_NUM_THREADS=2
              export OMP_STACKSIZE=1024M
              sbatch --parsable --ntasks-per-node=12 --nodes=1 --cpus-per-task=$OMP_NUM_THREADS \
               -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER_CH} \
               -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} run_pre-v14.chgres.sh ${MEMBER_CH}
            ;;
          v14)
              sbatch --parsable --ntasks-per-node=12 --nodes=1 -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER_CH} \
              -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} run_v14.chgres.sh ${MEMBER_CH}
            ;;
          v15)
              sbatch --parsable --ntasks-per-node=12 --nodes=1 -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER_CH} \
              -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} run_v15.chgres.sh ${MEMBER_CH}
            ;;
          v16)
              sbatch --parsable --ntasks-per-node=12 --nodes=1 -t $WALLT -A $PROJECT_CODE -q $QUEUE -J chgres_${MEMBER_CH} \
              -o log.${MEMBER_CH} -e log.${MEMBER_CH} ${DEPEND} run_v16.chgres.sh ${MEMBER_CH}
            ;;
        esac
        MEMBER=$(( $MEMBER + 1 ))
      done

    fi # v16 retro?

  fi  # which CDUMP?

fi  # run chgres?
