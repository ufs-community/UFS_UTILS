#!/bin/bash

#---------------------------------------------------------------------------
# Set directory names and file names for orog data.
#---------------------------------------------------------------------------

if [ ${CTAR} == 'C48' ] ; then
  OCNRES='500'
elif [ ${CTAR} == 'C96' ]; then
  OCNRES='500'
elif [ ${CTAR} == 'C192' ]; then
  OCNRES='050'
elif [ ${CTAR} == 'C384' ]; then
  OCNRES='025'
elif [ ${CTAR} == 'C768' ]; then
  OCNRES='025'
elif [ ${CTAR} == 'C1152' ]; then
  OCNRES='025'
else
  OCNRES='025'
fi

ORO_DIR="${CTAR}"
ORO_NAME="${CTAR}.mx${OCNRES}_oro_data"

if [ "$machine" = 'hera' ] ; then
  FIX_ORO_INPUT=/scratch1/NCEPDEV/global/glopara/fix/orog/20230615
elif [ "$machine" = 'wcoss2' ] ; then
  FIX_ORO_INPUT=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/orog/20230615
elif [ "$machine" = 'jet' ] ; then
  FIX_ORO_INPUT=/lfs4/HFIP/hfv3gfs/glopara/git/fv3gfs/fix/orog/20230615
elif [ "$machine" = 's4' ] ; then
  FIX_ORO_INPUT=/data/prod/glopara/fix/orog/20230615
else
  set +x
  echo ERROR machine $machine not supported.
  exit 3
fi
