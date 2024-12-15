#!/bin/bash

#---------------------------------------------------------------------------
# Set directory names and file names for the target grid orog data.
# A default ocean resolution (OCNRES) based on CTAR is used.
#---------------------------------------------------------------------------

if [ ${CTAR} == 'C48' ] ; then
  OCNRES='500'
elif [ ${CTAR} == 'C96' ]; then
  OCNRES='500'
elif [ ${CTAR} == 'C192' ]; then
  OCNRES='025'
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

#---------------------------------------------------------------------------
# When using the v15/v16 tiled warm restart data as input to the chgres
# step, the input grid orography is needed (there is no orography record
# in the restart files). Since the restart data was created before the
# latest orog version (20231027), need to use a previous version.
#---------------------------------------------------------------------------

if [ "$machine" = 'hera' ] ; then
  FIX_ORO_INPUT=/scratch1/NCEPDEV/global/glopara/fix/orog/20230615
elif [ "$machine" = 'wcoss2' ] ; then
  FIX_ORO_INPUT=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/orog/20230615
elif [ "$machine" = 'jet' ] ; then
  FIX_ORO_INPUT=/lfs5/HFIP/hfv3gfs/glopara/FIX/fix/orog/20230615
elif [ "$machine" = 's4' ] ; then
  FIX_ORO_INPUT=/data/prod/glopara/fix/orog/20230615
else
  set +x
  echo ERROR machine $machine not supported.
  exit 3
fi
