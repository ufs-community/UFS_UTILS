#!/bin/bash

#---------------------------------------------------------------------------
# Set directory names and file names for orog data
# The old and new (support fractional grid) orog data have different file names
#---------------------------------------------------------------------------

if [ "${FRAC_ORO:-"no"}" = "yes" ]; then
  if  [ ${CTAR} == 'C48' ] ; then
    OCNRES='500'
  elif [ ${CTAR} == 'C96' ] ; then
    OCNRES='100'
  elif [ ${CTAR} == 'C192' ] ; then
    OCNRES='050'
  elif [ ${CTAR} == 'C384' ] || [ ${CTAR} == 'C768' ] || [ ${CTAR} == 'C1152' ]; then
    OCNRES='025'
  fi
  ORO_DIR="${CTAR}.mx${OCNRES}_frac"
  ORO_NAME="oro_${CTAR}.mx${OCNRES}"
else
  ORO_DIR="${CTAR}"
  ORO_NAME="${CTAR}_oro_data"
fi
