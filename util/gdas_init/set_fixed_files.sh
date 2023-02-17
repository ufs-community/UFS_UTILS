#!/bin/bash

if [ "${FRAC_ORO:-"no"}" = "yes" ]; then
  if  [ ${CTAR} == 'C48' ] ; then
    OCNRES='400'
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
