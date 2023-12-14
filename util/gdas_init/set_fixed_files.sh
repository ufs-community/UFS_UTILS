#!/bin/bash

#---------------------------------------------------------------------------
# Set directory names and file names for orog data.
#---------------------------------------------------------------------------

if [ -z "${OCNRES}" ]; then  # for uncoupled runs.
  ORO_DIR="${CTAR}"
  ORO_NAME="${CTAR}_oro_data"
else                         # for coupled runs.
  ORO_DIR="${CTAR}"
  ORO_NAME="${CTAR}.mx${OCNRES}_oro_data"
fi
