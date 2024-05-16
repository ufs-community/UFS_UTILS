#!/bin/bash
set -eux

function edit_namelist {

    sed -e "s/SRCDIMS/$SRCDIMS/g" \
        -e "s/DSTDIMS/$DSTDIMS/g" \
	-e "s/FILETYPE/$FILETYPE/g" \
        -e "s|WTSDIR|$WTSDIR|g" \
	-e "s|GRDDIR|$GRDDIR|g" \
        -e "s/DO_DEBUG/$DO_DEBUG/g"
}

APRUN=${APRUN:-"srun --nodes=1 -A nems "}

# Two possible input files: ocean.nc and ice.nc
# One possible source grid, 1/4 deg (mx025: 1440,1080)
# Two possible destination grids, 1/2deg and 1deg (mx050: 720,576 and mx100: 360,320)
# The files produced will be ocean.mx[dest res].nc and ice.mx[des res].nc
# The program needs to execute twice, once for ocean and once for ice

# For the purposes of the RT, we can have staged input files retrieved from
# https://noaa-ufs-gefsv13replay-pds.s3.amazonaws.com/2021/03/2021032206/
# The required single mom file was created using
# ncks -O -v Temp,Salt,h,u replay-2021032206/MOM.res.nc ocean.nc
# ncks -v v,sfc -A replay-2021032206/MOM.res_1.nc ocean.nc
# I am assuming that the g-w will do the filename globbing and retrieval,
# process the NCO command and rename (timestamp) the output files at the end.

export SRCDIMS="1440,1080"
export DSTDIMS="360,320"
#export FILETYPE="ocean"
export FILETYPE="ice"
export WTSDIR="/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_GRIDGEN/rt_1191124/"
export GRDDIR="/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_GRIDGEN/rt_1191124/"
export DO_DEBUG=".false."

export OUTDIR_PATH="./"

cd ${OUTDIR_PATH}

edit_namelist < ocniceprep.nml.IN > ocniceprep.nml

$APRUN ./oiprep
