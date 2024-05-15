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

# the test needs to be once setting the filetype=ocean and once for ice
# one source 1/4deg grid (1440,1080)
# two destination grids: 1deg (360,320) and 1/2deg (720,576)
# only one would need to be tested in an RT

export SRCDIMS="1440,1080"
export DSTDIMS="360,320"
#export FILETYPE="ocean"
export FILETYPE="ice"
export WTSDIR="/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_GRIDGEN/rt_1191124/"
export GRDDIR="/scratch1/NCEPDEV/stmp4/Denise.Worthen/CPLD_GRIDGEN/rt_1191124/"
export DO_DEBUG=".false."

export OUTDIR_PATH="./"

cd ${OUTDIR_PATH}

# Stage input
# Use NCO to merge the first two MOM6 restarts to a single restart
# These two files will need to be available in the RT.

# The G-W will make these two files are available (having been
# retrieved from "somewhere") as MOM.res.nc and MOM.res_1.nc.
# Also assume that the ice restart is retrieved as ice.nc.
# The files produced will be ocean.mx[resolution].nc and
# ice.mx[resolution].nc
# Assume that the g-w will do the filename globbing and retrieval
# and rename (timestamp) the output files at the end.

#comment out for testing
#ncks -O -v Temp,Salt,h,u replay-2021032206/MOM.res.nc ocean.nc
#ncks -v v,sfc -A replay-2021032206/MOM.res_1.nc ocean.nc
#cp -f replay-2021032206/iced.2021-03-22-10800.nc ice.nc

edit_namelist < ocniceprep.nml.IN > ocniceprep.nml

$APRUN ./oiprep
