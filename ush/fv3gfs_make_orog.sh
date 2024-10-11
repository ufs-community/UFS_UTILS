#!/bin/bash

#-------------------------------------------------------------------
# Program Name: fv3gfs_make_orog
#
# Run the orography ('orog') program to create mask, terrain and
# GWD fields on the model tile.
#
# Author: GFDL Programmer
#
# History Log:
#   01/2018: Initial version.
#   04/2024: Some clean up.
#
# Usage:
#  Arguments:
#    res     - "C" Resolution of model grid - 48, 96, 768, etc.
#    tile    - Tile number.
#    griddir - Location of model 'grid' file.
#    outdir  - Location of the model orography file output by
#              the 'orog' program.
#    indir   - Location of input land mask and terrain data.
#
#  Input Files:
#    $GRIDFILE                         - The model 'grid' file 
#                                        containing georeference info.
#    topography.antarctica.ramp.30s.nc - RAMP terrain data.
#    landcover.umd.30s.nc              - Global land mask data.
#    topography.gmted2010.30s.nc       - Global USGS GMTED 2010
#                                        terrain data.
#
#  Output Files:
#    out.oro.nc - The model orography file (single tile).
#
# Condition codes:
#    0 - Normal termination.
#    1 - Incorrect number of script arguments.
#    2 - Program executable does not exits.
#    3 - Error running program.
#-------------------------------------------------------------------

set -eux

nargv=$#

if [ $nargv -eq 5 ]; then
  res=$1 
  tile=$2
  griddir=$3
  outdir=$4
  indir=$5
else
  set +x
  echo "FATAL ERROR: Number of arguments must be 5."
  echo "Usage: $0 resolution tile griddir outdir indir."
  set -x
  exit 1
fi

executable=$exec_dir/orog
if [ ! -s $executable ]; then
  set +x
  echo "FATAL ERROR, ${executable} does not exist."
  set -x
  exit 2 
fi

workdir=$TEMP_DIR/C${res}/orog/tile$tile

if [ ! -s $workdir ]; then mkdir -p $workdir ;fi
if [ ! -s $outdir ]; then mkdir -p $outdir ;fi

GRIDFILE="C${res}_grid.tile${tile}.nc"

# Make Orograraphy
set +x
echo "GRIDFILE = $GRIDFILE"
echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"
set -x

cd $workdir

ln -fs ${indir}/topography.antarctica.ramp.30s.nc .
ln -fs ${indir}/landcover.umd.30s.nc .
ln -fs ${indir}/topography.gmted2010.30s.nc .
ln -fs ${griddir}/$GRIDFILE .
ln -fs $executable .

#-------------------------------------------------------------------
# Set up program namelist. The entries are:
#
#  1 - GRIDFILE - model 'grid' file.
#  2 - Logical to output land mask only. When creating a grid
#      for the coupled model ("ocn" resolution is specified) 
#      this is true. The mask is then tweaked during the
#      ocean merge step before the 'orog' program is run again
#      (in fv3gfs_ocean_merge.sh) to create the full 'orog'
#      file. When false, the 'orog' program outputs the
#      full orography file.
#  3 - The input file from the ocean merge step. Defaults
#      to 'none' for this script.
#-------------------------------------------------------------------

echo $GRIDFILE > INPS
if [ -z ${ocn+x} ]; then
  echo ".false." >> INPS
else
  echo ".true." >> INPS
fi 
echo "none" >> INPS

cat INPS
time $executable < INPS
rc=$?

if [ $rc -ne 0 ]; then
  set +x
  echo "FATAL ERROR running $executable."
  set -x
  exit 3
else
  outfile=oro.C${res}.tile${tile}.nc
  mv ./out.oro.nc $outdir/$outfile
  set +x
  echo "Successfully ran ${executable}."
  echo "File $outdir/$outfile is created."
  set -x
  exit 0
fi
