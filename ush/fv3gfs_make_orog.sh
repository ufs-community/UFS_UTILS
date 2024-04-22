#!/bin/bash

set -eux

nargv=$#

if [ $nargv -eq 6 ]; then  # cubed-sphere grid
  res=$1 
  tile=$2
  griddir=$3
  outdir=$4
  script_dir=$5
  orogfile="none"
  indir=$6
  workdir=$TEMP_DIR/C${res}/orog/tile$tile
else
  echo "Number of arguments must be 6 for cubic sphere grid"
  echo "Usage for cubic sphere grid: $0 resolution tile griddir outdir script_dir indir"
  exit 1
fi

executable=$exec_dir/orog
if [ ! -s $executable ]; then
  echo "executable does not exist"
  exit 1 
fi

if [ ! -s $workdir ]; then mkdir -p $workdir ;fi
if [ ! -s $outdir ]; then mkdir -p $outdir ;fi

OUTGRID="C${res}_grid.tile${tile}.nc"

# Make Orograraphy
echo "OUTGRID = $OUTGRID"
echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"

cd $workdir

cp ${indir}/topography.antarctica.ramp.30s.nc .
cp ${indir}/landcover.umd.30s.nc .
cp ${indir}/topography.gmted2010.30s.nc .
cp ${griddir}/$OUTGRID .
cp $executable .

echo $OUTGRID >> INPS
echo $orogfile >> INPS
if [ -z ${ocn+x} ]; then
  echo ".false." >> INPS
else
  echo ".true." >> INPS
fi 
echo "none" >> INPS

cat INPS
time $executable < INPS

if [ $? -ne 0 ]; then
  echo "ERROR in running $executable "
  exit 1
else
  outfile=oro.C${res}.tile${tile}.nc
  mv ./out.oro.nc $outdir/$outfile
  echo "file $outdir/$outfile is created"
  echo "Successfully running $executable "
  exit 0
fi

exit
