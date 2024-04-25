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
else
  set +x
  echo "FATAL ERROR: Number of arguments must be 6 for cubic sphere grid."
  echo "Usage for cubic sphere grid: $0 resolution tile griddir outdir script_dir indir."
  set -x
  exit 1
fi

executable=$exec_dir/orog
if [ ! -s $executable ]; then
  set +x
  echo "FATAL ERROR, ${executable} does not exist."
  set -x
  exit 1 
fi

workdir=$TEMP_DIR/C${res}/orog/tile$tile

if [ ! -s $workdir ]; then mkdir -p $workdir ;fi
if [ ! -s $outdir ]; then mkdir -p $outdir ;fi

OUTGRID="C${res}_grid.tile${tile}.nc"

# Make Orograraphy
set +x
echo "OUTGRID = $OUTGRID"
echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"
set -x

cd $workdir

ln -fs ${indir}/topography.antarctica.ramp.30s.nc .
ln -fs ${indir}/landcover.umd.30s.nc .
ln -fs ${indir}/topography.gmted2010.30s.nc .
ln -fs ${griddir}/$OUTGRID .
ln -fs $executable .

echo $OUTGRID > INPS
echo $orogfile >> INPS
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
  exit 1
else
  outfile=oro.C${res}.tile${tile}.nc
  mv ./out.oro.nc $outdir/$outfile
  set +x
  echo "Successfully ran ${executable}."
  echo "File $outdir/$outfile is created."
  set -x
  exit 0
fi

exit
