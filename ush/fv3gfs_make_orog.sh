#!/bin/bash

set -eux

nargv=$#

inorogexist=0

if [ $nargv -eq 5 ];  then  # lat-lon grid
  lonb=$1
  latb=$2
  outdir=$3
  script_dir=$4
  is_latlon=1
  orogfile="none"
  hist_dir=$5
  workdir=$TEMP_DIR/latlon/orog/latlon_${lonb}x${latb}
elif [ $nargv -eq 6 ]; then  # cubed-sphere grid
  res=$1 
  lonb=$1
  latb=$1
  tile=$2
  griddir=$3
  outdir=$4
  script_dir=$5
  is_latlon=0
  orogfile="none"
  hist_dir=$6
  workdir=$TEMP_DIR/C${res}/orog/tile$tile
elif [ $nargv -eq 8 ]; then  # input your own orography files
  res=$1 
  lonb=$1
  latb=$1
  tile=$2
  griddir=$3
  outdir=$4
  is_latlon=0
  inputorog=$5
  script_dir=$6
  orogfile=$inputorog:t
  inorogexist=1
  hist_dir=$7
  workdir=$TEMP_DIR/C${res}/orog/tile$tile
else
  echo "Number of arguments must be 6 for cubic sphere grid"
  echo "Usage for cubic sphere grid: $0 resolution tile griddir outdir script_dir hist_dir"
  exit 1
fi

indir=$hist_dir
executable=$exec_dir/orog
if [ ! -s $executable ]; then
  echo "executable does not exist"
  exit 1 
fi

if [ ! -s $workdir ]; then mkdir -p $workdir ;fi
if [ ! -s $outdir ]; then mkdir -p $outdir ;fi

#jcap is for Gaussian grid
#jcap=`expr $latb - 2 `
jcap=0
NF1=0
NF2=0
mtnres=1
efac=0
blat=0
NR=0

if [ $is_latlon -eq 1 ]; then
  OUTGRID="none"
else
  OUTGRID="C${res}_grid.tile${tile}.nc"
fi

# Make Orograraphy
echo "OUTGRID = $OUTGRID"
echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"

cd $workdir

cp ${indir}/thirty.second.antarctic.new.bin fort.15
cp ${indir}/landcover30.fixed.nc .
cp ${indir}/gmted2010.30sec.nc  fort.235
if [ $inorogexist -eq 1 ]; then
   cp $inputorog .
fi   
     
if [ $is_latlon -eq 0 ]; then
   cp ${griddir}/$OUTGRID .
fi
cp $executable .

echo  $mtnres $lonb $latb $jcap $NR $NF1 $NF2 $efac $blat > INPS
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
  if [ $is_latlon -eq 1 ]; then
     outfile=oro.${lonb}x${latb}.nc
  else
     outfile=oro.C${res}.tile${tile}.nc
  fi

  mv ./out.oro.nc $outdir/$outfile
  echo "file $outdir/$outfile is created"
  echo "Successfully running $executable "
  exit 0
fi

exit
