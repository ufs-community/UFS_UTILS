#!/bin/ksh
set -ax

nargv=$#

export inorogexist=0

if [ $nargv -eq 5 ];  then  # lat-lon grid
  export lonb=$1
  export latb=$2
  export outdir=$3
  export is_latlon=1
  export ntiles=1
  export orogfile="none"
  export hist_dir=$4
  export TMPDIR=$5
  export workdir=$TMPDIR/latlon/orog/latlon_${lonb}x${latb}
elif [ $nargv -eq 6 ]; then  # cubed-sphere grid
  export res=$1 
  export tile_id=$2
  export griddir=$3
  export outdir=$4
  export ntiles=6
  export is_latlon=0
  export orogfile="none"
  export hist_dir=$5
  export TMPDIR=$6
elif [ $nargv -eq 8 ]; then  # input your own orography files
  export res=$1 
  export lonb=$1
  export latb=$1
  export tile=$2
  export griddir=$3
  export outdir=$4
  export ntiles=6
  export is_latlon=0
  export inputorog=$5
  export orogfile=$inputorog:t
  export inorogexist=1
  export hist_dir=$7
  export TMPDIR=$8
else
  echo "number of arguments must be 5 or 6 for cubic sphere grid and 4 for lat-lon grid"
  echo "Usage for cubic sphere grid: $0 resolution tile grid_dir out_dir hist_dir TMPDIR"
  exit 1
fi

export indir=$hist_dir

export exe_add_lake=$exec_dir/lakefrac
if [ ! -s $exe_add_lake ]; then
  echo "executable to add lake does not exist"
  exit 1 
fi

export exe_inland=$exec_dir/inland
if [ ! -s $exe_inland ]; then
  echo "executable to create inland mask does not exist"
  exit 1 
fi

export workdir=$TMPDIR/C${res}/orog/tiles
if [ ! -s $workdir ]; then mkdir -p $workdir ;fi

# Make lake data - 
# Create and add inland, lake_status, and lake_depth to the orography files
echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"

if [ $is_latlon -eq 1 ]; then
  echo "lake_frac has not been implemented for lat/lon grid"
  exit 0
else
  if [ $tile_id -eq 7 ]; then
    echo "lake_frac has not been implemented for FV3 stand-alone Regional (SAR) model"
    exit 0
  fi

  cd $workdir

# link all required files to the current work directory
  tile_beg=1
  tile_end=6
  tile=$tile_beg
  while [ $tile -le $tile_end ]; do
    outfile=oro.C${res}.tile${tile}.nc
    ln -s $outdir/$outfile .
    outgrid="C${res}_grid.tile${tile}.nc"
    ln -s $griddir/$outgrid .
    tile=$(( $tile + 1 ))
  done

# create inland mask and save it to the orography files
  cutoff=0.99
  rd=7
  $exe_inland $res $cutoff $rd

# create lake data for FV3 grid and save it to the orography files
  tile_beg=1
  tile_end=6
  tile=$tile_beg
  while [ $tile -le $tile_end ]; do
    outfile=oro.C${res}.tile${tile}.nc
    $exe_add_lake ${tile} ${res} ${indir} ${lake_cutoff}
    echo "lake fraction is added to $outfile"
    tile=$(( $tile + 1 ))
  done
  exit 0
fi
exit
