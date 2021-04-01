#!/bin/bash

echo
echo "CREATE AND ADD INLAND, LAKE_STATUS, AND LAKE_DEPTH TO THE OROGRAPHY FILES"
echo

set -eux

outdir=$orog_dir
indir=$topo

if [ $gtype != uniform ] && [ $gtype != regional_gfdl ]; then
  echo "lakefrac has only been implemented for 'uniform' and 'regional_gfdl'."
  exit 0
fi

exe_add_lake=$exec_dir/lakefrac
if [ ! -s $exe_add_lake ]; then
  echo "executable to add lake does not exist"
  exit 1 
fi

exe_inland=$exec_dir/inland
if [ ! -s $exe_inland ]; then
  echo "executable to create inland mask does not exist"
  exit 1 
fi

workdir=$TEMP_DIR/C${res}/orog/tiles
if [ ! -s $workdir ]; then mkdir -p $workdir ;fi

# Make lake data - 
# Create and add inland, lake_status, and lake_depth to the orography files

echo "workdir = $workdir"
echo "outdir = $outdir"
echo "indir = $indir"

cd $workdir

# link all required files to the current work directory


if [ $gtype == uniform ]; then
  tile_beg=1
  tile_end=6
  tile=$tile_beg
  while [ $tile -le $tile_end ]; do
    outfile=oro.C${res}.tile${tile}.nc
    ln -s $outdir/$outfile .
    outgrid="C${res}_grid.tile${tile}.nc"
    ln -s $grid_dir/$outgrid .
    tile=$(( $tile + 1 ))
  done
fi

if [ $gtype == regional_gfdl ]; then
  tile_beg=7
  tile_end=7
  tile=7
  outfile=oro.C${res}.tile${tile}.nc
  ln -s $outdir/$outfile .
  outgrid="C${res}_grid.tile${tile}.nc"
  ln -s $grid_dir/$outgrid .
fi

# create inland mask and save it to the orography files

cutoff=0.99
rd=7
if [ $gtype == uniform ]; then
  $APRUN $exe_inland $res $cutoff $rd g
fi
if [ $gtype == regional_gfdl ]; then
  $APRUN $exe_inland $res $cutoff $rd r
fi
err=$?
if [ $err != 0 ]; then
  set +x
  echo ERROR CREATING INLAND MASK
  exit $err
fi

# create lake data for FV3 grid and save it to the orography files

if [ $machine = WCOSS_C ]; then
  touch ./lake.txt
  tile=$tile_beg
  while [ $tile -le $tile_end ]; do
    echo "$exe_add_lake ${tile} ${res} ${indir} ${lake_cutoff}" >> ./lake.txt
    tile=$(( $tile + 1 ))
  done
  aprun -j 1 -n 6 -N 6 -d 1 -cc depth cfp ./lake.txt
  err=$?
  if [ $err != 0 ]; then
    set +x
    echo ERROR CREATING LAKE FRACTION
    exit $err
  fi
  rm ./lake.txt
else
  tile=$tile_beg
  while [ $tile -le $tile_end ]; do
    outfile=oro.C${res}.tile${tile}.nc
    $APRUN $exe_add_lake ${tile} ${res} ${indir} ${lake_cutoff}
    err=$?
    if [ $err != 0 ]; then
      set +x
      echo ERROR CREATING LAKE FRACTION FOR TILE $tile
      exit $err
    fi
    echo "lake fraction is added to $outfile"
    tile=$(( $tile + 1 ))
  done
fi

exit 0
