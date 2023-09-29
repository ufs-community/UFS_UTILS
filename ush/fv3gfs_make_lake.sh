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
echo "lake_data_srce = $lake_data_srce"
if [ $lake_data_srce == GLDBV3 ]; then
  lakestatusrc="GLDBV3"
  lakedepthsrc="GLDBV3"
fi
if [ $lake_data_srce == GLDBV2 ]; then
  lakestatusrc="GLDBV2"
  lakedepthsrc="GLDBV2"
fi
if [ $lake_data_srce == MODISP_GLOBATHY ]; then
  lakestatusrc="MODISP"
  lakedepthsrc="GLOBATHY"
fi
if [ $lake_data_srce == VIIRS_GLOBATHY ]; then
  lakestatusrc="VIIRS"
  lakedepthsrc="GLOBATHY"
fi
if [ $lake_data_srce == MODISP_GLDBV3 ]; then
  lakestatusrc="MODISP"
  lakedepthsrc="GLDBV3"
fi
if [ $lake_data_srce == VIIRS_GLDBV3 ]; then
  lakestatusrc="VIIRS"
  lakedepthsrc="GLDBV3"
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

if [ $gtype == regional_gfdl ] || [ $gtype == regional_esg ]; then
  tile_beg=7
  tile_end=7
  tile=7
  outfile=oro.C${res}.tile${tile}.nc
  ln -s $outdir/$outfile .
  outgrid="C${res}_grid.tile${tile}.nc"
  ln -s $grid_dir/$outgrid .
fi

# create inland mask and save it to the orography files

cutoff=0.75
rd=7
if [ $gtype == uniform ]; then
  $APRUN $exe_inland $res $cutoff $rd g
fi
if [ $gtype == regional_gfdl ]; then
  $APRUN $exe_inland $res $cutoff $rd r
fi
if [ $gtype == regional_esg ]; then
  $APRUN $exe_inland $res $cutoff $rd r
fi
err=$?
if [ $err != 0 ]; then
  set +x
  echo ERROR CREATING INLAND MASK
  exit $err
fi

# create fractional lake data for FV3 grid and save it to the orography files

tile=$tile_beg
while [ $tile -le $tile_end ]; do
  outfile=oro.C${res}.tile${tile}.nc
  $APRUN $exe_add_lake ${tile} ${res} ${indir} ${lakestatusrc} ${lakedepthsrc} ${lake_cutoff}
  err=$?
  if [ $err != 0 ]; then
    set +x
    echo ERROR CREATING LAKE FRACTION FOR TILE $tile
    exit $err
  fi
  echo "lake fraction is added to $outfile"
  tile=$(( $tile + 1 ))
done

exit 0
