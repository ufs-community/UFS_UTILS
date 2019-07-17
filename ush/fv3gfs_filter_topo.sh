#!/bin/ksh
set -ax

if [ $# -ne 9 ]; then
   set +x
   echo
   echo "FATAL ERROR: Usage: $0 resolution grid_dir orog_dir out_dir cd4 peak_fac max_slope n_del2_weak script_dir"
   echo
   set -x
   exit 1
fi

if [ $gtype = stretch ] || [ $gtype = regional ]; then
  stretch=$stretch_fac
else
  stretch=1.0
fi

if [ $gtype = regional ]; then
  refine_ratio=$refine_ratio
else
  refine_ratio=1
fi

export res=$1 
export griddir=$2
export orodir=$3
export outdir=$4
export script_dir=$9

executable=$exec_dir/filter_topo
if [ ! -s $executable ]; then
  set +x
  echo
  echo "FATAL ERROR: ${executable} does not exist"
  echo
  set -x
  exit 1 
fi

export mosaic_grid=C${res}_mosaic.nc
export topo_file=oro.C${res}

if [ ! -s $outdir ]; then mkdir -p $outdir ;fi
cd $outdir ||exit 8

cp $griddir/$mosaic_grid .
cp $griddir/C${res}_grid.tile?.nc .
cp $orodir/${topo_file}.tile?.nc .
cp $executable .

regional=.false.
if [ $gtype = regional ]; then
  regional=.true.
fi

cat > input.nml <<EOF
&filter_topo_nml
  grid_file = $mosaic_grid
  topo_file = $topo_file
  mask_field = "land_frac"    ! Defaults:
  cd4 = $5                    ! 0.15
  peak_fac =  $6              ! 1.0
  max_slope = $7              ! 0.12
  n_del2_weak = $8            ! 16
  n_del2_weak = $8            ! 16
  regional = $regional 
  stretch_fac = $stretch
  refine_ratio = $refine_ratio
  res = $res
  /
EOF

$APRUN $executable

if [ $? -ne 0 ]; then
  set +x
  echo
  echo "FATAL ERROR running filter topography for C$res "
  echo
  set -x
  exit 1
else
  set +x
  echo
  echo "Successfully ran filter topography for C$res"
  echo
  set -x
  exit 0
fi

exit
