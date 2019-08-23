#!/bin/ksh
set -ax

#-----------------------------------------------------------------------------------------
#
# Script name: fv3gfs_make_grid.sh
# -----------
#
# Description: Makes the 'grid' and mosaic file for an fv3 grid.
# -----------
#
# Shell variables:
# ---------------
#
#   APRUN            Command to invoke executables
#   exec_dir         Location of executables
#   gtype            Grid type.  'uniform' - global uniform; 'stretch' - global 
#                    stretched; 'nest' - global stretched with nest; 'regional' -
#                    stand alone regional nest.
#   halo             Lateral boundary halo size, regional grids only.
#   i/jstart_nest    Starting i/j index of nest within parent tile.
#   i/jend_nst       Ending i/j index of nest within parent tile.
#   outdir           Working directory.
#   refine_ratio     Nest refinement ratio.
#   res              Resolution, i.e, 96.
#   stretch_fac      Stretching factor
#   target_lat       Center latitude of highest resolution tile
#   target_lon       Center longitude of highest resolution tile
#
#-----------------------------------------------------------------------------------------

gtype=${gtype:?}
res=${res:?}
outdir=${outdir:-$1}
exec_dir=${exec_dir:?}
APRUN=${APRUN:-time}
nx=`expr $res \* 2 `

if [ ! -s $outdir ]; then  mkdir -p $outdir ;fi
cd $outdir

executable=$exec_dir/make_hgrid
if [ ! -s $executable ]; then
  set +x
  echo
  echo "FATAL ERROR: ${executable} does not exist"
  echo
  set -x
  exit 1 
fi

#-----------------------------------------------------------------------------------------
# Create 'grid' files - one for each tile in netcdf.
#-----------------------------------------------------------------------------------------

if [ $gtype = uniform ]; then
  ntiles=6
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid
elif  [ $gtype = stretch ]; then
  stretch_fac=${stretch_fac:?}
  target_lon=${target_lon:?}
  target_lat=${target_lat:?}
  ntiles=6
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid \
                     --do_schmidt --stretch_factor ${stretch_fac} --target_lon ${target_lon} --target_lat ${target_lat} 
elif  [ $gtype = nest ] || [ $gtype = regional ] ; then
  stretch_fac=${stretch_fac:?}
  target_lon=${target_lon:?}
  target_lat=${target_lat:?}
  refine_ratio=${refine_ratio:?}
  istart_nest=$2
  jstart_nest=$3
  iend_nest=$4
  jend_nest=$5
  halo=${halo:?}
  if  [ $gtype = regional ]; then
    ntiles=1
  else
    ntiles=7
  fi
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid \
                     --do_schmidt --stretch_factor ${stretch_fac} --target_lon ${target_lon} --target_lat ${target_lat} \
                     --nest_grid --parent_tile 6 --refine_ratio $refine_ratio --istart_nest $istart_nest --jstart_nest $jstart_nest \
                     --iend_nest $iend_nest --jend_nest $jend_nest --halo $halo --great_circle_algorithm
fi

if [ $? -ne 0 ]; then
  set +x
  echo
  echo "FATAL ERROR creating C$res grid."
  echo
  set -x
  exit 1
fi

#---------------------------------------------------------------------------------------
# Create mosaic file.
#
# For global nest grids, create three mosaic files - one for the nest, one for the
# global tiles and one for all tiles.  Needed because ESMF-based utilities can't
# process seven tiles at once.
#
# For regional grids, there is only one tile and it is tile number 7.
#---------------------------------------------------------------------------------------

executable=$exec_dir/make_solo_mosaic
if [ ! -s $executable ]; then
  set +x
  echo
  echo "FATAL ERROR: ${executable} does not exist."
  echo
  set -x
  exit 1 
fi

if [ $gtype = uniform ] || [ $gtype = stretch ] ; then

  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file \
     C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc

elif [ $gtype = nest ]; then

  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file \
     C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc,C${res}_grid.tile7.nc    

  $APRUN $executable --num_tiles 6 --dir $outdir --mosaic C${res}_coarse_mosaic --tile_file \
     C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc

  $APRUN $executable --num_tiles 1 --dir $outdir --mosaic C${res}_nested_mosaic --tile_file C${res}_grid.tile7.nc 

elif [ $gtype = regional ];then

  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file C${res}_grid.tile7.nc

fi

if [ $? -ne 0 ]; then
  set +x
  echo
  echo "FATAL ERROR creating mosaic file."
  echo
  set -x
  exit 1
fi

exit
