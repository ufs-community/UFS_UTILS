#!/bin/bash

set -eux

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
#                    stretched; 'nest' - global stretched with nest; 'regional_gfdl' -
#                    stand alone GFDL regional nest; 'regional_esg' - stand alone
#                    extended Schmidt gnomonic regional grid.
#   halo             Lateral boundary halo size, regional grids only.
#   i/jstart_nest    Starting i/j index of nest within parent tile.
#   i/jend_nest      Ending i/j index of nest within parent tile.
#   outdir           Working directory.
#   refine_ratio     Nest refinement ratio.
#   res              Global resolution, i.e, 96.  For regional grids, this
#                    redefined as a global equivalent resolution.
#   stretch_fac      Stretching factor
#   target_lat       Center latitude of highest resolution tile
#   target_lon       Center longitude of highest resolution tile
#
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Retrieve global equivalent resolution from grid file.
#-----------------------------------------------------------------------------------------

function get_res
{
  res=$( $NCDUMP -h $1 | grep -o ":RES_equiv = [0-9]\+" | grep -o "[0-9]" )
  res=${res//$'\n'/}
}

#-----------------------------------------------------------------------------------------
# Main script begins here.
#-----------------------------------------------------------------------------------------

gtype=${gtype:?}
res=${res:?}
outdir=${outdir:-$1}
exec_dir=${exec_dir:?}
APRUN=${APRUN:-time}
nx=`expr $res \* 2 `

if [ ! -s $outdir ]; then  mkdir -p $outdir ;fi
cd $outdir

#-----------------------------------------------------------------------------------------
# Create 'grid' files - one for each tile in netcdf.
#-----------------------------------------------------------------------------------------

if [ $gtype = regional_esg ]; then
  executable=$exec_dir/regional_esg_grid
else
  executable=$exec_dir/make_hgrid
fi

if [ ! -s $executable ]; then
  set +x
  echo
  echo "FATAL ERROR: ${executable} does not exist"
  echo
  set -x
  exit 1 
fi

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
elif  [ $gtype = nest ] || [ $gtype = regional_gfdl ] ; then
  stretch_fac=${stretch_fac:?}
  target_lon=${target_lon:?}
  target_lat=${target_lat:?}
  refine_ratio=${refine_ratio:?}
  istart_nest=$2
  jstart_nest=$3
  iend_nest=$4
  jend_nest=$5
  halo=${halo:?}
  if  [ $gtype = regional_gfdl ]; then
    ntiles=1
  else
    ntiles=7
  fi
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid \
                     --do_schmidt --stretch_factor ${stretch_fac} --target_lon ${target_lon} --target_lat ${target_lat} \
                     --nest_grid --parent_tile 6 --refine_ratio $refine_ratio --istart_nest $istart_nest --jstart_nest $jstart_nest \
                     --iend_nest $iend_nest --jend_nest $jend_nest --halo $halo --great_circle_algorithm

elif [ $gtype = regional_esg ] ; then

  (( halop2=halo+2 ))
  (( lx=idim+halop2*2 ))
  (( ly=jdim+halop2*2 ))

  cat > ./regional_grid.nml << EOF
    &regional_grid_nml
      plon = ${target_lon}
      plat = ${target_lat}
      delx = ${delx}
      dely = ${dely}
      lx   = -${lx}
      ly   = -${ly}
    /
EOF

  $APRUN $executable

fi

if [ $? -ne 0 ]; then
  set +x
  echo
  echo "FATAL ERROR creating grid files."
  echo
  set -x
  exit 1
fi

#---------------------------------------------------------------------------------------
# Compute the equivalent 'global' resolution for regional grids.  For GFDL
# regional grids, the CRES value input to this script is the resolution of the
# global grid in which it is embedded.  Compute a more realistic value.
# Program adds the equivalent resolution as a global attribute to the grid file.
# Below, this attribute is retrieved and used to rename the grid files.
#---------------------------------------------------------------------------------------

if [ $gtype = regional_gfdl ] || [ $gtype = regional_esg ]; then
  executable=$exec_dir/global_equiv_resol
  if [ ! -s $executable ]; then
    set +x
    echo
    echo "FATAL ERROR: ${executable} does not exist."
    echo
    set -x
    exit 1 
  fi
  if [ $gtype = regional_esg ]; then
    $APRUN $executable  regional_grid.nc
  elif [ $gtype = regional_gfdl ]; then
    $APRUN $executable  C${res}_grid.tile7.nc
  fi
  if [ $? -ne 0 ]; then
    set +x
    echo
    echo "FATAL ERROR running global_equiv_resol."
    echo
    set -x
    exit 2
  fi
fi

#---------------------------------------------------------------------------------------
# Create mosaic file.
#
# For global nest grids, create three mosaic files - one for the nest, one for the
# global tiles and one for all tiles.  Needed because ESMF-based utilities can't
# process seven tiles at once.
#
# For regional grids, there is only one tile and it is tile number 7.
#
# For regional grids, rename the grid files according to their global equivalent
# resolution.
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

elif [ $gtype = regional_gfdl ];then

  res_save=$res
  get_res C${res}_grid.tile7.nc
  mv C${res_save}_grid.tile7.nc  C${res}_grid.tile7.nc
  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file C${res}_grid.tile7.nc

elif [ $gtype = regional_esg ]; then

  get_res regional_grid.nc
  mv regional_grid.nc C${res}_grid.tile7.nc
  $APRUN $executable --num_tiles 1 --dir $outdir --mosaic C${res}_mosaic --tile_file C${res}_grid.tile7.nc

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
