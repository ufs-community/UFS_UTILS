#!/bin/bash
################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         fv3gfs_driver_grid.sh
# Script description:  Driver script to create cubic-sphere model grid
#
# Abstract: This driver script creates a cubic-sphere based model grid
#           and associated fixed fields. Supports:
#             1) Global uniform grid
#             2) Global stretched grid
#             3) Global stretched grid with nest
#             4) Stand-alone GFDL regional grid
#             5) Stand-alone Extended Schmidt Gnomonic (ESG) regional grid
#
#           Produces mosaic/grid files, orography files, optional GSL drag
#           suite orography files, and surface climatology fields.
#
#           Calls the following scripts:
#             1) fv3gfs_make_grid.sh (make grid files)
#             2) fv3gfs_make_lake.sh (add lakes)
#             3) fv3gfs_make_orog.sh (make land mask and terrain)
#             4) fv3gfs_ocean_merge.sh (merge ocean grid for uniform only)
#             5) fv3gfs_make_orog_gsl.sh (make GSL drag orog files)
#             6) fv3gfs_filter_topo.sh (filter topography)
#             7) sfc_climo_gen.sh (create surface climo fields)
#
# Usage:  fv3gfs_driver_grid.sh
#
#   Required Shell Variables:
#     (None - all have defaults)
#
#   Optional Shell Variables:
#     Grid configuration:
#     res               Resolution of tile (e.g., 48, 96, 192, 384, 768, 1152, 3072).
#                       Defaults to 96.
#     gtype             Grid type: 'uniform', 'stretch', 'nest', 'regional_gfdl',
#                       or 'regional_esg'. Defaults to 'uniform'.
#
#     Stretched/nested grid parameters (gtype=stretch, nest, regional_gfdl):
#     stretch_fac       Stretching factor. Defaults to 1.5.
#     target_lon        Center longitude of highest resolution tile. Defaults to -97.5.
#     target_lat        Center latitude of highest resolution tile. Defaults to 35.5.
#     refine_ratio      Refinement ratio (nest/regional_gfdl). Defaults to 3.
#     istart_nest       Starting i-index of nest in parent supergrid. Defaults to 27.
#     jstart_nest       Starting j-index of nest in parent supergrid. Defaults to 37.
#     iend_nest         Ending i-index of nest in parent supergrid. Defaults to 166.
#     jend_nest         Ending j-index of nest in parent supergrid. Defaults to 164.
#     halo              Halo size (regional grids). Defaults to 3.
#
#     ESG regional grid parameters (gtype=regional_esg):
#     idim              Grid dimension in i-direction. Defaults to 200.
#     jdim              Grid dimension in j-direction. Defaults to 200.
#     delx              Grid spacing in degrees (i-direction, supergrid). Defaults to 0.0585.
#     dely              Grid spacing in degrees (j-direction, supergrid). Defaults to 0.0585.
#     pazi              Azimuthal rotation angle. Defaults to 0.
#
#     Processing options:
#     add_lake          Add lake fraction and depth (uniform only). Defaults to false.
#     lake_cutoff       Lake fraction threshold. Defaults to 0.50.
#     binary_lake       Return 1 if lake_frac >= cutoff. Defaults to 1.
#     lake_data_srce    Lake data source. Defaults to "MODISP_GLDBV3".
#     make_gsl_orog     Create GSL drag suite orog files. Defaults to false.
#     vegsoilt_frac     Output fractional vegetation/soil type. Defaults to .false.
#     veg_type_src      Vegetation type source. Defaults to "modis.igbp.0.05".
#     soil_type_src     Soil type source. Defaults to "statsgo.0.05".
#     ocn               Ocean grid resolution (e.g., 025, 050, 100).
#                       When set, uses coupled model orog files.
#
#   Example:
#     export res=96
#     export gtype=uniform
#     ./fv3gfs_driver_grid.sh
#
# Remarks:
#   - This is a driver script typically run by machine-specific drivers
#     in ./driver_scripts
#   - sfc_climo_gen requires MPI task count that is multiple of 6
#   - Large grids may require tasks spread across multiple nodes
#   - For individual component details, see the called scripts
#
# Attributes:
#   Language: POSIX shell
#
################################################################################

set -eux

#----------------------------------------------------------------------------------
# Makes FV3 cubed-sphere grid
#----------------------------------------------------------------------------------


export veg_type_src=${veg_type_src:-modis.igbp.0.05}
export soil_type_src=${soil_type_src:-statsgo.0.05} 
export lake_data_srce=${lake_data_srce:-MODISP_GLDBV3}

export res=${res:-96}           # resolution of tile: 48, 96, 128, 192, 384, 768, 1152, 3072
export gtype=${gtype:-uniform}  # grid type: uniform, stretch, nest, regional_gfdl,
                                # or regional_esg

export add_lake=${add_lake:-false}      # add lake fraction and depth. uniform only.
export lake_cutoff=${lake_cutoff:-0.50} # return 0 if lake_frac <  lake_cutoff & add_lake=T
export binary_lake=${binary_lake:-1}    # return 1 if lake_frac >= lake_cutoff & add_lake=T

export make_gsl_orog=${make_gsl_orog:-false} # when true, create GSL drag suite orog files.
export vegsoilt_frac=${vegsoilt_frac:-.false.}

if [ $gtype = uniform ];  then
 # export ocn=${ocn:-025}
  echo "Creating global uniform grid"
elif [ $gtype = stretch ]; then
  export stretch_fac=${stretch_fac:-1.5}  # Stretching factor for the grid
  export target_lon=${target_lon:--97.5}  # Center longitude of the highest resolution tile
  export target_lat=${target_lat:-35.5}   # Center latitude of the highest resolution tile
  title=c${res}s
  echo "Creating global stretched grid"
elif [ $gtype = nest ] || [ $gtype = regional_gfdl ]; then
  export stretch_fac=${stretch_fac:-1.5}  # Stretching factor for the grid
  export target_lon=${target_lon:--97.5}  # Center longitude of the highest resolution tile
  export target_lat=${target_lat:-35.5}   # Center latitude of the highest resolution tile
  export refine_ratio=${refine_ratio:-3}  # The refinement ratio
  export istart_nest=${istart_nest:-27}   # Starting i-direction index of nest grid in parent tile supergrid
  export jstart_nest=${jstart_nest:-37}   # Starting j-direction index of nest grid in parent tile supergrid
  export iend_nest=${iend_nest:-166}      # Ending i-direction index of nest grid in parent tile supergrid
  export jend_nest=${jend_nest:-164}      # Ending j-direction index of nest grid in parent tile supergrid
  export halo=${halo:-3}                  # Halo size. Regional grids only.
  title=c${res}s
  if [ $gtype = nest ];then
   echo "Creating global nested grid"
  else
   echo "Creating gfdl regional grid"
  fi
elif [ $gtype = regional_esg ]; then
  echo "Creating esg regional grid"
  export target_lon=${target_lon:--97.5}  # Center longitude of grid
  export target_lat=${target_lat:-35.5}   # Center latitude of grid
  export idim=${idim:-200}                # Dimension of grid in 'i' direction
  export jdim=${jdim:-200}                # Dimension of grid in 'j' direction
  export delx=${delx:-0.0585}             # Grid spacing (in degrees) in the 'i' direction
                                          # on the SUPERGRID (which has twice the resolution of
                                          # the model grid).  The physical grid spacing in the 'i'
                                          # direction is related to delx as follows:
                                          #    distance = 2*delx*(circumf_Earth/360 deg)
  export dely=${dely:-0.0585}             # Grid spacing (in degrees) in the 'j' direction.
  export halo=${halo:-3}                  # Number of rows/cols for halo.
  title=esg
else
  echo "Error: please specify grid type with 'gtype' as uniform, stretch, nest, regional_gfdl or regional_esg"
  exit 9
fi

export TEMP_DIR=${TEMP_DIR:?}
export out_dir=${out_dir:?}
export home_dir=${home_dir:-"$PWD/../"}
export script_dir=$home_dir/ush
export exec_dir=${exec_dir:-"$home_dir/exec"}
export topo=$home_dir/fix/orog
export NCDUMP=${NCDUMP:-ncdump}


rm -fr $TEMP_DIR
mkdir -p $TEMP_DIR
cd $TEMP_DIR ||exit 8

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Make grid and orography.
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Uniform, stretch or nest grid.
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

if [ $gtype = uniform ] || [ $gtype = stretch ] || [ $gtype = nest ];  then

  if [ $gtype = uniform ] ; then
    export ntiles=6
    name=C${res}
  elif [ $gtype = stretch ]; then
    export ntiles=6
    rn=$( echo "$stretch_fac * 10" | bc | cut -c1-2 )
    name=C${res}r${rn}_${title}
  elif [ $gtype = nest ]; then
    export ntiles=7
    rn=$( echo "$stretch_fac * 10" | bc | cut -c1-2 )
    name=C${res}r${rn}n${refine_ratio}_${title}
  fi

  export grid_dir=$TEMP_DIR/$name/grid
  export orog_dir=$TEMP_DIR/$name/orog


	#if [ $gtype = uniform ]; then
	if declare -p ocn &>/dev/null;then
			out_dir=$out_dir/C$res.mx$ocn
                	readme_name=readme.C$res.mx$ocn.txt
	else

	out_dir=$out_dir/C$res
        readme_name=readme.C$res.txt
	fi         


  mkdir -p $out_dir
  

  if [ $gtype = nest ]; then
    filter_dir=$orog_dir   # nested grid topography will be filtered online
  else
    filter_dir=$TEMP_DIR/$name/filter_topo
  fi

  rm -rf $TEMP_DIR/$name                  
  mkdir -p $grid_dir $orog_dir $filter_dir

  set +x
  echo 
  echo "............ Execute fv3gfs_make_grid.sh ................."
  echo 
  set -x
  if [ $gtype = nest ]; then
    $script_dir/fv3gfs_make_grid.sh $grid_dir $istart_nest $jstart_nest $iend_nest $jend_nest
  else
    $script_dir/fv3gfs_make_grid.sh $grid_dir
  fi
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi
 
  echo "Begin uniform orography generation at `date`"

  tile=1
  while [ $tile -le $ntiles ]; do
    set +x
    echo
    echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
    echo
    set -x
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $topo
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi
    if [ $make_gsl_orog = true ]; then
      set +x
      echo
      echo "............ Execute fv3gfs_make_orog_gsl.sh for tile $tile .................."
      echo 
      set -x
      export halo_tmp="-999"  # no halo
      $script_dir/fv3gfs_make_orog_gsl.sh $res $tile $halo_tmp $grid_dir $orog_dir $topo
      err=$?
      if [ $err != 0 ]; then
        exit $err
      fi
    fi
    tile=$(( $tile + 1 ))
  done

  if [ $add_lake = true ]; then
    $script_dir/fv3gfs_make_lake.sh
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi
  fi

	if [ $gtype = uniform ]; then
	  if declare -p ocn &>/dev/null;then 
 	     $script_dir/fv3gfs_ocean_merge.sh
	     err=$?
             if [ $err != 0 ]; then
                exit $err
     	     fi
           fi
	fi

  set +x
  echo "End uniform orography generation at `date`"
  set -x

#----------------------------------------------------------------------------------
# Topo filtering for uniform and stretched grids only.
#----------------------------------------------------------------------------------

  if [ $gtype = uniform ] || [ $gtype = stretch ]; then
 
    set +x
    echo 
    echo "............ Execute fv3gfs_filter_topo.sh .............."
    echo
    set -x
    $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi

  fi # run topo filtering

  echo "Copy grid and orography files to output directory"

  tile=1
  while [ $tile -le $ntiles ]; do
	
  	if declare -p ocn &>/dev/null;then
	cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/C${res}.mx${ocn}_oro_data.tile${tile}.nc
   	cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}_grid.tile${tile}.nc
        else
	cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/C${res}_oro_data.tile${tile}.nc
        cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}_grid.tile${tile}.nc
	fi

	 if [ $make_gsl_orog = true ]; then
      		cp $orog_dir/C${res}_oro_data*.tile${tile}*.nc $out_dir/  # gsl drag suite oro_data files
   	 fi
    		tile=`expr $tile + 1 `
  done

  cp $grid_dir/C${res}_*mosaic.nc             $out_dir

  echo "Grid and orography files are now prepared."

#exit 0

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Regional grid (gfdl or esg)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

elif [ $gtype = regional_gfdl ] || [ $gtype = regional_esg ]; then
 
#----------------------------------------------------------------------------------
# We are now creating only 1 tile and it is tile 7
#----------------------------------------------------------------------------------
 
  export ntiles=1
  halop1=$(( halo + 1 ))
  tile=7
  name=regional
  export grid_dir=$TEMP_DIR/${name}/grid
  export orog_dir=$TEMP_DIR/${name}/orog
  filter_dir=$TEMP_DIR/$name/filter_topo
  rm -rf $TEMP_DIR/$name
  mkdir -p $grid_dir $orog_dir $filter_dir
  readme_name=readme.$gtype.txt

#----------------------------------------------------------------------------------
# Create regional gfdl grid files.
#----------------------------------------------------------------------------------

  if [ $gtype = regional_gfdl ]; then

    set +x # don't echo all the computation to figure out how many points to add/subtract from start/end nest values
 
    nptsx=`expr $iend_nest - $istart_nest + 1`  # parent points
    nptsy=`expr $jend_nest - $jstart_nest + 1`
 
    idim=`expr $nptsx  \* $refine_ratio / 2`    # number of compute points
    jdim=`expr $nptsy  \* $refine_ratio / 2`
 
#----------------------------------------------------------------------------------
# Figure out how many columns/rows to add in each direction so we have at least 
# 5 halo points for make_hgrid and the orography program.
#----------------------------------------------------------------------------------
 
    index=0
    add_subtract_value=0
    while (test "$index" -le "0")
     do
      add_subtract_value=`expr $add_subtract_value + 1`
      iend_nest_halo=`expr $iend_nest + $add_subtract_value`
      istart_nest_halo=`expr $istart_nest - $add_subtract_value`
      newpoints_i=`expr $iend_nest_halo - $istart_nest_halo + 1`
      newpoints_cg_i=`expr $newpoints_i  \* $refine_ratio / 2`
      diff=`expr $newpoints_cg_i - $idim`
      if [ $diff -ge 10 ]; then 
       index=`expr $index + 1`
      fi
     done
    jend_nest_halo=`expr $jend_nest + $add_subtract_value`
    jstart_nest_halo=`expr $jstart_nest - $add_subtract_value`

    echo "================================================================================== "
    echo "For refine_ratio= $refine_ratio" 
    echo " iend_nest= $iend_nest iend_nest_halo= $iend_nest_halo istart_nest= $istart_nest istart_nest_halo= $istart_nest_halo"
    echo " jend_nest= $jend_nest jend_nest_halo= $jend_nest_halo jstart_nest= $jstart_nest jstart_nest_halo= $jstart_nest_halo"
    echo "================================================================================== "
 
    set +x
    echo
    echo "............ Execute fv3gfs_make_grid.sh ................."
    echo
    set -x
    $script_dir/fv3gfs_make_grid.sh $grid_dir $istart_nest_halo $jstart_nest_halo $iend_nest_halo $jend_nest_halo
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi

#----------------------------------------------------------------------------------
# Create regional esg grid files.
#----------------------------------------------------------------------------------

  elif [ $gtype = regional_esg ]; then

    set +x
    echo
    echo "............ Execute fv3gfs_make_grid.sh ................."
    echo
    set -x
    $script_dir/fv3gfs_make_grid.sh $grid_dir
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi

  fi

#----------------------------------------------------------------------------------
# Redefine resolution for regional grids as a global equivalent resolution.
#----------------------------------------------------------------------------------

  res=$( $NCDUMP -h ${grid_dir}/C*_grid.tile7.nc | grep -o ":RES_equiv = [0-9]\+" | grep -o "[0-9]" )
  res=${res//$'\n'/}
  out_dir=$out_dir/C${res}
  mkdir -p $out_dir

#----------------------------------------------------------------------------------
# Create orography.
#----------------------------------------------------------------------------------
 
  echo "Begin orography generation at `date`"

  set +x
  echo
  echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
  echo
  set -x
  $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $topo
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi

# add lake data to the orography file, if $add_lake is true
 
  if [ $add_lake = true ]; then
    $script_dir/fv3gfs_make_lake.sh
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi
  fi

  echo "Grid and orography files are now prepared."

  set +x
  echo
  echo "............ Execute  fv3gfs_filter_topo.sh .............."
  echo
  set -x
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi

#----------------------------------------------------------------------------------
# For regional grids, shave the orography file and then the grid file, the echo 
# creates the file that contains the number of required points in x and y and the 
# input and output file names.This first run of shave uses a halo of 4.
# This is necessary so that chgres will create BC's with 4 rows/columns which is 
# necessary for pt.
#----------------------------------------------------------------------------------

  set +x
  echo
  echo "............ Execute shave to reduce grid and orography files to required compute size .............."
  echo
  set -x

  cd $filter_dir

  echo $idim $jdim $halop1 \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog
  echo $idim $jdim $halop1 \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid

  $APRUN $exec_dir/shave <input.shave.orog
  $APRUN $exec_dir/shave <input.shave.grid

  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc   $out_dir/C${res}_oro_data.tile${tile}.halo${halop1}.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo${halop1}.nc
 
#----------------------------------------------------------------------------------
# Now shave the orography file and then the grid file with a halo of 3. 
# This is necessary for running the model.
#----------------------------------------------------------------------------------

  echo $idim $jdim $halo \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog.halo$halo
  echo $idim $jdim $halo \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid.halo$halo

  $APRUN $exec_dir/shave <input.shave.orog.halo$halo
  $APRUN $exec_dir/shave <input.shave.grid.halo$halo
 
  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc $out_dir/C${res}_oro_data.tile${tile}.halo${halo}.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo${halo}.nc
 
#----------------------------------------------------------------------------------
# Now shave the orography file and then the grid file with a halo of 0. 
# This is handy for running chgres.
#----------------------------------------------------------------------------------

  echo $idim $jdim 0 \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog.halo0
  echo $idim $jdim 0 \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid.halo0

  $APRUN $exec_dir/shave <input.shave.orog.halo0
  $APRUN $exec_dir/shave <input.shave.grid.halo0

  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc   $out_dir/C${res}_oro_data.tile${tile}.halo0.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo0.nc
 
  cp $grid_dir/C${res}_*mosaic.nc                   $out_dir


#----------------------------------------------------------------------------------
# Now that C${res}_grid.tile${tile}.halo0.nc has been created, we can use it
# to generate gsl drag suite oro_data files, which are generated only for halo0
# Note:  This is carried out only if $make_gsl_orog = true
#----------------------------------------------------------------------------------

  if [ $make_gsl_orog = true ]; then
    export halo_tmp="0"
    ln -sf $out_dir/C${res}_grid.tile${tile}.halo0.nc $grid_dir/
    set +x 
    echo
    echo "............ Execute fv3gfs_make_orog_gsl.sh for tile $tile .................."
    echo
    set -x
    $script_dir/fv3gfs_make_orog_gsl.sh $res $tile $halo_tmp $grid_dir $orog_dir $topo
    err=$?
    if [ $err != 0 ]; then
      exit $err
    fi
    cp $orog_dir/C${res}_oro_data_*.tile${tile}*.nc $out_dir/  # gsl drag suite oro_data files

  fi

  echo "Grid and orography files are now prepared for regional grid"

fi

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Create surface static fields - vegetation type, soil type, etc.
#
# For global grids with a nest, the program is run twice.  First
# to create the fields for the six global tiles.  Then to create
# the fields on the high-res nest.  This is done because the
# ESMF libraries can not interpolate to seven tiles at once.
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

export WORK_DIR=$TEMP_DIR/sfcfields
export SAVE_DIR=$out_dir/sfc
export BASE_DIR=$home_dir
export FIX_FV3=$out_dir
export input_sfc_climo_dir=$home_dir/fix/sfc_climo


if [ $gtype = regional_gfdl ] || [ $gtype = regional_esg ]; then
  export HALO=$halop1
  ln -fs $out_dir/C${res}_grid.tile${tile}.halo${HALO}.nc $out_dir/C${res}_grid.tile${tile}.nc
  ln -fs $out_dir/C${res}_oro_data.tile${tile}.halo${HALO}.nc $out_dir/C${res}_oro_data.tile${tile}.nc

  export GRIDTYPE=regional
elif [ $gtype = nest ]; then
  export mosaic_file=$out_dir/C${res}_coarse_mosaic.nc
fi

$script_dir/sfc_climo_gen.sh
err=$?
if [ $err != 0 ]; then
  echo error in sfc_climo_gen
  exit $err
fi

if [ $gtype = regional_gfdl ] || [ $gtype = regional_esg ]; then
  rm -f $out_dir/C${res}_grid.tile${tile}.nc
  rm -f $out_dir/C${res}_oro_data.tile${tile}.nc
fi

#------------------------------------------------------------------------------------
# Run for the global nest - tile 7.
#------------------------------------------------------------------------------------

if [ $gtype = nest ]; then
  export mosaic_file=$out_dir/C${res}_nested_mosaic.nc
  export GRIDTYPE=nest
  $script_dir/sfc_climo_gen.sh
  err=$?
  if [ $err != 0 ]; then
    echo error in sfc_climo_gen
    exit $err
  fi
fi



#------------------------------------------------------------------------------------
# Make the README files with all relevant info to reproduce the outputs
#------------------------------------------------------------------------------------

cd $home_dir

commit_string=$(git log -1 --oneline)
commit_num=$(echo $commit_string | cut -c1-7)
cd $out_dir

if [ $gtype = uniform ] || [ $gtype = stretch ]; then

cat <<EOF > $readme_name
The following parameters were used
	commit_num=$commit_num
	creation date=$(date +%Y-%m-%d)
        gtype=$gtype
        make_gsl_orog=$make_gsl_orog
        vegsoilt_frac=$vegsoilt_frac
        veg_type=$veg_type_src
        soil_type=$soil_type_src
        add_lake=$add_lake
	lake_data_srce=$lake_data_srce
        binary_lake=$binary_lake
	lake_cutoff=$lake_cutoff
EOF
elif [ $gtype = nest ] || [ $gtype = regional_gfdl ]; then


cat <<EOF > $readme_name
The following parameters were used
        commit_num=$commit_num
	creation date=$(date +%Y-%m-%d)
        gtype=$gtype
        vegsoilt_frac=$vegsoilt_frac
        veg_type=$veg_type_src
        soil_type=$soil_type_src
        make_gsl_orog=$make_gsl_orog
        vegsoilt_frac=$vegsoilt_frac
        veg_type=$veg_type_src
        soil_type=$soil_type_src
        add_lake=$add_lake
	lake_data_srce=$lake_data_srce
        lake_cutoff=$lake_cutoff
        binary_lake=$binary_lake
        stretch_fac=$stretch_fac        # Stretching factor for the grid
        target_lon=$target_lon          # Center longitude of the highest resolution tile
        target_lat=$target_lat          # Center latitude of the highest resolution tile
        refine_ratio=$refine_ratio      # The refinement ratio
        istart_nest=$istart_nest        # Starting i-direction index of nest grid in parent tile supergrid
        jstart_nest=$jstart_nest        # Starting j-direction index of nest grid in parent tile supergrid
        iend_nest=$iend_nest            # Ending i-direction index of nest grid in parent tile supergrid
        jend_nest=$jend_nest            # Ending j-direction index of nest grid in parent tile supergrid
        halo=$halo                      # Lateral boundary halo
EOF
elif [ $gtype = regional_esg ] ; then

cat <<EOF > $readme_name
The following parameters were used
        commit_num=$commit_num
        creation date=$(date +%Y-%m-%d)
	gtype=$gtype
        res=-999                        # equivalent resolution is computed
        vegsoilt_frac=$vegsoilt_frac
        veg_type=$veg_type_src
        soil_type=$soil_type_src
	lake_data_srce=$lake_data_srce
        target_lon=$target_lon          # Center longitude of grid
        target_lat=$target_lat          # Center latitude of grid
        idim=$idim                      # Dimension of grid in 'i' direction
        jdim=$jdim                      # Dimension of grid in 'j' direction
        delx=$delx                      # Grid spacing (in degrees) in the 'i' direction
                                        # on the SUPERGRID (which has twice the resolution of
                                        # the model grid).  The physical grid spacing in the 'i'
                                        # direction is related to delx as follows:
        dely=$dely                      # Grid spacing (in degrees) in the 'j' direction.
        halo=$halo                      # number of row/cols for halo
EOF
fi




exit
