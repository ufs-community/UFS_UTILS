#!/bin/bash
#
#-----------------------------------------------------------------------
# Driver script to create a cubic-sphere based model grid.
#
# Supports the following grids:
#   1) global uniform
#   2) global stretched
#   3) global stretched with nest
#   4) stand-alone GFDL regional
#   5) stand-alone extended Schmidt gnonomic (ESG) regional
#
# Produces the following files (netcdf, each tile in separate file):
#   1) 'mosaic' and 'grid' files containing lat/lon and other
#      records that describe the model grid.
#   2) 'oro' files containing land mask, terrain and gravity
#      wave drag fields.
#   3) 'oro' files ('oro_data_ls' and 'oro_data_ss') specific to
#      the GSL drag suite physics parameterization (only if
#      flag make_gsl_orog = true)
#   4) surface climo fields, such as soil type, vegetation
#      greenness and albedo.
#
# Calls the following scripts
#   1) fv3gfs_make_grid.sh (make 'grid' files)
#   2) fv3gfs_make_lake.sh (adds lakes)
#   3) fv3gfs_make_orog.sh (make land mask and exits for uniform grid type only, generates oro for other grid types)
#   4) fv3gfs_ocean_merge.sh (Reads pre-generated ocean grid and merges them for uniform grid type only)
#   5) fv3gfs_make_orog.sh (reads merged masks and makes 'oro' files for uniform grid type only)
#   6) fv3gfs_make_orog_gsl.sh (make gsl drag 'oro' files for uniform grid type only)
#   7) fv3gfs_filter_topo.sh (filter topography)
#   8) sfc_climo_gen.sh (create surface climo fields)
# 
#    Note: The sfc_climo_gen program only runs with an
#       mpi task count that is a multiple of six.  This is
#       an ESMF library requirement.  Large grids may require
#       tasks spread across multiple nodes.
#
# This script is run by its machine-specific driver script in
# ./driver_scripts.
#-----------------------------------------------------------------------

set -eux

#----------------------------------------------------------------------------------
# Makes FV3 cubed-sphere grid
#----------------------------------------------------------------------------------

export res=${res:-96}           # resolution of tile: 48, 96, 128, 192, 384, 768, 1152, 3072
export gtype=${gtype:-uniform}  # grid type: uniform, stretch, nest, regional_gfdl,
                                # or regional_esg

export add_lake=${add_lake:-false}      # add lake fraction and depth. uniform only.
export lake_cutoff=${lake_cutoff:-0.50} # return 0 if lake_frac <  lake_cutoff & add_lake=T
export binary_lake=${binary_lake:-1}    # return 1 if lake_frac >= lake_cutoff & add_lake=T

export make_gsl_orog=${make_gsl_orog:-false} # when true, create GSL drag suite orog files.

if [ $gtype = uniform ];  then
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
export topo=$home_dir/fix/orog_raw

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


	if [ $gtype = uniform ]; then
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
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo
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
  		 $script_dir/fv3gfs_ocean_merge.sh
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
	cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/oro_C${res}.mx${ocn}.tile${tile}.nc
   	cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}.mx${ocn}_grid.tile${tile}.nc
        else
	cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/oro_C${res}.tile${tile}.nc
        cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}_grid.tile${tile}.nc
	fi

	 if [ $make_gsl_orog = true ]; then
      		cp $orog_dir/C${res}_oro_data*.tile${tile}*.nc $out_dir/  # gsl drag suite oro_data files
   	 fi
    		tile=`expr $tile + 1 `
  done

  cp $grid_dir/C${res}_*mosaic.nc             $out_dir

  echo "Grid and orography files are now prepared."

#exit 8

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
  $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo
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

cd $out_dir


if [ $gtype = uniform ] || [ $gtype = stretch ]; then

cat <<EOF > $readme_name
The following # was used
https://github.com/sanatcumar/UFS_UTILS/tree/single_step
The following parameters were used
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
The following # was used
https://github.com/sanatcumar/UFS_UTILS/tree/single_step
The following parameters were used
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
The following # was used
https://github.com/sanatcumar/UFS_UTILS/tree/single_step
The following parameters were used
        gtype=$gtype
        res=-999                        # equivalent resolution is computed
        vegsoilt_frac=$vegsoilt_frac
        veg_type=$veg_type_src
        soil_type=$soil_type_src
        target_lon=$target_lon          # Center longitude of grid
        target_lat=target_lat           # Center latitude of grid
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
