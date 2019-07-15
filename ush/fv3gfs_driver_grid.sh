#!/bin/sh
#
#-----------------------------------------------------------------------
# Driver script to create a cubic-sphere based model grid.
#
# Supports the following grids:
#   1) global uniform
#   2) global stretched
#   3) global stretched with nest
#   4) stand-alone regional
#
# Produces the following files (netcdf, each tile in separate file):
#   1) 'mosaic' and 'grid' files containing lat/lon and other
#      records that describe the model grid.
#   2) 'oro' files containing land mask, terrain and gravity
#      wave drag fields.
#   3) surface climo fields, such as soil type, vegetation
#      greenness and albedo.
#
# Note: The sfc_climo_gen program only runs with an
#       mpi task count that is a multiple of six.  This is
#       an ESMF library requirement.  Large grids may require
#       tasks spread across multiple nodes.
#
# To run, do the following:
#   1) Uncomment workload management directives for your machine.
#      Script runs on Theia, WCOSS-Cray and WCOSS-Dell.
#   2) Set "machine" - choices are "THEIA", "WCOSS_C"
#      and "WCOSS_DELL_P3".
#   3) Set "C" resolution, "res" - Example: res=96.
#   4) Set grid type ("gtype").  Valid choices are
#         "uniform"  - global uniform grid
#         "stretch"  - global stretched grid
#         "nest"     - global stretched grid with nest
#         "regional" - stand-alone regional grid
#   5) For "stretch" and "nest" grids, set the stretching factor -
#       "stretch_fac", and center lat/lon of highest resolution
#      tile - "target_lat" and "target_lon".
#   6) For "nest" grids, set the refinement ratio - "refine_ratio", 
#      the starting/ending i/j index location within the parent
#      tile - "istart_nest", "jstart_nest", "iend_nest", "jend_nest"
#   7) For "regional" grids, set the "halo".  Default is three
#      rows/columns.
#   8) Submit script.  On Theia: "sbatch $script".  On Cray and Dell:
#      "cat $script | bsub".
#   9) All files will be placed in "out_dir".
#   
#-----------------------------------------------------------------------

#---- WCOSS DELL JOBCARD
#---- Submit script as "cat $script | bsub"
#BSUB -oo log.grid.%J
#BSUB -eo log.grid.%J
#BSUB -q debug
#BSUB -P FV3GFS-T2O
#BSUB -J grid_fv3
#BSUB -W 0:30
#BSUB -x                 # run not shared
#BSUB -n 24              # total tasks
#BSUB -R span[ptile=24]   # tasks per node
#BSUB -R affinity[core(1):distribute=balance]
machine=WCOSS_DELL_P3

#---- WCOSS_CRAY JOBCARD
#---- Submit script as "cat $script | bsub"
##BSUB -L /bin/sh
##BSUB -P FV3GFS-T2O
##BSUB -oo log.grid.%J
##BSUB -eo log.grid.%J
##BSUB -J grid_fv3
##BSUB -q debug
##BSUB -M 2400
##BSUB -W 00:30
##BSUB -extsched 'CRAYLINUX[]'
#machine=WCOSS_C

#---- THEIA JOBCARD
#---- Submit script as 'sbatch $script'
##SBATCH -J fv3_grid_driver
##SBATCH -A fv3-cpu
##SBATCH --open-mode=truncate
##SBATCH -o log.fv3_grid_driver
##SBATCH -e log.fv3_grid_driver
##SBATCH --nodes=1 --ntasks-per-node=24
##SBATCH -q debug
##SBATCH -t 00:30:00
#machine=THEIA

set -ax

export machine=${machine:?}

#----------------------------------------------------------------------------------
# Makes FV3 cubed-sphere grid
#----------------------------------------------------------------------------------

export res=96              # resolution of tile: 48, 96, 128, 192, 384, 768, 1152, 3072
export gtype=uniform       # grid type: uniform, stretch, nest or regional

if [ $gtype = uniform ];  then
  echo "Creating uniform ICs"
elif [ $gtype = stretch ]; then
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
  title=c${res}s               # Identifier based on refined location
  echo "Creating stretched grid"
elif [ $gtype = nest ] || [ $gtype = regional ]; then
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
  export refine_ratio=3        # The refinement ratio
  export istart_nest=27        # Starting i-direction index of nest grid in parent tile supergrid
  export jstart_nest=37        # Starting j-direction index of nest grid in parent tile supergrid
  export iend_nest=166         # Ending i-direction index of nest grid in parent tile supergrid
  export jend_nest=164         # Ending j-direction index of nest grid in parent tile supergrid
  export halo=3                # Halo size. Regional grids only.
  title=c${res}s               # Identifier based on nest location
  if [ $gtype = nest ];then
   echo "Creating global nested grid"
  else
   echo "Creating regional grid"
  fi
else
  echo "Error: please specify grid type with 'gtype' as uniform, stretch, nest or regional"
  exit 9
fi

if [ $machine = WCOSS_C ]; then
  set +x
  . $MODULESHOME/init/sh
  module load PrgEnv-intel cfp-intel-sandybridge/1.1.0
  module list
  set -x
  export NODES=1
  export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
  export APRUN_SFC="aprun -j 1 -n 24 -N 24"
  export home_dir=$LS_SUBCWD/..
  export TMPDIR=/gpfs/hps3/stmp/$LOGNAME/fv3_grid.$gtype
  export out_dir=/gpfs/hps3/stmp/$LOGNAME/C${res}
#  On Cray, the orography code is optimized for six threads.
#  Do not change.
  export OMP_NUM_THREADS=6
  export OMP_STACKSIZE=2048m
  export KMP_AFFINITY=disabled
elif [ $machine = WCOSS_DELL_P3 ]; then
# On dell, load the EnvVar module before setting any environment variables.
  set +x
  module purge
  module load EnvVars/1.0.2
  module load lsf/10.1
  module load ips/18.0.1.163
  module load impi/18.0.1
  module load NetCDF/4.5.0
  module load HDF5-serial/1.10.1
  module list
  set -x
  export APRUN=time
  export APRUN_SFC="mpirun -l"
  export home_dir=$LS_SUBCWD/..
  export TMPDIR=/gpfs/dell1/stmp/$LOGNAME/fv3_grid.$gtype
  export out_dir=/gpfs/dell1/stmp/$LOGNAME/C${res}
  export OMP_NUM_THREADS=24 # orog code worked best with 24 threads.
  export OMP_STACKSIZE=2048m
elif [ $machine = THEIA ]; then
  . /apps/lmod/lmod/init/sh
  set +x
  module purge
  module load intel/16.1.150
  module load impi
  module load hdf5/1.8.14
  module load netcdf/4.3.0
  module list
  set -x
  export APRUN=time
  export APRUN_SFC=srun
  export home_dir=$SLURM_SUBMIT_DIR/..
  export TMPDIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/fv3_grid.$gtype
  export out_dir=/scratch3/NCEPDEV/stmp1/$LOGNAME/C${res}
  export OMP_NUM_THREADS=24
  export OMP_STACKSIZE=2048m
else
  echo FATAL ERROR. $machine not supported.  Exiting.
  exit 1
fi

export script_dir=$home_dir/ush
export exec_dir=$home_dir/exec
export topo=$home_dir/fix/fix_orog

ulimit -a
ulimit -s unlimited

rm -fr $TMPDIR
mkdir -p $out_dir $TMPDIR
cd $TMPDIR ||exit 8

#----------------------------------------------------------------------------------------
# filter_topo parameters. C192->50km, C384->25km, C768->13km, C1152->8.5km, C3072->3.2km
#----------------------------------------------------------------------------------------

if [ $res -eq 48 ]; then 
  cd4=0.12;  max_slope=0.12; n_del2_weak=4; peak_fac=1.1  
elif [ $res -eq 96 ]; then 
  cd4=0.12;  max_slope=0.12; n_del2_weak=8; peak_fac=1.1  
elif [ $res -eq 128 ]; then
  cd4=0.13;  max_slope=0.12; n_del2_weak=8;  peak_fac=1.1
elif [ $res -eq 192 ]; then 
  cd4=0.15;  max_slope=0.12; n_del2_weak=12; peak_fac=1.05  
elif [ $res -eq 384 ]; then 
  cd4=0.15;  max_slope=0.12; n_del2_weak=12; peak_fac=1.0  
elif [ $res -eq 768 ]; then 
  cd4=0.15;  max_slope=0.12; n_del2_weak=16; peak_fac=1.0  
elif [ $res -eq 1152 ]; then 
  cd4=0.15;  max_slope=0.16; n_del2_weak=20; peak_fac=1.0  
elif [ $res -eq 3072 ]; then 
  cd4=0.15;  max_slope=0.30; n_del2_weak=24; peak_fac=1.0  
else
 echo "grid C$res not supported, exit"
 exit 2
fi

#----------------------------------------------------------------------------------
# Make grid and orography.
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# Uniform grid.
#----------------------------------------------------------------------------------

if [ $gtype = uniform ];  then

  export ntiles=6
  name=C${res}
  grid_dir=$TMPDIR/$name/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$TMPDIR/$name/filter_topo
  rm -rf $TMPDIR/$name                  
  mkdir -p $grid_dir $orog_dir $filter_dir

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
 
  echo "Begin uniform orography generation at `date`"

#----------------------------------------------------------------------------------
# On WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#----------------------------------------------------------------------------------

  if [ $machine = WCOSS_C ]; then
    echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
    rm $TMPDIR/orog.file1
  elif [ $machine = THEIA ] || [ $machine = WCOSS_DELL_P3 ]; then
    for tile in 1 2 3 4 5 6 ; do
      set +x
      echo
      echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
      echo
      set -x
      $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
    done
  fi

  set +x
  echo "End uniform orography generation at `date`"
 
  echo 
  echo "............ Execute fv3gfs_filter_topo.sh .............."
  echo
  set -x
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir $cd4 $peak_fac $max_slope $n_del2_weak $script_dir
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi
  echo "Grid and orography files are now prepared for uniform grid"

#----------------------------------------------------------------------------------
# Stretched grid.
#----------------------------------------------------------------------------------

elif [ $gtype = stretch ]; then

  export ntiles=6
  rn=$( echo "$stretch_fac * 10" | bc | cut -c1-2 )
  name=C${res}r${rn}_${title}
  grid_dir=$TMPDIR/${name}/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$TMPDIR/${name}/filter_topo
  rm -rf $TMPDIR/$name                  
  mkdir -p $grid_dir $orog_dir $filter_dir

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
 
  echo "Begin stretch orography generation at `date`"

#----------------------------------------------------------------------------------
# On WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#----------------------------------------------------------------------------------

  if [ $machine = WCOSS_C ]; then
    echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
    rm $TMPDIR/orog.file1
  elif [ $machine = THEIA ] || [ $machine = WCOSS_DELL_P3 ]; then
    for tile in 1 2 3 4 5 6 ; do
      set +x
      echo
      echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
      echo
      set -x
      $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
    done
  fi

  set +x
  echo "End stretch orography generation at `date`"
 
  echo 
  echo "............ Execute fv3gfs_filter_topo.sh .............."
  echo
  set -x
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir $cd4 $peak_fac $max_slope $n_del2_weak $script_dir
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi
  echo "Grid and orography files are now prepared for stretched grid"

#----------------------------------------------------------------------------------
# Nested grid.
#----------------------------------------------------------------------------------

elif [ $gtype = nest ]; then

  export ntiles=7
  rn=$( echo "$stretch_fac * 10" | bc | cut -c1-2 )
  name=C${res}r${rn}n${refine_ratio}_${title}
  grid_dir=$TMPDIR/${name}/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$orog_dir   # nested grid topography will be filtered online
  rm -rf $TMPDIR/$name                  
  mkdir -p $grid_dir $orog_dir $filter_dir

  set +x
  echo 
  echo "............ Execute fv3gfs_make_grid.sh ................."
  echo
  set -x
  $script_dir/fv3gfs_make_grid.sh $grid_dir $istart_nest $jstart_nest $iend_nest $jend_nest
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi

  echo "Begin stretch nest orography generation at `date`"
 
#----------------------------------------------------------------------------------
# On WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#----------------------------------------------------------------------------------
 
  if [ $machine = WCOSS_C ]; then
    echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    echo "$script_dir/fv3gfs_make_orog.sh $res 7 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
    rm $TMPDIR/orog.file1
  elif [ $machine = THEIA ] || [ $machine = WCOSS_DELL_P3 ]; then
    for tile in 1 2 3 4 5 6 7; do
      set +x
      echo
      echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
      echo
      set -x
      $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
    done
  fi

  echo "Grid and orography files are now prepared for nested grid"

#----------------------------------------------------------------------------------
# Regional grid.
#----------------------------------------------------------------------------------

elif [ $gtype = regional ]; then
 
#----------------------------------------------------------------------------------
# We are now creating only 1 tile and it is tile 7
#----------------------------------------------------------------------------------
 
  export ntiles=1
  halop1=$(( halo + 1 ))
  tile=7
  set +x # don't echo all the computation to figure out how many points to add/subtract from start/end nest values
 
#----------------------------------------------------------------------------------
# Number of parent points
#----------------------------------------------------------------------------------
 
  nptsx=`expr $iend_nest - $istart_nest + 1`
  nptsy=`expr $jend_nest - $jstart_nest + 1`
 
#----------------------------------------------------------------------------------
# Number of compute grid points
#----------------------------------------------------------------------------------
 
  npts_cgx=`expr $nptsx  \* $refine_ratio / 2`
  npts_cgy=`expr $nptsy  \* $refine_ratio / 2`
 
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
    diff=`expr $newpoints_cg_i - $npts_cgx`
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
  set -x
 
  export ntiles=1
  tile=7
  rn=$( echo "$stretch_fac * 10" | bc | cut -c1-2 )
  name=C${res}r${rn}n${refine_ratio}_${title}
  grid_dir=$TMPDIR/${name}/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$orog_dir   # nested grid topography will be filtered online
  rm -rf $TMPDIR/$name
  mkdir -p $grid_dir $orog_dir $filter_dir

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

  echo "Begin regional orography generation at `date`"
 
#----------------------------------------------------------------------------------
# On WCOSS_C use cfp to run multiple tiles simulatneously for the orography.
# For now we only have one tile but in the future we will have more.
#----------------------------------------------------------------------------------
 
  if [ $machine = WCOSS_C ]; then
    echo "$script_dir/fv3gfs_make_orog.sh $res 7 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
    aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
    rm $TMPDIR/orog.file1
  elif [ $machine = THEIA ] || [ $machine = WCOSS_DELL_P3 ]; then
    set +x
    echo
    echo "............ Execute fv3gfs_make_orog.sh for tile $tile .................."
    echo
    set -x
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
  fi

  set +x
  echo
  echo "............ Execute fv3gfs_filter_topo.sh .............."
  echo
  set -x
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir $cd4 $peak_fac $max_slope $n_del2_weak $script_dir
  err=$?
  if [ $err != 0 ]; then
    exit $err
  fi
  set +x
  echo
  echo "............ Execute shave to reduce grid and orography files to required compute size .............."
  echo
  set -x
  cd $filter_dir

#----------------------------------------------------------------------------------
# Shave the orography file and then the grid file, the echo creates the input 
# file that contains the number of required points in x and y and the input
# and output file names.This first run of shave uses a halo of 4.
# This is necessary so that chgres will create BC's with 4 rows/columns which is 
# necessary for pt.
#----------------------------------------------------------------------------------

  echo $npts_cgx $npts_cgy $halop1 \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog
  echo $npts_cgx $npts_cgy $halop1 \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid

  $APRUN $exec_dir/shave.x <input.shave.orog
  $APRUN $exec_dir/shave.x <input.shave.grid

  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc   $out_dir/C${res}_oro_data.tile${tile}.halo${halop1}.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo${halop1}.nc
#
# Now shave the orography file and then the grid file with a halo of 3. This is necessary for running the model.
#
  echo $npts_cgx $npts_cgy $halo \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog.halo$halo
  echo $npts_cgx $npts_cgy $halo \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid.halo$halo

  $APRUN $exec_dir/shave.x <input.shave.orog.halo$halo
  $APRUN $exec_dir/shave.x <input.shave.grid.halo$halo
 
  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc $out_dir/C${res}_oro_data.tile${tile}.halo${halo}.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo${halo}.nc
#
# Now shave the orography file and then the grid file with a halo of 0. This is handy for running chgres.
#
  echo $npts_cgx $npts_cgy 0 \'$filter_dir/oro.C${res}.tile${tile}.nc\' \'$filter_dir/oro.C${res}.tile${tile}.shave.nc\' >input.shave.orog.halo0
  echo $npts_cgx $npts_cgy 0 \'$filter_dir/C${res}_grid.tile${tile}.nc\' \'$filter_dir/C${res}_grid.tile${tile}.shave.nc\' >input.shave.grid.halo0

  $APRUN $exec_dir/shave.x <input.shave.orog.halo0
  $APRUN $exec_dir/shave.x <input.shave.grid.halo0

  cp $filter_dir/oro.C${res}.tile${tile}.shave.nc   $out_dir/C${res}_oro_data.tile${tile}.halo0.nc
  cp $filter_dir/C${res}_grid.tile${tile}.shave.nc  $out_dir/C${res}_grid.tile${tile}.halo0.nc
 
  echo "Grid and orography files are now prepared for regional grid"

fi # if check on type of grid

#----------------------------------------------------------------------------------
# For non-regional grids, copy grid and orography files to output directory.
#----------------------------------------------------------------------------------

if [ $gtype != regional ]; then

  echo "Copy grid and orography files to output directory"

  tile=1
  while [ $tile -le $ntiles ]; do
    cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/C${res}_oro_data.tile${tile}.nc
    cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}_grid.tile${tile}.nc
    tile=`expr $tile + 1 `
  done

fi

#----------------------------------------------------------------------------------
# Copy mosaic file(s) to output directory.
#----------------------------------------------------------------------------------

cp $grid_dir/C${res}_*mosaic.nc $out_dir

#------------------------------------------------------------------------------------
# Create surface static fields - vegetation type, soil type, etc.
#
# For global grids with a nest, the program is run twice.  First
# to create the fields for the six global tiles.  Then to create
# the fields on the high-res nest.  This is done because the
# ESMF libraries can not interpolate to seven tiles at once.
#------------------------------------------------------------------------------------

export WORK_DIR=$TMPDIR/sfcfields
export SAVE_DIR=$out_dir/fix_sfc
export BASE_DIR=$home_dir
export FIX_FV3=$out_dir
export input_sfc_climo_dir=$home_dir/fix/fix_sfc_climo

if [ $gtype = regional ]; then
  export HALO=$halop1
  ln -fs $out_dir/C${res}_grid.tile${tile}.halo${HALO}.nc $out_dir/C${res}_grid.tile${tile}.nc
  ln -fs $out_dir/C${res}_oro_data.tile${tile}.halo${HALO}.nc $out_dir/C${res}_oro_data.tile${tile}.nc
  export GRIDTYPE=regional
elif [ $gtype = nest ]; then
  export mosaic_file=$out_dir/C${res}_coarse_mosaic.nc
fi

$script_dir/sfc_climo_gen.ksh
err=$?
if [ $err != 0 ]; then
  echo error in sfc_climo_gen
  exit $err
fi

if [ $gtype = regional ]; then
  rm -f $out_dir/C${res}_grid.tile${tile}.nc
  rm -f $out_dir/C${res}_oro_data.tile${tile}.nc
fi

#------------------------------------------------------------------------------------
# Run for the global nest - tile 7.
#------------------------------------------------------------------------------------

if [ $gtype = nest ]; then
  export mosaic_file=$out_dir/C${res}_nested_mosaic.nc
  export GRIDTYPE=nest
  $script_dir/sfc_climo_gen.ksh
  err=$?
  if [ $err != 0 ]; then
    echo error in sfc_climo_gen
    exit $err
  fi
fi

exit
