#!/bin/bash

#BSUB -L /bin/sh
#BSUB -P GFS-DEV
#BSUB -oo log.grid.%J
#BSUB -eo log.grid.%J
#BSUB -J grid_fv3
#BSUB -q debug
#BSUB -M 2400
#BSUB -W 00:30
#BSUB -extsched 'CRAYLINUX[]'

#-----------------------------------------------------------------------
# Driver script to create a cubic-sphere based model grid on Cray.
#
# Produces the following files (netcdf, each tile in separate file):
#   1) 'mosaic' and 'grid' files containing lat/lon and other
#      records that describe the model grid.
#   2) 'oro' files containing land mask, terrain and gravity
#      wave drag fields.
#   3) Surface climo fields, such as soil type, vegetation
#      greenness and albedo.
#
# Note: The sfc_climo_gen program only runs with an
#       mpi task count that is a multiple of six.  This is
#       an ESMF library requirement.  Large grids may require
#       tasks spread across multiple nodes.
#
# To run, do the following:
#
#   1) Set "C" resolution, "res" - Example: res=96.
#   2) Set grid type ("gtype").  Valid choices are
#         "uniform"  - global uniform grid
#         "stretch"  - global stretched grid
#         "nest"     - global stretched grid with nest
#         "regional_gfdl" - stand-alone gfdl regional grid
#         "regional_esg"  - stand-alone extended Schmidt gnomonic
#                           (esg) regional grid
#   3) For "uniform" and "regional_gfdl" grids - to include lake
#      fraction and depth, set "add_lake" to true, and the
#      "lake_cutoff" value.
#   4) For "stretch" and "nest" grids, set the stretching factor -
#       "stretch_fac", and center lat/lon of highest resolution
#      tile - "target_lat" and "target_lon".
#   5) For "nest" grids, set the refinement ratio - "refine_ratio",
#      the starting/ending i/j index location within the parent
#      tile - "istart_nest", "jstart_nest", "iend_nest", "jend_nest"
#   6) For "regional_gfdl" grids, set the "halo".  Default is three
#      rows/columns.
#   7) For "regional_esg" grids, set center lat/lon of grid,
#      - "target_lat/lon" - the i/j dimensions - "i/jdim", the
#      x/y grid spacing - "delx/y", and halo.
#   8) Set working directory - TEMP_DIR - and path to the repository
#      clone - home_dir.
#   9) Check settings for 'make_gsl_orog' and 'veg_type_src'
#      below.
#  10) Submit script: "cat $script | bsub".
#  11) All files will be placed in "out_dir".
#
#-----------------------------------------------------------------------

source ../sorc/machine-setup.sh > /dev/null 2>&1
module use ../modulefiles
module load build.$target.intel
module list

#-----------------------------------------------------------------------
# Set grid specs here.
#-----------------------------------------------------------------------

export gtype=uniform           # 'uniform', 'stretch', 'nest', 
                               # 'regional_gfdl', 'regional_esg'
export make_gsl_orog=false     # 'true' if user needs 'oro' files for GSL
                               # orographic drag suite
export veg_type_src="modis.igbp.0.05"  #  veg type data.
                                # For viirs-based vegetation type data, set to:
                                # 1) "viirs.igbp.0.05" for global 5km data
                                # 2) "viirs.igbp.0.1" for global 10km data
                                # 3) "viirs.igbp.0.03" for global 3km data
                                # 4) "viirs.igbp.conus.0.01" for regional 1km data
                                # For the modis-based data, set to:
                                # 1) "modis.igbp.0.05" for global 5km data
                                # 2) "modis.igbp.0.03" for global 3km data
                                # 3) "modis.igbp.conus.0.01" for regional 1km data

if [ $gtype = uniform ]; then
  export res=96
  export add_lake=false        # Add lake frac and depth to orography data.
  export lake_cutoff=0.20      # lake frac < lake_cutoff ignored when add_lake=T
elif [ $gtype = stretch ]; then
  export res=96
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
elif [ $gtype = nest ] || [ $gtype = regional_gfdl ]; then
  export add_lake=false        # Add lake frac and depth to orography data.
  export lake_cutoff=0.20      # lake frac < lake_cutoff ignored when add_lake=T
  export res=768
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=38.5       # Center latitude of the highest resolution tile
  export refine_ratio=3        # The refinement ratio
  export istart_nest=123       # Starting i-direction index of nest grid in parent tile supergrid
  export jstart_nest=331       # Starting j-direction index of nest grid in parent tile supergrid
  export iend_nest=1402        # Ending i-direction index of nest grid in parent tile supergrid
  export jend_nest=1194        # Ending j-direction index of nest grid in parent tile supergrid
  export halo=3                # Lateral boundary halo
elif [ $gtype = regional_esg ] ; then
  export res=-999              # equivalent resolution is computed
  export target_lon=-97.5      # Center longitude of grid
  export target_lat=35.5       # Center latitude of grid
  export idim=301              # Dimension of grid in 'i' direction
  export jdim=200              # Dimension of grid in 'j' direction
  export delx=0.0585           # Grid spacing (in degrees) in the 'i' direction
                               # on the SUPERGRID (which has twice the resolution of
                               # the model grid).  The physical grid spacing in the 'i'
                               # direction is related to delx as follows:
                               #    distance = 2*delx*(circumf_Earth/360 deg)
  export dely=0.0585           # Grid spacing (in degrees) in the 'j' direction.
  export halo=3                # number of row/cols for halo
fi

#-----------------------------------------------------------------------
# Check paths.
#   home_dir - location of repository.
#   TEMP_DIR - working directory.
#   out_dir  - where files will be placed upon completion.
#-----------------------------------------------------------------------

export home_dir=$LS_SUBCWD/..
export TEMP_DIR=/gpfs/hps3/stmp/$LOGNAME/fv3_grid.$gtype
export out_dir=/gpfs/hps3/stmp/$LOGNAME/my_grids

#-----------------------------------------------------------------------
# Should not need to change anything below here.
#-----------------------------------------------------------------------

export NODES=1
export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
export APRUN_SFC="aprun -j 1 -n 24 -N 24"
# The orography code is optimized for six threads.
export OMP_NUM_THREADS=6
export OMP_STACKSIZE=2048m
export KMP_AFFINITY=disabled
export machine=WCOSS_C
export NCDUMP=/gpfs/hps/usrx/local/prod/NetCDF/4.2/intel/sandybridge/bin/ncdump

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

$home_dir/ush/fv3gfs_driver_grid.sh

exit
