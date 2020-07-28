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
#   3) For "stretch" and "nest" grids, set the stretching factor -
#       "stretch_fac", and center lat/lon of highest resolution
#      tile - "target_lat" and "target_lon".
#   4) For "nest" grids, set the refinement ratio - "refine_ratio",
#      the starting/ending i/j index location within the parent
#      tile - "istart_nest", "jstart_nest", "iend_nest", "jend_nest"
#   5) For "regional_gfdl" grids, set the "halo".  Default is three
#      rows/columns.
#   6) For "regional_esg" grids, set center lat/lon of grid,
#      - "target_lat/lon" - the i/j dimensions - "i/jdim", the
#      x/y grid spacing - "delx/y", and halo.
#   7) Set working directory - TMPDIR - and path to the repository
#      clone - home_dir.
#   8) Submit script: "cat $script | bsub".
#   9) All files will be placed in "out_dir".
#
#-----------------------------------------------------------------------

source ../sorc/machine-setup.sh > /dev/null 2>&1
source ../modulefiles/build.$target
module list

#-----------------------------------------------------------------------
# Set grid specs here.
#-----------------------------------------------------------------------

export gtype=uniform        # 'uniform', 'stretch', 'nest'
                           # 'regional_gfdl', 'regional_esg'

if [ $gtype = uniform ]; then
  export res=96
elif [ $gtype = stretch ]; then
  export res=96
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
elif [ $gtype = nest ] || [ $gtype = regional_gfdl ]; then
  export res=96
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
  export refine_ratio=3        # The refinement ratio
  export istart_nest=27        # Starting i-direction index of nest grid in parent tile supergrid
  export jstart_nest=37        # Starting j-direction index of nest grid in parent tile supergrid
  export iend_nest=166         # Ending i-direction index of nest grid in parent tile supergrid
  export jend_nest=164         # Ending j-direction index of nest grid in parent tile supergrid
  export halo=3
elif [ $gtype = regional_esg ] ; then
  export res=-999              # equivalent res is computed.
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
#   TMPDIR   - working directory.
#   out_dir  - where files will be placed upon completion.
#-----------------------------------------------------------------------

export home_dir=$LS_SUBCWD/..
export TMPDIR=/gpfs/hps3/stmp/$LOGNAME/fv3_grid.$gtype
export out_dir=/gpfs/hps3/stmp/$LOGNAME/my_grids

export NODES=1
export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
export APRUN_SFC="aprun -j 1 -n 24 -N 24"
# The orography code is optimized for six threads.
export OMP_NUM_THREADS=6
export OMP_STACKSIZE=2048m
export KMP_AFFINITY=disabled
export machine=WCOSS_C

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

$home_dir/ush/fv3gfs_driver_grid.sh

exit
