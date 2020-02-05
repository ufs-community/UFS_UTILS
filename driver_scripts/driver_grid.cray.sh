#!/bin/bash

#BSUB -L /bin/sh
#BSUB -P FV3GFS-T2O
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
#         "regional" - stand-alone regional grid
#   3) For "stretch" and "nest" grids, set the stretching factor -
#       "stretch_fac", and center lat/lon of highest resolution
#      tile - "target_lat" and "target_lon".
#   4) For "nest" grids, set the refinement ratio - "refine_ratio",
#      the starting/ending i/j index location within the parent
#      tile - "istart_nest", "jstart_nest", "iend_nest", "jend_nest"
#   5) For "regional" grids, set the "halo".  Default is three
#      rows/columns.
#   6) Set working directory - TMPDIR - and path to the repository
#      clone - home_dir.
#   7) Submit script: "cat $script | bsub".
#   8) All files will be placed in "out_dir".
#
#-----------------------------------------------------------------------

. $MODULESHOME/init/sh
module load PrgEnv-intel cfp-intel-sandybridge/1.1.0
module list

#-----------------------------------------------------------------------
# Set grid specs here.
#-----------------------------------------------------------------------

export res=96
export gtype=regional   # 'uniform', 'stretch', 'nest', or 'regional'

if [ $gtype = stretch ]; then
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
elif [ $gtype = nest ] || [ $gtype = regional ]; then
  export stretch_fac=1.5       # Stretching factor for the grid
  export target_lon=-97.5      # Center longitude of the highest resolution tile
  export target_lat=35.5       # Center latitude of the highest resolution tile
  export refine_ratio=3        # The refinement ratio
  export istart_nest=27        # Starting i-direction index of nest grid in parent tile supergrid
  export jstart_nest=37        # Starting j-direction index of nest grid in parent tile supergrid
  export iend_nest=166         # Ending i-direction index of nest grid in parent tile supergrid
  export jend_nest=164         # Ending j-direction index of nest grid in parent tile supergrid
  export halo=3
fi

#-----------------------------------------------------------------------
# Check paths.
#   home_dir - location of repository.
#   TMPDIR   - working directory.
#   out_dir  - where files will be placed upon completion.
#-----------------------------------------------------------------------

export home_dir=$LS_SUBCWD/..
export TMPDIR=/gpfs/hps3/stmp/$LOGNAME/fv3_grid.$gtype
export out_dir=/gpfs/hps3/stmp/$LOGNAME/C${res}

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
