#!/bin/bash

#BSUB -oo log.grid.%J
#BSUB -eo log.grid.%J
#BSUB -q debug
#BSUB -P GFS-DEV
#BSUB -J grid_fv3
#BSUB -W 0:30
#BSUB -x                 # run not shared
#BSUB -n 24              # total tasks
#BSUB -R span[ptile=24]   # tasks per node
#BSUB -R affinity[core(1):distribute=balance]

#-----------------------------------------------------------------------
# Driver script to create a cubic-sphere based model grid on Dell.
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
#       tasks spread across multiple nodes. The orography code
#       benefits from threads.
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

module purge
module load EnvVars/1.0.2
module load lsf/10.1
module load ips/18.0.1.163
module load impi/18.0.1
module load NetCDF/4.5.0
module load HDF5-serial/1.10.1
module list

#-----------------------------------------------------------------------
# Set grid specs here.
#-----------------------------------------------------------------------

export res=96
export gtype=regional  # 'uniform', 'stretch', 'nest', or 'regional'

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
export TMPDIR=/gpfs/dell1/stmp/$LOGNAME/fv3_grid.$gtype
export out_dir=/gpfs/dell1/stmp/$LOGNAME/C${res}

#-----------------------------------------------------------------------
# Should not need to change anything below here.
#-----------------------------------------------------------------------

export APRUN=time
export APRUN_SFC="mpirun -l"
export OMP_NUM_THREADS=24 # orog code worked best with 24 threads.
export OMP_STACKSIZE=2048m
export machine=WCOSS_DELL_P3

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

$home_dir/ush/fv3gfs_driver_grid.sh

exit
