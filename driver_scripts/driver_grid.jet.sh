#!/bin/bash

#SBATCH -J fv3_grid_driver
#SBATCH -A hfv3gfs
#SBATCH --open-mode=truncate
#SBATCH -o log.fv3_grid_driver
#SBATCH -e log.fv3_grid_driver
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --partition=xjet
#SBATCH -q batch
#SBATCH -t 00:10:00

#-----------------------------------------------------------------------
# Driver script to create a cubic-sphere based model grid on Jet.
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
#       tasks spread across multiple nodes.  The orography
#       code benefits from threads.
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
#   3) For "uniform" grids - to include lake fraction and
#      depth, set "add_lake" to true, and the "lake_cutoff" value.
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
#   8) Set working directory - TMPDIR - and path to the repository
#      clone - home_dir.
#   9) Submit script: "sbatch $script".
#  10) All files will be placed in "out_dir".
#
#-----------------------------------------------------------------------

set -x

source ../sorc/machine-setup.sh > /dev/null 2>&1
source ../modulefiles/build.$target
module list

#-----------------------------------------------------------------------
# Set grid specs here.
#-----------------------------------------------------------------------

export gtype=uniform       # 'uniform', 'stretch', 'nest'
                           # 'regional_gfdl', 'regional_esg'

if [ $gtype = uniform ]; then
  export res=96
  export add_lake=false   # Add lake frac and depth to orography data.
                          # Uniform grids only.
  export lake_cutoff=0.20 # lake frac less than lake_cutoff is ignored
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
#   TMPDIR   - working directory.
#   out_dir  - where files will be placed upon completion. 
#-----------------------------------------------------------------------

export home_dir=$SLURM_SUBMIT_DIR/..
export TMPDIR=/lfs4/HFIP/emcda/$LOGNAME/stmp/fv3_grid.$gtype
export out_dir=/lfs4/HFIP/emcda/$LOGNAME/stmp/my_grids

#-----------------------------------------------------------------------
# Should not need to change anything below here.
#-----------------------------------------------------------------------

export APRUN=time
export APRUN_SFC=srun
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=2048m
export machine=JET

ulimit -a
ulimit -s unlimited

#-----------------------------------------------------------------------
# Start script.
#-----------------------------------------------------------------------

$home_dir/ush/fv3gfs_driver_grid.sh

exit
