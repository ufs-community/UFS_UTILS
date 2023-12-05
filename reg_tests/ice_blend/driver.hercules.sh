#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run ice_blend consistency test on Hercules.
#
# Set $DATA to your working directory.  Set the project code (SBATCH -A)
# and queue (SBATCH -q) as appropriate.
#
# Invoke the script as follows:  sbatch $script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline file is
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

#SBATCH -J ice_blend
#SBATCH -A fv3-cpu
#SBATCH --open-mode=truncate
#SBATCH -o consistency.log
#SBATCH -e consistency.log
#SBATCH --ntasks=1
#SBATCH -q debug
#SBATCH -t 00:03:00

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

ulimit -s unlimited

export DATA="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
export DATA="${DATA}/reg-tests/ice-blend"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

export WGRIB=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/grib-util-1.3.0-wenl3in/bin/wgrib
export WGRIB2=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/wgrib2-3.1.1-v7xhwos/bin/wgrib2
export COPYGB=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/grib-util-1.3.0-wenl3in/bin/copygb
export COPYGB2=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/grib-util-1.3.0-wenl3in/bin/copygb2
export CNVGRIB=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/grib-util-1.3.0-wenl3in/bin/cnvgrib

#TODO Update to official location when testing is complete
#export HOMEreg=/work/noaa/nems/role-nems/ufs_utils/hercules/reg_tests/ice_blend
export HOMEreg=/work/noaa/global/dhuber/noscrub/ufs_utils/reg_tests/ice_blend
export HOMEgfs=$PWD/../..

rm -fr $DATA

./ice_blend.sh

exit 0
