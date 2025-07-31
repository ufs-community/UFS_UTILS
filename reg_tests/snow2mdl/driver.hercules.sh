#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run snow2mdl consistency tests on Hercules.
#
# Set $DATA_ROOT to your working directory.  Set the project code
# and queue as appropriate.
#
# Invoke the script as follows:  ./$script
#
# Log output is placed in consistency.log.  A summary is
# placed in summary.log
#
# The test fails when its output does not match the baseline file
# as determined by the 'cmp' command.  The baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intelllvm
module load prod_util/2.1.1
module list

ulimit -s unlimited

export DATA_ROOT="${WORK_DIR:-/work2/noaa/stmp/$LOGNAME}"
export DATA_ROOT="${DATA_ROOT}/reg-tests/snow2mdl"

PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
QUEUE="${QUEUE:-batch}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

rm -fr $DATA_ROOT

export OMP_NUM_THREADS=1

export HOMEreg=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/snow2mdl
export HOMEgfs=$PWD/../..
export WGRIB=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/grib-util-1.3.0-wenl3in/bin/wgrib
export WGRIB2=/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env/install/intel/2021.9.0/wgrib2-3.1.1-v7xhwos/bin/wgrib2

# The first test uses the hemispheric air force/afwa data, which was used in OPS.

export DATA="${DATA_ROOT}/test.hemi"
TEST1=$(sbatch --parsable -J snow.hemi -A $PROJECT_CODE -o consistency.log \
        -e consistency.log --ntasks=1 -q $QUEUE -t 00:03:00 ./snow2mdl.hemi.sh)

# This tests the afwa global grib2 data, which is used in OPS.

export DATA="${DATA_ROOT}/test.global"
TEST2=$(sbatch --parsable -J snow.global -A $PROJECT_CODE -o consistency.log \
        -e consistency.log --ntasks=1 -q $QUEUE -t 00:03:00 --open-mode=append \
        -d afterok:$TEST1 ./snow2mdl.global.sh)

# Create the summary file.

sbatch --ntasks=1 -t 0:01:00 -A $PROJECT_CODE -J snow_summary -o consistency.log -e consistency.log \
       --open-mode=append -q $QUEUE -d afterok:$TEST2 << EOF
#!/bin/bash
grep -a '<<<' consistency.log  > summary.log
EOF

exit 0
