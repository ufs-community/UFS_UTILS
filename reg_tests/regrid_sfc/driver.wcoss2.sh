#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run regrid_sfc consistency test on WCOSS2.
#
# Set $WORK_DIR to your working directory.  Set the project code 
# and queue as appropriate.
#
# Invoke the script from the command line as follows:  ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility.  This baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

compiler=${compiler:-"intel"}

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.$compiler
module load netcdf
module load nccmp
module list

WORK_DIR="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"

PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
QUEUE="${QUEUE:-dev}"

#-----------------------------------------------------------------------------
# Should not have to change anything below.
#-----------------------------------------------------------------------------

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [ "$UPDATE_BASELINE" = "TRUE" ]; then
  source ../get_hash.sh
fi

DATA_DIR="${WORK_DIR}/reg-tests/regrid_sfc"

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/regrid_sfc

export OMP_NUM_THREADS_CY=2
export OMP_PLACES=cores

export NWPROD=$PWD/../..

reg_dir=$PWD

LOG_FILE=consistency.log
rm -f ${LOG_FILE}*

export DATA="${DATA_DIR}/test1"
TEST1=$(qsub -V -o ${LOG_FILE}01 -e ${LOG_FILE}01 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N gauss2fv3incr -l select=1:ncpus=6 $PWD/gauss2fv3incr.sh)

qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N cycle_summary -l select=1:ncpus=1:mem=100MB -W depend=afterok:$TEST1 << EOF
#!/bin/bash
cd $reg_dir
grep -a '<<<' ${LOG_FILE}?? | grep -v echo > summary.log
EOF

exit
