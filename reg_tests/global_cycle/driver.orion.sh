#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle consistency tests on Orion.
#
# Set $WORK_DIR to your working directory.  Set the project code and
# and queue as appropriate.
#
# Invoke the script from command line as follows:  ./$script
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

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

export WORK_DIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"

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

export DATA_DIR="${WORK_DIR}/reg-tests/global-cycle"

export HOMEreg=/work/noaa/nems/role-nems/ufs_utils/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2

export APRUNCY="srun"

export NWPROD=$PWD/../..

reg_dir=$PWD

LOG_FILE=consistency.log01
export DATA="${DATA_DIR}/test1"
export COMOUT=$DATA
TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c768.fv3gfs \
      -o $LOG_FILE -e $LOG_FILE ./C768.fv3gfs.sh)

LOG_FILE=consistency.log02
export DATA="${DATA_DIR}/test2"
export COMOUT=$DATA
TEST2=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c192.lndincsoilnoahmp \
      -o $LOG_FILE -e $LOG_FILE ./C192.lndincsoilnoahmp.sh)

LOG_FILE=consistency.log03
export DATA="${DATA_DIR}/test3"
export COMOUT=$DATA
TEST3=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J c768.lndincsnow \
      -o $LOG_FILE -e $LOG_FILE ./C768.lndincsnow.sh)

LOG_FILE=consistency.log
sbatch --nodes=1 -t 0:01:00 -A $PROJECT_CODE -J chgres_summary -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q $QUEUE -d\
      afterok:$TEST1:$TEST2:$TEST3 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}*  > summary.log
EOF

exit
