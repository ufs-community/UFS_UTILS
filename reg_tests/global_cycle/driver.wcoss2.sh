#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run global_cycle consistency test on WCOSS2.
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

DATA_DIR="${WORK_DIR}/reg-tests/global-cycle"

export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/global_cycle

export OMP_NUM_THREADS_CY=2
export OMP_PLACES=cores

export APRUNCY="mpiexec -n 6 -ppn 6 --cpu-bind core --depth ${OMP_NUM_THREADS_CY}"

export NWPROD=$PWD/../..

reg_dir=$PWD

export NCCMP=/lfs/h2/emc/global/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp

LOG_FILE=consistency.log
rm -f ${LOG_FILE}*

export DATA="${DATA_DIR}/test1"
export COMOUT=$DATA
TEST1=$(qsub -V -o ${LOG_FILE}01 -e ${LOG_FILE}01 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c768.fv3gfs -l select=1:ncpus=12:mem=12GB $PWD/C768.fv3gfs.sh)

export DATA="${DATA_DIR}/test2"
export COMOUT=$DATA
TEST2=$(qsub -V -o ${LOG_FILE}02 -e ${LOG_FILE}02 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c192.gsi_lndincsoilnoahmp -l select=1:ncpus=12:mem=8GB $PWD/C192.gsi_lndincsoilnoahmp.sh)

export DATA="${DATA_DIR}/test3"
export COMOUT=$DATA
TEST3=$(qsub -V -o ${LOG_FILE}03 -e ${LOG_FILE}03 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c768.lndincsnow -l select=1:ncpus=12:mem=8GB $PWD/C768.lndincsnow.sh)

export DATA="${DATA_DIR}/test4"
export COMOUT=$DATA
TEST4=$(qsub -V -o ${LOG_FILE}04 -e ${LOG_FILE}04 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c48.noahmp.frac -l select=1:ncpus=12:mem=8GB $PWD/C48.noahmp.fracgrid.sh)

export DATA="${DATA_DIR}/test5"
export COMOUT=$DATA
TEST5=$(qsub -V -o ${LOG_FILE}05 -e ${LOG_FILE}05 -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N c192.jedi_lndincsoilnoahmp -l select=1:ncpus=12:mem=8GB $PWD/C192.jedi_lndincsoilnoahmp.sh)

qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N cycle_summary -l select=1:ncpus=1:mem=100MB -W depend=afterok:$TEST1:$TEST2:$TEST3:$TEST4:$TEST5 << EOF
#!/bin/bash
cd $reg_dir
grep -a '<<<' ${LOG_FILE}?? | grep -v echo > summary.log
EOF

exit
