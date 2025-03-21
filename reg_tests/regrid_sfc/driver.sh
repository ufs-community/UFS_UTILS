#!/bin/bash

#-----------------------------------------------------------------------------
#
# Run regrid_sfc consistency tests.
#
# Set $WORK_DIR to your working directory. 
# Set the $PROJECT_CODE and $QUEUE as appropriate.
#
# Invoke the script from command line as follows:  ./$script
#
# Log output is placed in consistency.log??.  A summary is
# placed in summary.log
#
# A test fails when its output does not match the baseline files
# as determined by the 'nccmp' utility. The baseline files are
# stored in HOMEreg.
#
#-----------------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
compiler=${compiler:-intelllvm}
if [[ "$compiler" == "intelllvm" ]]; then
  if [[ ! -f ../../modulefiles/build.$target.$compiler.lua ]];then
     set +x
     echo "IntelLLVM not available. Will use Intel Classic."
     set -x
    compiler=intel
  fi
fi
module load build.$target.$compiler
set +x
module list
set -x

export UPDATE_BASELINE="FALSE"
#export UPDATE_BASELINE="TRUE"

if [[ "$UPDATE_BASELINE" == "TRUE" ]]; then
  source ../get_hash.sh
fi

if [[ "$target" == "jet" ]];then
  export WORK_DIR="${WORK_DIR:-/lfs5/HFIP/emcda/$LOGNAME/stmp}"
  PROJECT_CODE="${PROJECT_CODE:-hfv3gfs}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/regrid_sfc
  export APRUN_REGRID=srun
  PARTITION=xjet
elif [[ "$target" == "hera" ]];then
  WORK_DIR="${WORK_DIR:-/scratch2/NCEPDEV/stmp1/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/regrid_sfc/
  export APRUN_REGRID=srun
  PARTITION=hera
elif [[ "$target" == "orion" ]];then
  WORK_DIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/work/noaa/nems/role-nems/ufs_utils/reg_tests/regrid_sfc
  export APRUN_REGRID=srun
  PARTITION=orion
elif [[ "$target" == "hercules" ]];then
  WORK_DIR="${WORK_DIR:-/work2/noaa/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/regrid_sfc
  export APRUN_REGRID=srun
  PARTITION=hercules
fi

DATA_DIR="${WORK_DIR}/reg-tests/regrid_sfc"
export NWPROD=$PWD/../..

LOG_FILE=consistency.log01
rm -f $LOG_FILE
export DATA="${DATA_DIR}/test1"
TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J gauss2fv3incr \
      --partition=$PARTITION -o $LOG_FILE -e $LOG_FILE ./gauss2fv3incr.sh)

LOG_FILE=consistency.log
rm -f $LOG_FILE summary.log
sbatch --partition=$PARTITION --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J summary -o $LOG_FILE -e $LOG_FILE \
       --open-mode=append -q $QUEUE -d\
       afterok:$TEST1 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}* > ./summary.log
EOF

exit
