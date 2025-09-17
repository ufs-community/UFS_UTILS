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
if [[ "$target" == "wcoss2" ]];then
  module load nccmp-D/1.9.0.1
fi
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
  PARTITION="--partition=xjet"
elif [[ "$target" == "ursa" ]];then
  WORK_DIR="${WORK_DIR:-/scratch4/NCEPDEV/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/regrid_sfc/
  export APRUN_REGRID=srun
  PARTITION=''
elif [[ "$target" == "orion" ]];then
  WORK_DIR="${WORK_DIR:-/work/noaa/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/work/noaa/nems/role-nems/ufs_utils/reg_tests/regrid_sfc
  export APRUN_REGRID=srun
  PARTITION=''
  ulimit -a
elif [[ "$target" == "hercules" ]];then
  WORK_DIR="${WORK_DIR:-/work2/noaa/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-fv3-cpu}"
  QUEUE="${QUEUE:-batch}"
  export HOMEreg=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/regrid_sfc
  export APRUN_REGRID=srun
  PARTITION=''
elif [[ "$target" == "wcoss2" ]];then
  WORK_DIR="${WORK_DIR:-/lfs/h2/emc/stmp/$LOGNAME}"
  PROJECT_CODE="${PROJECT_CODE:-GFS-DEV}"
  QUEUE="${QUEUE:-dev}"
  export HOMEreg=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/regrid_sfc
  export APRUN_REGRID="mpiexec -n 6 -ppn 6 --cpu-bind core"
elif [[ "$target" == "noaacloud" ]];then
  WORK_DIR="${WORK_DIR:-/contrib/$LOGNAME}/dev/UFS_UTILS"
  PROJECT_CODE="${PROJECT_CODE:-${USER}}"
  QUEUE="${QUEUE:-process}"
  export HOMEreg=/contrib/ufs_utils/reg_tests/regrid_sfc/
  export APRUN_REGRID="srun --mpi=pmi2 -l -n 6"
  PARTITION='--partition process'
fi

DATA_DIR="${WORK_DIR}/reg-tests/regrid_sfc"
export NWPROD=$PWD/../..

LOG_FILE=consistency.log01
rm -f $LOG_FILE
export DATA="${DATA_DIR}/test1"
if [[ "$target" == "wcoss2" ]];then
  TEST1=$(qsub -V -o $LOG_FILE -e $LOG_FILE -q $QUEUE -A $PROJECT_CODE -l walltime=00:05:00 \
        -N gauss2fv3incr -l select=1:ncpus=6:ompthreads=1:mem=10GB ./gauss2fv3incr.sh)
else
  echo "PARTITION: $PARTITION"
  TEST1=$(sbatch --parsable --ntasks-per-node=6 --nodes=1 -t 0:05:00 -A $PROJECT_CODE -q $QUEUE -J gauss2fv3incr \
      $PARTITION -o $LOG_FILE -e $LOG_FILE ./gauss2fv3incr.sh)
fi

LOG_FILE=consistency.log
rm -f $LOG_FILE summary.log

if [[ "$target" == "wcoss2" ]];then

this_dir=$PWD
qsub -V -o ${LOG_FILE} -e ${LOG_FILE} -q $QUEUE -A $PROJECT_CODE -l walltime=00:01:00 \
        -N summary -l select=1:ncpus=1:mem=100MB -W depend=afterok:$TEST1 << EOF
#!/bin/bash
cd $this_dir
grep -a '<<<' ${LOG_FILE}?? | grep -v echo > ./summary.log
EOF

else

  echo "PARTITION: $PARTITION"
sbatch --nodes=1  -t 0:01:00 -A $PROJECT_CODE -J summary -o $LOG_FILE -e $LOG_FILE \
       $PARTITION --open-mode=append -q $QUEUE -d afterok:$TEST1 << EOF
#!/bin/bash
grep -a '<<<' ${LOG_FILE}* > ./summary.log
EOF

fi

exit
