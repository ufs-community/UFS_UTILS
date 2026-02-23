#!/bin/bash
set -eu

error() {
    set +x
    echo
    echo "$@" 1>&2
    exit 1
}

usage() {
    echo
    echo "Usage: $program [-c] [-m] [-h] [-b]"
    echo
    echo "  -b build the executable"
    echo
    echo "  -c create a new baseline"
    echo
    echo "  -m compare against the new baseline"
    echo
    echo "  -h display this help and exit"
    echo
    echo "  Examples"
    echo
    echo "    './rt.sh -b'  build exe file. compare against the existing baseline"
    echo "    './rt.sh -bc' build exe file. create a new baseline"
    echo "    './rt.sh -m'  do not build exe file. compare against the new baseline"
    echo
}

usage_and_exit() {
    set +x
    usage
    exit $1
}

# Execution starts here.

set -x
test_name="ocnice_prep"
RT_DIR=${RT_DIR:-${PWD}/..}

notlocal=${notlocal:-false}
waitlocal=false
if [[ ${notlocal} == "false" ]]; then
  waitlocal=true
fi

if [[ -f "${RT_DIR}/rt.control" ]]; then
    source "${RT_DIR}/rt.control"
else
    echo "ERROR: Cannot find rt.control script"
    exit 1
fi

program=$(basename $0)
readonly program

# PATHRT - Path to regression tests directory

PATHRT="$(cd "$(dirname $0)" && pwd -P)"
readonly PATHRT
export PATHRT

# PATHTR - Path to the UFS UTILS directory

PATHTR="$(cd "${PATHRT}/../.." && pwd)"
readonly PATHTR

TESTS_FILE="./rt.conf"

source ${HOMEUFSUTILS}/sorc/machine-setup.sh > /dev/null 2>&1
set +x
echo
echo "Machine: ${MACHINE_ID}"
echo
set -x

BASELINE_ROOT=${HOMEreg}/${test_name}/baseline_data
WEIGHTS_ROOT=${HOMEreg}/cpld_gridgen/baseline_data
INPUT_ROOT=${HOMEreg}/${test_name}/input_data
STMP=${WORK_DIR}
ACCOUNT=${PROJECT_CODE}

if [[ ${MACHINE_ID} = wcoss2 ]]; then
    #STMP=${STMP:-/lfs/h2/emc/stmp/$USER}
                  
    #BASELINE_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/ocnice_prep/baseline_data
    
    #WEIGHTS_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/cpld_gridgen/baseline_data
    
    #INPUT_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/ocnice_prep/input_data
    
    #ACCOUNT=${ACCOUNT:-GFS-DEV}
    export APRUN="mpiexec -n 1 -ppn 1 --cpu-bind core"
    #QUEUE=${QUEUE:-dev}
    WLCLK=15
    export NCCMP=nccmp
elif [[ ${MACHINE_ID} = ursa ]]; then
    #STMP=${STMP:-/scratch4/NCEPDEV/stmp/$USER}
    #BASELINE_ROOT=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/ocnice_prep/baseline_data
    #WEIGHTS_ROOT=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    #INPUT_ROOT=/scratch3/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/ocnice_prep/input_data
    #ACCOUNT=${ACCOUNT:-fv3-cpu}
    #QUEUE=${QUEUE:-batch}
    WLCLK=10
    export NCCMP=nccmp
    PARTITION='u1-compute'
elif [[ ${MACHINE_ID} = orion ]]; then
    #STMP=${STMP:-/work/noaa/stmp/$USER}
    #BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/ocnice_prep/baseline_data
    #WEIGHTS_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    #INPUT_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/ocnice_prep/input_data
    #ACCOUNT=${ACCOUNT:-nems}
    #QUEUE=${QUEUE:-batch}
    WLCLK=15
    export NCCMP=nccmp
    PARTITION='orion'
    ulimit -a
elif [[ ${MACHINE_ID} = hercules ]]; then
    #STMP=${STMP:-/work2/noaa/stmp/$USER}
    #BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/ocnice_prep/baseline_data
    #WEIGHTS_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/cpld_gridgen/baseline_data
    #INPUT_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/ocnice_prep/input_data
    #ACCOUNT=${ACCOUNT:-fv3-cpu}
    #QUEUE=${QUEUE:-batch}
    WLCLK=10
    export NCCMP=nccmp
    PARTITION='hercules'
    ulimit -s unlimited
elif [[ ${MACHINE_ID} = jet ]]; then
    #STMP=${STMP:-/lfs5/HFIP/h-nems/$USER}
    #BASELINE_ROOT=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/ocnice_prep/baseline_data
    #WEIGHTS_ROOT=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    #INPUT_ROOT=/lfs5/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/ocnice_prep/input_data
    #ACCOUNT=${ACCOUNT:-h-nems}
    #QUEUE=${QUEUE:-batch}
    WLCLK=10
    export NCCMP=nccmp
    PARTITION="--partition=xjet"
    ulimit -s unlimited
fi

NEW_BASELINE_ROOT=$STMP/ocnice_prep/baseline_data
RUNDIR_ROOT=$STMP/ocnice_prep/rt_$$

BUILD_EXE=false
CREATE_BASELINE=false
while getopts :bcmh opt; do
    case $opt in
	b)
	    BUILD_EXE=true
	    ;;
	c)
	    CREATE_BASELINE=true
	    ;;
	m)
	    BASELINE_ROOT=$NEW_BASELINE_ROOT
	    ;;
	h)
	    usage_and_exit 0
	    ;;
	'?')
            set +x
	    error "$program: invalid option"
	    ;;
    esac
done

compiler=${compiler:-intelllvm}
if [[ "$compiler" == "intelllvm" ]]; then
  if [[ ! -f ${PATHTR}/modulefiles/build.${MACHINE_ID}.$compiler.lua ]];then
    set +x
    echo "IntelLLVM not available. Will use Intel Classic."
    set -x
    compiler=intel
  fi
fi
set +x
echo "Compiler: $compiler"
set -x
export compiler

# Build the executable file
if [[ $BUILD_EXE = true ]]; then
    COMPILE_LOG=compile.log
    cd $PATHTR
    rm -rf $COMPILE_LOG $PATHTR/build $PATHTR/exec $PATHTR/lib
    ./build_all.sh >$PATHRT/$COMPILE_LOG 2>&1 && d=$? || d=$?
    if [[ d -ne 0 ]]; then
        set +x
	error "Build did not finish successfully. Check $COMPILE_LOG"
    else
        set +x
	echo "Build was successful"
        set -x
    fi
    cd $PATHRT
else
    if [[ ! -f $PATHTR/exec/oiprep ]]; then
      set +x
      error "oiprep exe file is not found in $PATHTR/exec/. Try -b to build or -h for help."
    fi
fi

# module use $PATHTR/modulefiles
module use ${HOMEUFSUTILS}/modulefiles
module load build.${MACHINE_ID}.$compiler
if [[ ${MACHINE_ID} = wcoss2 ]]; then
  module load nccmp-D/1.9.0.1
fi
set +x
module list
set -x

export CREATE_BASELINE
if [[ $CREATE_BASELINE = true ]]; then
    rm -rf $NEW_BASELINE_ROOT
    mkdir -p $NEW_BASELINE_ROOT
fi

declare -A tests
all_tests=""

rm -f nccmp_*.log summary.log run_*log RegressionTests_${MACHINE_ID}.$compiler.*.log

# Run tests specified in $TESTS_FILE
i=0
while read -r line || [ "$line" ]; do

    line="${line#"${line%%[![:space:]]*}"}"
    [[ ${#line} == 0 ]] && continue
    [[ $line =~ \# ]] && continue

    LINEVAL=$(echo $line | cut -d'|' -f1 | sed -e 's/^ *//' -e 's/ *$//')
    TEST_SORC=${LINEVAL##mx}
    TEST_FTYP=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    LINEVAL=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
    TEST_DEST=${LINEVAL##mx}
    TEST_NAME=${TEST_SORC}_${TEST_FTYP}_${TEST_DEST}

    RUNDIR=$RUNDIR_ROOT/$TEST_NAME
    BASELINE=$BASELINE_ROOT/$TEST_NAME
    export BASELINE
    NEW_BASELINE=$NEW_BASELINE_ROOT/$TEST_NAME
    export NEW_BASELINE
    mkdir -p $RUNDIR

    export SRCRES=$TEST_SORC
    export DSTRES=$TEST_DEST
    export FTYPE=$TEST_FTYP
    export WEIGHTS=$WEIGHTS_ROOT
    export REGRESSIONTEST_LOG=RegressionTests_${MACHINE_ID}.$compiler.${TEST_NAME}.log

    cp $PATHTR/exec/oiprep $RUNDIR
    cp ./ocnice_prep.sh $RUNDIR
    cp ./parm/ocniceprep.nml.IN $RUNDIR
    cp ./parm/$FTYPE.csv $RUNDIR
    cp $INPUT_ROOT/$FTYPE.nc $RUNDIR
    export RUNDIR
    export TEST_NAME

    if [[ ${MACHINE_ID} = wcoss2 ]]; then

      tests[$i]=$(qsub -V -o run_${TEST_NAME}.log -e run_${TEST_NAME}.log -q $QUEUE  -A $ACCOUNT \
            -l walltime=00:${WLCLK}:00 -N $TEST_NAME -l select=1:ncpus=1:mem=24GB -v RESNAME=$TEST_NAME ./ocnice_prep.sh)

    else

      tests[$i]=$(sbatch --parsable --ntasks-per-node=1 --nodes=1 --mem=24g -t 00:${WLCLK}:00 -A $ACCOUNT -q $QUEUE -J $TEST_NAME \
                -p $PARTITION -o run_${TEST_NAME}.log -e run_${TEST_NAME}.log ./ocnice_prep.sh "$TEST_NAME")

    fi

    all_tests=${all_tests}":"${tests[$i]}

    ((i=i+1))

done <$TESTS_FILE

export target

if [[ ${MACHINE_ID} = wcoss2 ]]; then

  (qsub -V -o /dev/null -e /dev/null -q $QUEUE -A $ACCOUNT -l walltime=00:01:00 \
        -N summary -l select=1:ncpus=1:mem=100MB \
        -W depend=afterok${all_tests} ./rt.summary.sh) &

else

  (sbatch --ntasks=1 --mem=25m -t 0:01:00 -A $ACCOUNT -J summary -o /dev/null -e /dev/null \
       -p $PARTITION --open-mode=append -q $QUEUE -d afterok${all_tests} ./rt.summary.sh) &

fi

if [[ "${waitlocal}" == "true" ]]; then
  sleep_time=0
  echo "Waiting for OCNICE_PREP tests to complete..."
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       mail -s "UFS_UTILS Consistency Test OCNICE_PREP timed out on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
       exit 1
    fi
  done
  mail -s "UFS_UTILS Consistency Test OCNICE_PREP COMPLETED on ${MACHINE_ID}" "${MAILTO}" < "./summary.log"
fi
exit 0
