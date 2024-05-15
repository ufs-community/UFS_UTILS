#!/bin/bash
set -eu

SECONDS=0

error() {
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
    usage
    exit $1
}

check_results() {

    [ -o xtrace ] && set_x='set -x' || set_x='set +x'
    set +x

    local test_status=PASS
    # verification run
    if [[ $CREATE_BASELINE = false ]]; then

	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo "Baseline dir = $BASELINE" | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo "Checking test $TEST_NAME results ...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

	for file in $BASELINE/*.nc; do
	    printf %s "Comparing " $(basename ${file}) "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

	    if [[ ! -f $RUNDIR/$(basename ${file}) ]]; then
		echo "....MISSING file" | tee -a $PATHRT/$REGRESSIONTEST_LOG
		test_status=FAIL
	    else
		$NCCMP -dmfqS -w format $(basename ${file}) $file >>${PATHRT}/nccmp_${TEST_NAME}.log 2>&1 && d=$? || d=$?
		if [[ $d -ne 0 ]]; then
		    echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
		    test_status=FAIL
		else
		    echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
		fi
	    fi
	done
	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

	# baseline creation run
    else

	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo "Moving baseline files to $NEW_BASELINE ...." | tee -a $PATHRT/$REGRESSIONTEST_LOG
	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

	mkdir -p $NEW_BASELINE

	for file in *.nc; do
	    printf %s "Moving " $file "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

	    cp $file $NEW_BASELINE/$file && d=$? || d=$?
	    if [[ $d -ne 0 ]]; then
		echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
		test_status=FAIL
	    else
		echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
	    fi
	done
	echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

    fi

    if [[ $test_status == FAIL ]]; then
	echo "$TEST_NAME failed" >> $PATHRT/fail_test_$TEST_NAME
    fi
}

readonly program=$(basename $0)
# PATHRT - Path to regression tests directory
readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
export PATHRT
# PATHTR - Path to the UFS UTILS directory
readonly PATHTR="$(cd $PATHRT/../.. && pwd)"
export PATHTR
TESTS_FILE="$PATHRT/rt.conf"
export TEST_NAME=

# for C3072 on hera, use WLCLK=60 and MEM="--exclusive"
WLCLK_dflt=35
export WLCLK=$WLCLK_dflt
MEM_dflt="--mem=16g"
export MEM=$MEM_dflt

cd $PATHRT
export compiler=${compiler:-intel}
source $PATHTR/sorc/machine-setup.sh >/dev/null 2>&1
echo "Machine: $target"
echo "Compiler: $compiler"

COMPILE_LOG=compile.log
REGRESSIONTEST_LOG=RegressionTests_$target.$compiler.log
rm -f fail_test* $COMPILE_LOG run_*.log nccmp_*.log summary.log

if [[ $target = wcoss2 ]]; then
    STMP=${STMP:-/lfs/h2/emc/stmp/$USER}
    export MOM6_FIXDIR=/lfs/h2/emc/global/noscrub/emc.global/FIX/fix/mom6/20220805
    BASELINE_ROOT=/lfs/h2/emc/nems/noscrub/emc.nems/UFS_UTILS/reg_tests/cpld_gridgen/baseline_data
    ACCOUNT=${ACCOUNT:-GFS-DEV}
    export APRUN="mpiexec -n 1 -ppn 1 --cpu-bind core"
    QUEUE=${QUEUE:-dev}
    SBATCH_COMMAND="./cpld_gridgen.sh"
    NCCMP=/lfs/h2/emc/global/noscrub/George.Gayno/util/nccmp/nccmp-1.8.5.0/src/nccmp
elif [[ $target = hera ]]; then
    STMP=${STMP:-/scratch1/NCEPDEV/stmp4/$USER}
    export MOM6_FIXDIR=/scratch1/NCEPDEV/global/glopara/fix/mom6/20220805
    BASELINE_ROOT=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    ACCOUNT=${ACCOUNT:-nems}
    QUEUE=${QUEUE:-batch}
    NCCMP=nccmp
    PARTITION=hera
    SBATCH_COMMAND="./cpld_gridgen.sh"
elif [[ $target = orion ]]; then
    STMP=${STMP:-/work/noaa/stmp/$USER}
    export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/20220805
    BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    ACCOUNT=${ACCOUNT:-nems}
    QUEUE=${QUEUE:-batch}
    NCCMP=nccmp
    PARTITION=orion
    ulimit -s unlimited
    SBATCH_COMMAND="./cpld_gridgen.sh"
elif [[ $target = hercules ]]; then
    STMP=${STMP:-/work/noaa/stmp/$USER}
    export MOM6_FIXDIR=/work/noaa/global/glopara/fix/mom6/20220805
    BASELINE_ROOT=/work/noaa/nems/role-nems/ufs_utils.hercules/reg_tests/cpld_gridgen/baseline_data
    ACCOUNT=${ACCOUNT:-nems}
    QUEUE=${QUEUE:-batch}
    NCCMP=nccmp
    PARTITION=hercules
    ulimit -s unlimited
    SBATCH_COMMAND="./cpld_gridgen.sh"
elif [[ $target = jet ]]; then
    STMP=${STMP:-/lfs4/HFIP/h-nems/$USER}
    export MOM6_FIXDIR=/lfs4/HFIP/hfv3gfs/glopara/git/fv3gfs/fix/mom6/20220805
    BASELINE_ROOT=/lfs4/HFIP/hfv3gfs/emc.nemspara/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
    ACCOUNT=${ACCOUNT:-h-nems}
    QUEUE=${QUEUE:-batch}
    NCCMP=nccmp
    PARTITION=xjet
    ulimit -s unlimited
    SBATCH_COMMAND="./cpld_gridgen.sh"
fi
NEW_BASELINE_ROOT=$STMP/CPLD_GRIDGEN/BASELINE
RUNDIR_ROOT=$STMP/CPLD_GRIDGEN/rt_$$

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
	    error "$program: invalid option"
	    ;;
    esac
done

# Build the executable file
if [[ $BUILD_EXE = true ]]; then
    cd $PATHTR
    rm -rf $PATHTR/build $PATHTR/exec $PATHTR/lib
    ./build_all.sh >$PATHRT/$COMPILE_LOG 2>&1 && d=$? || d=$?
    if [[ d -ne 0 ]]; then
	error "Build did not finish successfully. Check $COMPILE_LOG"
    else
	echo "Build was successful"
    fi
fi

if [[ ! -f $PATHTR/exec/cpld_gridgen ]]; then
    error "cpld_gridgen exe file is not found in $PATHTR/exe/. Try -b to build or -h for help."
else
    echo "cpld_gridgen exe file is found in $PATHTR/exec/"
fi

module use $PATHTR/modulefiles
module load build.$target.$compiler
module list

if [[ $CREATE_BASELINE = true ]]; then
    rm -rf $NEW_BASELINE_ROOT
    mkdir -p $NEW_BASELINE_ROOT
fi

date > $PATHRT/$REGRESSIONTEST_LOG
echo "Start Regression test" | tee -a  $PATHRT/$REGRESSIONTEST_LOG

# Run tests specified in $TESTS_FILE
while read -r line || [ "$line" ]; do

    line="${line#"${line%%[![:space:]]*}"}"
    [[ ${#line} == 0 ]] && continue
    [[ $line =~ \# ]] && continue

    TEST_NAME=$(echo $line | cut -d'|' -f1 | sed -e 's/^ *//' -e 's/ *$//')
    MOSAICRES=${TEST_NAME%_*}
    TEST_NAME=${TEST_NAME##*_}

    cd $PATHRT
    RUNDIR=$RUNDIR_ROOT/$TEST_NAME
    BASELINE=$BASELINE_ROOT/$TEST_NAME
    NEW_BASELINE=$NEW_BASELINE_ROOT/$TEST_NAME
    mkdir -p $RUNDIR

    # OUTDIR_PATH is passed down to $PATHTR/ush/cpld_gridgen.sh
    # It MUST be set
    export OUTDIR_PATH=$RUNDIR
    export MOSAICRES=$MOSAICRES

    cp $PATHTR/exec/cpld_gridgen $RUNDIR
    cp $PATHTR/ush/cpld_gridgen.sh $RUNDIR
    cp $PATHRT/parm/grid.nml.IN $RUNDIR
    cd $RUNDIR

    if [[ $target = wcoss2 ]]; then

	#   rm -f $RUNDIR/bad.${TEST_NAME}

	TEST=$(qsub -V -o $PATHRT/run_${TEST_NAME}.log -e $PATHRT/run_${TEST_NAME}.log -q $QUEUE  -A $ACCOUNT \
	    -Wblock=true -l walltime=00:${WLCLK}:00 -N $TEST_NAME -l select=1:ncpus=1:mem=12GB -v RESNAME=$TEST_NAME $SBATCH_COMMAND)

	#   qsub -o $PATHRT/run_${TEST_NAME}.log -e $PATHRT/run_${TEST_NAME}.log -q $QUEUE  -A $ACCOUNT \
	    # -Wblock=true -l walltime=00:01:00 -N chgres_summary -l select=1:ncpus=1:mem=100MB -W depend=afternotok:$TEST << EOF
	#!/bin/bash
	#   touch $RUNDIR/bad.${TEST_NAME}
	#EOF
	#   if [[ -f $RUNDIR/bad.${TEST_NAME} ]]; then
	#     error "Batch job for test $TEST_NAME did not finish successfully. Refer to run_${TEST_NAME}.log"
	#   fi

    else
	sbatch --wait --ntasks-per-node=1 --nodes=1 ${MEM} -t 00:${WLCLK}:00 -A $ACCOUNT -q $QUEUE -J $TEST_NAME \
	    --partition=$PARTITION -o $PATHRT/run_${TEST_NAME}.log -e $PATHRT/run_${TEST_NAME}.log \
	    --wrap "time $SBATCH_COMMAND $TEST_NAME" && d=$? || d=$?

	if [[ d -ne 0 ]]; then
	    error "Batch job for test $TEST_NAME did not finish successfully. Refer to run_${TEST_NAME}.log"
	fi

    fi

    check_results

done <$TESTS_FILE
if [[ $? -ne 0 ]]; then
    error "Run test while loop did not finish properly"
fi

cd $PATHRT
FAIL_FILES="fail_test_*"
for file in $FAIL_FILES; do
    if [[ -f "$file" ]]; then
	cat "$file" >> fail_test
    fi
done

if [[ -e fail_test ]]; then
    echo | tee -a $REGRESSIONTEST_LOG
    for file in fail_test_*; do
	cat $file >>$REGRESSIONTEST_LOG
	cat $file >>summary.log
    done

    echo | tee -a $REGRESSIONTEST_LOG
    echo "REGRESSION TEST FAILED" | tee -a $REGRESSIONTEST_LOG
else
    echo | tee -a $REGRESSIONTEST_LOG
    echo "REGRESSION TEST WAS SUCCESSFUL" | tee -a $REGRESSIONTEST_LOG
    echo "All tests passed" >>summary.log
fi
date >> $REGRESSIONTEST_LOG

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )
echo "Elapsed time: ${elapsed_time}. Have a nice day!" >> ${REGRESSIONTEST_LOG}
echo "Elapsed time: ${elapsed_time}. Have a nice day!"
