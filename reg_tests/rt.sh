#!/bin/bash

export MAILTO=

# Directory to download UFS_UTILS to and run the regression tests
export WORK_DIR=

export PROJECT_CODE=
export QUEUE=

mkdir -p ${WORK_DIR}
cd ${WORK_DIR}
rm -rf UFS_UTILS

git clone --recursive https://github.com/NOAA-EMC/UFS_UTILS.git
cd UFS_UTILS

source sorc/machine-setup.sh

current_hash=$(git rev-parse HEAD)

if [[ -f "${WORK_DIR}/prev_hash.txt" ]]; then
    prev_hash=$(cat ${WORK_DIR}/prev_hash.txt)
    if [[ "$current_hash" == "$prev_hash" ]]; then
        echo `date`
        echo ""
        echo "UFS_UTILS has not changed since last time. Not building."
        echo "UFS_UTILS hash: ${current_hash}"
        exit 0
    fi
fi

./build_all.sh

if [[ $target == "wcoss_dell_p3" ]] || [[ $target == "wcoss_cray" ]]; then
    prod_machine=`cat /etc/prod`
    prod_letter=${prod_machine:0:1}

    this_machine=`hostname -s`
    this_letter=${this_machine:0:1}

    # Mars (m), Venus (v)
    if [[ "${this_letter}" == "${prod_letter}" ]]; then
        exit 0
    fi
fi

# Set machine_id variable for running link_fixdirs
if [[ $target == "wcoss_dell_p3" ]]; then
    machine_id=dell
elif [[ $target == "wcoss_cray" ]]; then
    machine_id=cray
else
    machine_id=$target
fi

cd fix
./link_fixdirs.sh emc $machine_id

cd ../reg_tests

for dir in chgres_cube grid_gen; do
    cd $dir
    ./driver.$target.sh
    # Wait for job to complete
    while [ ! -f "summary.log" ]; do
        sleep 10
    done
    cd ..
done

for dir in snow2mdl global_cycle ice_blend; do
    cd $dir
    if [[ $target == "hera" ]] || [[ $target == "jet" ]] || [[ $target == "orion" ]]; then
        sbatch -A ${PROJECT_CODE} ./driver.$target.sh
    elif [[ $target == "wcoss_dell_p3" ]] || [[ $target == "wcoss_cray" ]]; then
        cat ./driver.$target.sh | bsub -P ${PROJECT_CODE}
    fi
    # Wait for job to complete
    while [ ! -f "summary.log" ]; do
        sleep 10
    done
    cd ..
done

# Wait chgres_cube and grid_gen to finish before submitting more jobs
while [ ! -f "snow2mdl/summary.log" ] || [ ! -f "global_cycle/summary.log" ] || [ ! -f "ice_blend/summary.log" ]; do
    sleep 10
done

echo "Commit hash: ${current_hash}" >> reg_test_results.txt
echo "" >> reg_test_results.txt

for dir in chgres_cube grid_gen global_cycle ice_blend snow2mdl; do
    success=true
    if grep -qi "FAILED" ${dir}/summary.log; then
        success=false
        echo "${dir} regression tests FAILED" >> reg_test_results.txt
    else
        echo "${dir} regression tests PASSED" >> reg_test_results.txt
    fi
done

if [[ "$success" == true ]]; then
    mail -s "UFS_UTILS Regression Tests PASSED on ${target}" ${MAILTO} < reg_test_results.txt
else
    mail -s "UFS_UTILS Regression Tests FAILED on ${target}" ${MAILTO} < reg_test_results.txt
fi

# Save current hash as previous hash for next time
echo $current_hash > ${WORK_DIR}/prev_hash.txt




