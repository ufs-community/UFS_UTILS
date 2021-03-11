#!/bin/bash

export MAILTO="kyle.gerheiser@noaa.gov"

export UFS_UTILS_WORKING_DIR=/home/Kyle.Gerheiser/UFS_UTILS-cron
export UFS_UTILS_HOME_DIR=$PWD/..

mkdir -p ${UFS_UTILS_WORKING_DIR}
cd ${UFS_UTILS_WORKING_DIR}
rm -rf UFS_UTILS

export MACHINE_ID=hera

cd ${UFS_UTILS_WORKING_DIR}

git clone --recursive https://github.com/NOAA-EMC/UFS_UTILS.git >> /dev/null 2>&1
cd UFS_UTILS

source sorc/machine-setup.sh

current_hash=$(git rev-parse HEAD)

if [[ -f "${UFS_UTILS_WORKING_DIR}/prev_hash.txt" ]]; then
    prev_hash=$(cat ${HPC_HOMEDIR}/prev_hash.txt)
    if [[ "$current_hash" == "$prev_hash" ]]; then
        echo `date`
        echo ""
        echo "UFS_UTILS has not changed since last time. Not building."
        echo "UFS_UTILS hash: ${current_hash}"
        exit 0
    fi
fi

./build_all.sh

cd fix
./link_fixdirs.sh emc $MACHINE_ID

cd ../reg_tests

for dir in chgres_cube grid_gen ice_blend snow2mdl; do
    cd $dir
    ./driver.$target.sh
    cd ..
done


time=0
should_wait=true
while [ ${should_wait} == "true" ]; do
    sleep 10
    for dir in chgres_cube grid_gen ice_blend snow2mdl; do
        should_wait=false
        if [[ ! -f "${dir}/summary.log" ]]; then
            should_wait=true
            break
        fi
    done
done

for dir in chgres_cube grid_gen ice_blend snow2mdl; do
    if grep -qi "FAILED" ${dir}/summary.log; then
        echo "${dir} regression tests FAILED" >> reg_test_results.txt
        mail -s "UFS_UTILS Regression Tests FAILED" ${MAILTO} < reg_test_results.txt
    else
        echo "${dir} regression tests PASSED" >> reg_test_results.txt
        mail -s "UFS_UTILS Regression Tests PASSED" ${MAILTO} < reg_test_results.txt
    fi
done


