#!/bin/bash -l

ulimit -s unlimited

export MAILTO=george.gayno@noaa.gov

# Directory to download UFS_UTILS to and run the consistency tests
export WORK_DIR=/lfs/h2/emc/stmp/george.gayno/reg.tests.branch

export PROJECT_CODE=GFS-DEV
export QUEUE=dev
TIMEOUT_LIMIT=21600

mkdir -p ${WORK_DIR}
cd ${WORK_DIR}
rm -f reg_test_results.txt
rm -rf UFS_UTILS

git clone --recursive https://github.com/GeorgeGayno-NOAA/UFS_UTILS.git
cd UFS_UTILS
git checkout feature/new_fix

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

echo "Started on " `hostname -s` >> ${WORK_DIR}/reg_test_results.txt

./build_all.sh

if [[ $target == "wcoss2" ]]; then
    this_machine=`cat /etc/cluster_name`
    prod_machine=`grep primary /lfs/h1/ops/prod/config/prodmachinefile`
    prod_machine=`echo ${prod_machine/primary:}`
    if [[ "${this_machine}" == "${prod_machine}" ]]; then
        exit 0
    fi
fi

machine_id=$target

cd fix
./link_fixdirs.sh emc $machine_id

cd ../reg_tests

#if [[ $target == "orion" ]] || [[ $target == "jet" ]] || [[ $target == "hera" ]] || [[ $target == "wcoss2" ]] ; then
if [[ $target == "orion" ]] || [[ $target == "jet" ]] || [[ $target == "hera" ]] ; then

  cd cpld_gridgen
  export ACCOUNT=$PROJECT_CODE
  export STMP=$WORK_DIR/reg-tests

  ./rt.sh 2>/dev/null &

  set -x

  sleep_time=0
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       kill -9 %1
       mail -s "UFS_UTILS Consistency Tests timed out on ${target}" ${MAILTO} < ${WORK_DIR}/reg_test_results.txt
       exit 1
    fi
  done
  cd ..

fi

sleep_time=0
for dir in snow2mdl global_cycle chgres_cube grid_gen; do
    cd $dir
    ./driver.$target.sh
    # Wait for job to complete
    while [ ! -f "summary.log" ]; do
        sleep 10
        sleep_time=$((sleep_time+10))
        if (( sleep_time > TIMEOUT_LIMIT )); then
             mail -s "UFS_UTILS Consistency Tests timed out on ${target}" ${MAILTO} < ${WORK_DIR}/reg_test_results.txt
            exit 1
        fi
    done
    cd ..
done

for dir in ice_blend; do
    cd $dir
    if [[ $target == "hera" ]] || [[ $target == "jet" ]] || [[ $target == "orion" ]] || [[ $target == "s4" ]] ; then
        sbatch -A ${PROJECT_CODE} ./driver.$target.sh
    elif [[ $target == "wcoss2" ]] ; then
        qsub -v WORK_DIR ./driver.$target.sh
    fi
    
    # Wait for job to complete
    sleep_time=0
    while [ ! -f "summary.log" ]; do
        sleep 10
        sleep_time=$((sleep_time+10))
        if (( sleep_time > TIMEOUT_LIMIT )); then
            mail -s "UFS_UTILS Consistency Tests timed out on ${target}" ${MAILTO} < ${WORK_DIR}/reg_test_results.txt
            exit 1
        fi
    done
    cd ..
done


echo "Commit hash: ${current_hash}" >> ${WORK_DIR}/reg_test_results.txt
echo "" >> ${WORK_DIR}/reg_test_results.txt

success=true
for dir in cpld_gridgen chgres_cube grid_gen global_cycle ice_blend snow2mdl; do
    if grep -qi "FAILED" ${dir}/summary.log; then
        success=false
        echo "${dir} consistency tests FAILED" >> ${WORK_DIR}/reg_test_results.txt
    else
        echo "${dir} consistency tests PASSED" >> ${WORK_DIR}/reg_test_results.txt
    fi
done

if [[ "$success" == true ]]; then
    mail -s "UFS_UTILS Consistency Tests PASSED on ${target}" ${MAILTO} < ${WORK_DIR}/reg_test_results.txt
else
    mail -s "UFS_UTILS Consistency Tests FAILED on ${target}" ${MAILTO} < ${WORK_DIR}/reg_test_results.txt
fi

# Save current hash as previous hash for next time
echo $current_hash > ${WORK_DIR}/prev_hash.txt
