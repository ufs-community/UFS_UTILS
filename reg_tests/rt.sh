#!/bin/bash

wait_for_fin() {
  set -xu
  sleep_time=0
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       mail -s "UFS_UTILS Consistency Tests timed out on ${MACHINE_ID}" "${MAILTO}" < "${WORK_DIR}/reg_test_results.txt"
       exit 1
    fi
  done
}

if [[ $# -gt 0 ]] && [[ $# -lt 2 ]]; then
    SCRIPT_LOGGER=$1
else
    SCRIPT_LOGGER="rtsh.log"
fi
exec 2> "${SCRIPT_LOGGER}"

start_time=$SECONDS
if [[ "$(hostname)" =~ "Orion" || "$(hostname)" =~ "orion" ]]; then
  ulimit -a
else
  ulimit -s unlimited
fi

if [[ "$(hostname)" =~ "WCOSS2" ]]; then
  compiler="intel"
else
  compiler="intelllvm"
fi
export compiler

# shellcheck source=./rt.control
RT_DIR=${PWD}
export RT_DIR
source "${RT_DIR}/rt.control"


mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || { echo "Can't change directory to '${WORK_DIR}'.. exiting"; exit; }
rm -f reg_test_results.txt
rm -rf UFS_UTILS

git clone -b "${REPO_BRANCH}" "${REPO_LOC}"
rc=$?

### Check to see if the clone was successful.
### Previously, it has failed due to lack of disk space.

if [[ $rc == 0 ]] && [[ -d UFS_UTILS ]];then
  echo "Clone Successful"
else
  echo "Clone Failed" | mail -s "UFS_UTILS Consistency Tests failed on ${MACHINE_ID}" "${MAILTO}"
  exit
fi
###

cd UFS_UTILS || { echo "Can't change directory into 'UFS_UTILS'.. exiting"; exit; }

# shellcheck source=../sorc/machine-setup.sh
source sorc/machine-setup.sh

current_hash=$(git rev-parse HEAD)

if [[ -f "${WORK_DIR}/prev_hash.txt" ]]; then
    prev_hash=$(cat "${WORK_DIR}/prev_hash.txt")
    if [[ "${current_hash}" == "${prev_hash}" ]]; then
        date
        echo ""
        echo "UFS_UTILS has not changed since last time. Not building."
        echo "UFS_UTILS hash: ${current_hash}"
        exit 0
    fi
fi

echo "Started on $(hostname -s)" >> "${WORK_DIR}/reg_test_results.txt"

./build_all.sh

if [[ ${MACHINE_ID} == "wcoss2" ]]; then
    this_machine=$(cat /etc/cluster_name)
    prod_machine=$(grep primary /lfs/h1/ops/prod/config/prodmachinefile)
    prod_machine=${prod_machine/primary:}
    if [[ "${this_machine}" == "${prod_machine}" ]]; then
        echo "ERROR: Cannot use the production machine"; exit 0
    fi
fi

cd fix || { echo "Can't change directory into 'fix'.. exiting"; exit; }

./link_fixdirs.sh emc "${MACHINE_ID}"

cd ../reg_tests || { echo "Can't change directory into '../reg_tests'.. exiting"; exit; }

set -x
PID_LIST=()
notlocal=true
export notlocal
export ACCOUNT=$PROJECT_CODE
export STMP=$WORK_DIR/reg-tests

for dir in ocnice_prep cpld_gridgen; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    echo "Running ${dir} tests"
    cd ${dir} || { echo "Can't change directory into '${dir}'.. exiting"; exit; }
    if [[ ${UPDATE_BASELINE} == "TRUE" ]]; then
        (bash ./rt.sh -c > ${dir}_rt.out 2>&1) &
    else
      (bash ./rt.sh > ${dir}_rt.out 2>&1) &
    fi
  PID_LIST+=($!)
  cd ..
  fi
done

for dir in regrid_sfc global_cycle chgres_cube snow2mdl grid_gen weight_gen ice_blend; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    echo "Running ${dir} tests"
    cd "${dir}" || { echo "Can't change directory into '${dir}'.. exiting"; exit; }
    (bash "./driver.sh" > "${dir}_rt.out" 2>&1) &
    PID_LIST+=($!)
    cd ..
  fi
done

echo "SUBMITTED ALL TASKS. Waiting for them to finish.."
wait "${PID_LIST[@]}"

echo "Commit hash: ${current_hash}" >> "${WORK_DIR}/reg_test_results.txt"
echo "" >> "${WORK_DIR}/reg_test_results.txt"

success=true
for dir in regrid_sfc weight_gen ocnice_prep cpld_gridgen chgres_cube grid_gen global_cycle ice_blend snow2mdl; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    if [[ ! -f ${dir}/summary.log ]]; then
        success=false
        echo "${dir} consistency tests FAILED: no summary.log" >> "${WORK_DIR}/reg_test_results.txt"
    elif grep -qi "FAILED" ${dir}/summary.log; then
        success=false
        echo "${dir} consistency tests FAILED" >> "${WORK_DIR}/reg_test_results.txt"
    else
        echo "${dir} consistency tests PASSED" >> "${WORK_DIR}/reg_test_results.txt"
    fi
  fi
done

end_time=$SECONDS
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds" >> "${WORK_DIR}/reg_test_results.txt"
echo "Finished on ${MACHINE_ID}" >> "${WORK_DIR}/reg_test_results.txt"

# shellcheck disable=SC2086
if [[ "$success" == true ]]; then
    mail -s "UFS_UTILS Consistency Tests PASSED on ${MACHINE_ID}" ${MAILTO} < "${WORK_DIR}/reg_test_results.txt"
else
    mail -s "UFS_UTILS Consistency Tests FAILED on ${MACHINE_ID}" ${MAILTO} < "${WORK_DIR}/reg_test_results.txt"
fi

# Save current hash as previous hash for next time
echo "${current_hash}" > "${WORK_DIR}/prev_hash.txt"


