#!/bin/bash

wait_for_fin() {
  set -xu
  sleep_time=0
  while [ ! -f "summary.log" ]; do
    sleep 10
    sleep_time=$((sleep_time+10))
    if (( sleep_time > TIMEOUT_LIMIT )); then
       mail -s "UFS_UTILS Consistency Tests timed out on ${target}" "${MAILTO}" < "${WORK_DIR}/reg_test_results.txt"
       exit 1
    fi
  done
}

if [[ "$(hostname)" =~ "Orion" || "$(hostname)" =~ "orion" ]]; then
  ulimit -a
else
  ulimit -s unlimited
fi

# export MAILTO=

# # Directory to download UFS_UTILS to and run the consistency tests
# export WORK_DIR=

# export PROJECT_CODE=

# export QUEUE=

# shellcheck source=./rt.control
source ./rt.control

# TIMEOUT_LIMIT=3600

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || { echo "Can't change directory to '${WORK_DIR}'.. exiting"; exit; }
rm -f reg_test_results.txt
rm -rf UFS_UTILS

#git clone https://github.com/ufs-community/UFS_UTILS.git
git clone "${REPO_TO_RUN}"
rc=$?

# Check to see if the clone was successful. Previously, it has
# failed due to lack of disk space.

if [[ $rc == 0 ]] && [[ -d UFS_UTILS ]];then
  echo "Clone Successful"
else
  target=$(hostname -s)
  if [[ -d /lfs3 ]] ; then
    target=Jet
  elif [[ -d /lfs/h1 ]] ; then
    target=WCOSS2
  elif [[ -d /scratch3 ]] ; then
    target=Ursa
  fi
  echo "Clone Failed" | mail -s "UFS_UTILS Consistency Tests failed on ${target}" "${MAILTO}"
fi

cd UFS_UTILS || { echo "Can't change directory into 'UFS_UTILS'.. exiting"; exit; }

# shellcheck source=../sorc/machine-setup.sh
source ../sorc/machine-setup.sh

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

if [[ ${target} == "wcoss2" ]]; then
    this_machine=$(cat /etc/cluster_name)
    prod_machine=$(grep primary /lfs/h1/ops/prod/config/prodmachinefile)
    prod_machine=${prod_machine/primary:}
    if [[ "${this_machine}" == "${prod_machine}" ]]; then
        echo "ERROR: Cannot use the production machine"; exit 0
    fi
fi

machine_id=${target}

cd fix || { echo "Can't change directory into 'fix'.. exiting"; exit; }

./link_fixdirs.sh emc "${machine_id}"

cd ../reg_tests || { echo "Can't change directory into '../reg_tests'.. exiting"; exit; }

set -x

if [[ " ${RUN_SET[*]} " =~ " RUN_REGRID_SFC " ]]; then
  cd regrid_sfc || { echo "Can't change directory into 'regrid_sfc'.. exiting"; exit; }
  ./driver.sh && wait_for_fin &
  cd ..
fi
export ACCOUNT=$PROJECT_CODE
export STMP=$WORK_DIR/reg-tests

if [[ " ${RUN_SET[*]} " =~ " RUN_OCNICE_PREP " ]]; then
  cd ocnice_prep || { echo "Can't change directory into 'ocnice_prep'.. exiting"; exit; }
  ./rt.sh && wait_for_fin &
  cd ..
fi

if [[ " ${RUN_SET[*]} " =~ " RUN_CPLD_GRIDDEN " ]]; then
  cd cpld_gridgen || { echo "Can't change directory into 'cpld_gridgen'.. exiting"; exit; }
  ./rt.sh && wait_for_fin &
  cd ..
fi

for dir in snow2mdl global_cycle chgres_cube grid_gen; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    cd "${dir}" || { echo "Can't change directory into '${dir}'.. exiting"; exit; }
    bash "./driver.${target}.sh" && wait_for_fin &
    cd ..
  fi
done

for dir in weight_gen ice_blend; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    cd "${dir}" || { echo "Can't change directory into '${dir}'.. exiting"; exit; }
    if [[ ${target} == "ursa" ]] || [[ ${target} == "jet" ]] || [[ ${target} == "orion" ]] || [[ ${target} == "hercules" ]] ; then
        sbatch -A "${PROJECT_CODE}" "./driver.${target}.sh" && wait_for_fin &
    elif [[ ${target} == "wcoss2" ]] ; then
        qsub -v WORK_DIR "./driver.${target}.sh" && wait_for_fin &
    fi
    cd ..
  fi
done

wait

echo "Commit hash: ${current_hash}" >> "${WORK_DIR}/reg_test_results.txt"
echo "" >> "${WORK_DIR}/reg_test_results.txt"

success=true
for dir in regrid_sfc weight_gen ocnice_prep cpld_gridgen chgres_cube grid_gen global_cycle ice_blend snow2mdl; do
  RUN_CHECK=RUN_${dir^^}
  if [[ " ${RUN_SET[*]} " =~ ${RUN_CHECK} ]]; then
    if grep -qi "FAILED" ${dir}/summary.log; then
        success=false
        echo "${dir} consistency tests FAILED" >> "${WORK_DIR}/reg_test_results.txt"
    else
        echo "${dir} consistency tests PASSED" >> "${WORK_DIR}/reg_test_results.txt"
    fi
  fi
done

if [[ "$success" == true ]]; then
    mail -s "UFS_UTILS Consistency Tests PASSED on ${target}" "${MAILTO}" < "${WORK_DIR}/reg_test_results.txt"
else
    mail -s "UFS_UTILS Consistency Tests FAILED on ${target}" "${MAILTO}" < "${WORK_DIR}/reg_test_results.txt"
fi

# Save current hash as previous hash for next time
echo "${current_hash}" > "${WORK_DIR}/prev_hash.txt"
