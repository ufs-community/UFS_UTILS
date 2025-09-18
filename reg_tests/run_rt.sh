#!/bin/bash

export MAILTO="wei.huang@noaa.edu"
#repo=https://github.com/ufs-community/UFS_UTILS.git
repo=git@github.com:NOAA-EPIC/UFS_UTILS-cloud.git
branch=feature/reg_tests_on_cloud

target=`hostname -s`
if [[ -d /lfs3 ]] ; then
  target=Jet
elif [[ -d /lfs/h1 ]] ; then
  target=WCOSS2
elif [[ -d /scratch3 ]] ; then
  target=Ursa
elif [[ -d /contrib ]] ; then
  target=noaacloud
fi
echo "Run rt.sh on ${target}"

if [[ $target == "Ursa" ]]; then
  export WORK_DIR=/scratch4/NAGAPE/epic/${USER}/debug
  export PROJECT_CODE=epic
  export QUEUE=batch
elif [[ $target == "Ursa" ]]; then
  export WORK_DIR=/contrib/Wei.Huang/dev/UFS_UTILS
  export PROJECT_CODE=${USER}
  export QUEUE=compute
fi

TIMEOUT_LIMIT=3600

rt.sh -w ${WORK_DIR} -m ${MAILTO} -r ${repo} -p ${PROJECT_CODE} -q ${QUEUE} -b ${branch} -v
