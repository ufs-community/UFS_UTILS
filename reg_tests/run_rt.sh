#!/bin/bash

export MAILTO="wei.huang@noaa.edu"
export WORK_DIR=/contrib/Wei.Huang/dev/UFS_UTILS
export PROJECT_CODE=Wei.Huang
export QUEUE=compute
TIMEOUT_LIMIT=3600
#repo=https://github.com/ufs-community/UFS_UTILS.git
repo=git@github.com:NOAA-EPIC/UFS_UTILS-cloud.git
branch=feature/reg_tests_on_cloud

rt.sh -w ${WORK_DIR} -m ${MAILTO} -r ${repo} -b ${branch} -v
