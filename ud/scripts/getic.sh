#!/bin/ksh

#########################################################################
# get global ic from info under
# https://github.com/NOAA-EMC/global-workflow/wiki/Run-Global-Workflow#warm-starts-from-production
##################################

set -ex

YYYY=${1}
YYYYMM=${YYYY}${2}
YYYYMMDD=${YYYYMM}${3} 
CC=${4}

dir=/NCEPPROD/hpssprod/runhistory
com=gpfs_dell1_nco_ops_com_gfs_prod
#mdl=gfs
 mdl=gdas

htar -xvf ${dir}/rh${YYYY}/${YYYYMM}/${YYYYMMDD}/${com}_${mdl}.${YYYYMMDD}_${CC}.${mdl}_restart.tar
