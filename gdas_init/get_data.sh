#!/bin/bash

set -x

mkdir -p $EXTRACT_DIR
cd $EXTRACT_DIR

date10=`$NDATE -6 $yy$mm$dd$hh`

echo $date10
yy=$(echo $date10 | cut -c1-4)
mm=$(echo $date10 | cut -c5-6)
dd=$(echo $date10 | cut -c7-8)
hh=$(echo $date10 | cut -c9-10)

directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}

# the hires file
file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy}${mm}${dd}_${hh}.gdas_restart.tar
htar -xvf $directory/$file

for group in 'grp1' 'grp2' 'grp3' 'grp4' 'grp5' 'grp6' 'grp7' 'grp8'
do
  file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy}${mm}${dd}_${hh}.enkfgdas_restart_${group}.tar
  htar -xvf $directory/$file 
done

echo DONE

exit 0
