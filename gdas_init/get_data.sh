#!/bin/bash

if [ $yy$mm$dd$hh -lt 2019061206 ]; then
  set +x
  echo "ERROR: SCRIPTS DO NOT SUPPORT PRE-GFS.V15 DATA"
  echo "CHOOSE DATE STARTING 2019061206"
  exit 1
fi

set -x

rm -fr $EXTRACT_DIR
mkdir -p $EXTRACT_DIR
cd $EXTRACT_DIR

date10=`$NDATE -6 $yy$mm$dd$hh`

echo $date10
yy_m6=$(echo $date10 | cut -c1-4)
mm_m6=$(echo $date10 | cut -c5-6)
dd_m6=$(echo $date10 | cut -c7-8)
hh_m6=$(echo $date10 | cut -c9-10)

#----------------------------------------------------------------------
# Get the hires tiled restart files.  Need to use the 6-hour forecast files from
# the previous cycle.  
#----------------------------------------------------------------------

directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.gdas_restart.tar

rm -f ./list*
touch ./list3
htar -tvf  $directory/$file > ./list1
grep ${yy}${mm}${dd}.${hh} ./list1 > ./list2
while read -r line
do 
  echo ${line##*' '} >> ./list3
done < "./list2"

htar -xvf $directory/$file -L ./list3

#----------------------------------------------------------------------
# Get the 'abias' and 'radstat' files from current cycle
#----------------------------------------------------------------------

directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy}${mm}${dd}_${hh}.gdas.tar

htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.radstat
htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias
htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_air
htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_int
htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_pc

for group in 'grp1' 'grp2' 'grp3' 'grp4' 'grp5' 'grp6' 'grp7' 'grp8'
do

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
  file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${group}.tar

  rm -f ./list*
  htar -tvf  $directory/$file > ./list1
  grep ${yy}${mm}${dd}.${hh} ./list1 > ./list2
  while read -r line
  do 
    echo ${line##*' '} >> ./list3
  done < "./list2"
  htar -xvf $directory/$file  -L ./list3

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy}${mm}${dd}_${hh}.enkfgdas_restart_${group}.tar

  rm -f ./list*
  htar -tvf  $directory/$file > ./list1
  grep radstat ./list1 > ./list2
  grep abias ./list1 >> ./list2
    while read -r line
  do
    echo ${line##*' '} >> ./list3
  done < "./list2"

  htar -xvf $directory/$file -L ./list3

done

echo DONE

exit 0
