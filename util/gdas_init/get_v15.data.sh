#!/bin/bash

#----------------------------------------------------------------------
# Retrieve gfs v15 data from hpss.
#
# Data available after 2019061206.
#----------------------------------------------------------------------

bundle=$1

set -x

cd $EXTRACT_DIR

date10_m6=`$NDATE -6 $yy$mm$dd$hh`

echo $date10_m6
yy_m6=$(echo $date10_m6 | cut -c1-4)
mm_m6=$(echo $date10_m6 | cut -c5-6)
dd_m6=$(echo $date10_m6 | cut -c7-8)
hh_m6=$(echo $date10_m6 | cut -c9-10)

#----------------------------------------------------------------------
# Get the hires tiled restart files.  Need to use the 6-hour forecast files from
# the previous cycle.  
#----------------------------------------------------------------------

if [ $bundle = 'hires' ]; then

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
  if [ $date10_m6 -lt 2020022600 ]; then
    file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.gdas_restart.tar
  else
    file=com_gfs_prod_gdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.gdas_restart.tar
  fi

  rm -f ./list.hires*
  touch ./list.hires3
  htar -tvf  $directory/$file > ./list.hires1
  grep ${yy}${mm}${dd}.${hh} ./list.hires1 > ./list.hires2
  while read -r line
  do 
    echo ${line##*' '} >> ./list.hires3
  done < "./list.hires2"

  htar -xvf $directory/$file -L ./list.hires3
  rc=$?
  [ $rc != 0 ] && exit $rc

#----------------------------------------------------------------------
# Get the 'abias' and 'radstat' files from current cycle
#----------------------------------------------------------------------

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  if [ ${yy}${mm}${dd}${hh} -lt 2020022600 ]; then
    file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy}${mm}${dd}_${hh}.gdas.tar
  else
    file=com_gfs_prod_gdas.${yy}${mm}${dd}_${hh}.gdas.tar
  fi

  htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.radstat
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_air
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_int
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.abias_pc
  rc=$?
  [ $rc != 0 ] && exit $rc

  rm -f ./list.hires*

#----------------------------------------------------------------------
# Get the enkf tiled restart files for all members.
#----------------------------------------------------------------------

else

  for group in $bundle
  do

    directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
    if [ $date10_m6 -lt 2020022600 ]; then
      file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${group}.tar
    else
      file=com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${group}.tar
    fi

    rm -f ./list*.${group}
    htar -tvf  $directory/$file > ./list1.${group}
    grep ${yy}${mm}${dd}.${hh} ./list1.${group} > ./list2.${group}
    while read -r line
    do 
      echo ${line##*' '} >> ./list3.${group}
    done < "./list2.${group}"
    htar -xvf $directory/$file  -L ./list3.${group}
    rc=$?
    [ $rc != 0 ] && exit $rc
    rm -f ./list*.${group}

  done

fi

set +x
echo DATA PULL FOR $bundle DONE

exit 0
