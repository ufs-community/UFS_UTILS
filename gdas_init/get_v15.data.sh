#!/bin/bash

bundle=$1

set -x

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

if [ $bundle = 'hires' ]; then

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
  file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.gdas_restart.tar

  rm -f ./list.hires.*
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
  file=gpfs_dell1_nco_ops_com_gfs_prod_gdas.${yy}${mm}${dd}_${hh}.gdas.tar

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

  rm -f ./list.hires.*

  touch ./gdas.${yy}${mm}${dd}/${hh}/gdas.t${hh}z.loginc.txt

  set +x
  echo DATA PULL FOR $bundle DONE

  exit 0

fi

#----------------------------------------------------------------------
# Get the enkf tiled restart files for all members.
#----------------------------------------------------------------------

#mkdir -p $EXTRACT_DIR/$bundle
#cd $EXTRACT_DIR/$bundle

for group in $bundle
do

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
  file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${group}.tar

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

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy}${mm}${dd}_${hh}.enkfgdas_restart_${group}.tar

  rm -f ./list*.${group}
  htar -tvf  $directory/$file > ./list1.${group}
  grep radstat ./list1.${group} > ./list2.${group}
  grep abias ./list1.${group} >> ./list2.${group}
    while read -r line
  do
    echo ${line##*' '} >> ./list3.${group}
  done < "./list2.${group}"

  htar -xvf $directory/$file -L ./list3.${group}
  rc=$?
  [ $rc != 0 ] && exit $rc

  case $bundle in
    grp1)
     members="mem001 mem002 mem003 mem004 mem005 mem006 mem007 mem008 mem009 mem010"
     ;;
    grp2)
     members="mem011 mem012 mem013 mem014 mem015 mem016 mem017 mem018 mem019 mem020"
     ;;
    grp3)
     members="mem021 mem022 mem023 mem024 mem025 mem026 mem027 mem028 mem029 mem030"
     ;;
    grp4)
     members="mem031 mem032 mem033 mem034 mem035 mem036 mem037 mem038 mem039 mem040"
     ;;
    grp5)
     members="mem041 mem042 mem043 mem044 mem045 mem046 mem047 mem048 mem049 mem050"
     ;;
    grp6)
     members="mem051 mem052 mem053 mem054 mem055 mem056 mem057 mem058 mem059 mem060"
     ;;
    grp7)
     members="mem061 mem062 mem063 mem064 mem065 mem066 mem067 mem068 mem069 mem070"
     ;;
    grp8)
     members="mem071 mem072 mem073 mem074 mem075 mem076 mem077 mem078 mem079 mem080"
     ;;
  esac

  for member in $members
  do
    touch ./enkfgdas.${yy}${mm}${dd}/${hh}/${member}/enkfgdas.t${hh}z.loginc.txt
  done

  rm -f ./list*.${group}

done

set +x
echo DATA PULL FOR $bundle DONE

exit 0
