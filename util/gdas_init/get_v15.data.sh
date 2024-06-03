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
# Because the use of nemsio data is being phased out from chgres_cube,
# use the GDAS tiled restart files (which are netcdf) for both
# the GDAS high-res and the GFS free forecast runs.
#
# Note: Need to use the 6-hour forecast files from the previous
# cycle as they are not saved at the current cycle.
#----------------------------------------------------------------------

if [ "$bundle" == "gdas" ] || [ "$bundle" == "gfs" ] ; then

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

  rm -f ./list.hires*

  [ "$bundle" == "gfs" ] && exit 0

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

#----------------------------------------------------------------------
# Get the enkf tiled restart files for all members.  They are not
# stored for the current cycle, so use the 6-hr old tarball.
#----------------------------------------------------------------------

else

    directory=/NCEPPROD/hpssprod/runhistory/rh${yy_m6}/${yy_m6}${mm_m6}/${yy_m6}${mm_m6}${dd_m6}
    if [ $date10_m6 -lt 2020022600 ]; then
      file=gpfs_dell1_nco_ops_com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${bundle}.tar
    else
      file=com_gfs_prod_enkfgdas.${yy_m6}${mm_m6}${dd_m6}_${hh_m6}.enkfgdas_restart_${bundle}.tar
    fi

    rm -f ./list*.${bundle}
    htar -tvf  $directory/$file > ./list1.${bundle}
    grep ${yy}${mm}${dd}.${hh} ./list1.${bundle} > ./list2.${bundle}
    while read -r line
    do 
      echo ${line##*' '} >> ./list3.${bundle}
    done < "./list2.${bundle}"
    htar -xvf $directory/$file  -L ./list3.${bundle}
    rc=$?
    [ $rc != 0 ] && exit $rc
    rm -f ./list*.${bundle}

fi

set +x
echo DATA PULL FOR $bundle DONE

exit 0
