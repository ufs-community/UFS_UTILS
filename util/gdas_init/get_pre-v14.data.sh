#!/bin/bash

#-----------------------------------------------------------
# Retrieve data prior to v14 (the old sigio/sfcio data
# from the spectral gfs).
#
# Script works for data between 00z May 21, 2012
# and 06z July 19, 2017.
#-----------------------------------------------------------

bundle=$1

set -x

if [ $bundle = 'hires' ]; then

  mkdir -p $EXTRACT_DIR/gdas.${yy}${mm}${dd}/${hh}
  cd $EXTRACT_DIR/gdas.${yy}${mm}${dd}/${hh}

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  if [ $gfs_ver = 'v12' ]; then
    file=com_gfs_prod_gdas.${yy}${mm}${dd}${hh}.tar
  else
    file=com2_gfs_prod_gdas.${yy}${mm}${dd}${hh}.tar
  fi

  htar -xvf $directory/$file ./gdas1.t${hh}z.radstat
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas1.t${hh}z.abias_air
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas1.t${hh}z.abias
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas1.t${hh}z.abias_pc
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas1.t${hh}z.sanl
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas1.t${hh}z.sfcanl
  rc=$?
  [ $rc != 0 ] && exit $rc

elif [ $bundle = 'enkf' ]; then

#----------------------------------------------------------------------
# Get the enkf tiled restart files for all members.
#----------------------------------------------------------------------

  mkdir -p $EXTRACT_DIR/enkf.${yy}${mm}${dd}/${hh}
  cd $EXTRACT_DIR/enkf.${yy}${mm}${dd}/${hh}

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  if [ $gfs_ver = 'v12' ]; then
    file=com_gfs_prod_enkf.${yy}${mm}${dd}_${hh}.anl.tar
  else
    file=com2_gfs_prod_enkf.${yy}${mm}${dd}_${hh}.anl.tar
  fi

  rm -f ./list*.${bundle}
  htar -tvf $directory/$file > ./list1.${bundle}
  grep siganl ./list1.${bundle} > ./list2.${bundle}
  grep sfcanl ./list1.${bundle} >> ./list2.${bundle}
  while read -r line
  do
    echo ${line##*' '} >> ./list3.${bundle}
  done < "./list2.${bundle}"
  htar -xvf $directory/$file  -L ./list3.${bundle}
  rc=$?
  [ $rc != 0 ] && exit $rc

  MEMBER=1
  while [ $MEMBER -le 80 ]; do
    if [ $MEMBER -lt 10 ]; then
      MEMBER_CH="00${MEMBER}"
    else
      MEMBER_CH="0${MEMBER}"
    fi
    mkdir -p mem${MEMBER_CH}
    mv *_mem${MEMBER_CH}* ./mem${MEMBER_CH}
    MEMBER=$(( $MEMBER + 1 ))
  done

  rm -f *ensmean

fi

set +x
echo DATA PULL FOR $bundle DONE

exit 0
