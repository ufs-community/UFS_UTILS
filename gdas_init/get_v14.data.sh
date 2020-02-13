#!/bin/bash

#-----------------------------------------------------------
# Retrieve gfs v14 data from hpss.
#
# v14 data starts July 19, 2017 at 12z
#-----------------------------------------------------------

bundle=$1

set -x

if [ $bundle = 'hires' ]; then

  mkdir -p $EXTRACT_DIR/gdas.${yy}${mm}${dd}/${hh}
  cd $EXTRACT_DIR/gdas.${yy}${mm}${dd}/${hh}

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  file=gpfs_hps_nco_ops_com_gfs_prod_gdas.${yy}${mm}${dd}${hh}.tar

  htar -xvf $directory/$file ./gdas.t${hh}z.radstat
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.abias_air
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.abias
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.abias_pc
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.atmanl.nemsio
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.nstanl.nemsio
  rc=$?
  [ $rc != 0 ] && exit $rc
  htar -xvf $directory/$file ./gdas.t${hh}z.sfcanl.nemsio
  rc=$?
  [ $rc != 0 ] && exit $rc

  set +x
  echo DATA PULL FOR $bundle DONE

  exit 0

elif [ $bundle = 'enkf' ]; then

#----------------------------------------------------------------------
# Get the enkf tiled restart files for all members.
#----------------------------------------------------------------------

  mkdir -p $EXTRACT_DIR/enkf.${yy}${mm}${dd}/${hh}
  cd $EXTRACT_DIR/enkf.${yy}${mm}${dd}/${hh}

  directory=/NCEPPROD/hpssprod/runhistory/rh${yy}/${yy}${mm}/${yy}${mm}${dd}
  file=gpfs_hps_nco_ops_com_gfs_prod_enkf.${yy}${mm}${dd}_${hh}.anl.tar

  htar -xvf $directory/$file
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
    mv *.mem${MEMBER_CH}* ./mem${MEMBER_CH}
    MEMBER=$(( $MEMBER + 1 ))
  done

  rm -f gdas.*

fi

set +x
echo DATA PULL FOR $bundle DONE

exit 0
