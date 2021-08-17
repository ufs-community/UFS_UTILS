#!/bin/bash

#----------------------------------------------------------------------
# Retrieve gfs v16 retrospective parallel data from hpss.
#----------------------------------------------------------------------

set -x

MEMBER=$1

cd $EXTRACT_DIR

if [ "$MEMBER" = "gfs" ]; then

  if [ ${yy}${mm}${dd}${hh} -lt 2019050106 ]; then
    set +x
    echo NO DATA FOR ${yy}${mm}${dd}${hh}
    exit 2
  elif [ ${yy}${mm}${dd}${hh} -lt 2019060106 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro0e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -lt 2019090100 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro1e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -lt 2019111000 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro2e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -le 2020122200 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -le 2021032506 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2n/${yy}${mm}${dd}${hh}
  else
    set +x
    echo NO DATA FOR ${yy}${mm}${dd}${hh}
    exit 3
  fi

  file=gfs_netcdfa.tar

  touch ./list.hires3
  htar -tvf  $directory/$file > ./list.hires1
  grep sfcanl ./list.hires1 > ./list.hires2
  grep atmanl ./list.hires1 >> ./list.hires2
  while read -r line
  do
    echo ${line##*' '} >> ./list.hires3
  done < "./list.hires2"

  htar -xvf $directory/$file -L ./list.hires3
  rc=$?
  [ $rc != 0 ] && exit $rc
  
  rm -f ./list.hires?

else

  date10_m6=`$NDATE -6 $yy$mm$dd$hh`

  echo $date10_m6
  yy_m6=$(echo $date10_m6 | cut -c1-4)
  mm_m6=$(echo $date10_m6 | cut -c5-6)
  dd_m6=$(echo $date10_m6 | cut -c7-8)
  hh_m6=$(echo $date10_m6 | cut -c9-10)

  if [ $date10_m6 -lt 2019050100 ]; then
   set +x
   echo NO DATA FOR $date10_m6
   exit 2
  elif [ $date10_m6 -le 2019060100 ]; then
   directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro0e/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
  elif [ $date10_m6 -lt 2019090100 ]; then
   directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro1e/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
  elif [ $date10_m6 -lt 2019101700 ]; then
   directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro2e/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
  elif [ $date10_m6 -lt 2020122212 ]; then
   directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
  elif [ $date10_m6 -le 2021032500 ]; then
   directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2n/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
  else
    set +x
    echo NO DATA FOR $date10_m6
    exit 3
  fi

#----------------------------------------------------------------------
# Pull restart files.
#----------------------------------------------------------------------

  file=gdas_restartb.tar

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
# Pull abias and radstat files.
#----------------------------------------------------------------------

  rm -f ./list.hires*

  if [ ${yy}${mm}${dd}${hh} -lt 2019050106 ]; then
    set +x
    echo NO DATA FOR ${yy}${mm}${dd}${hh}
    exit 2
  elif [ ${yy}${mm}${dd}${hh} -lt 2019060106 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro0e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -lt 2019090106 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro1e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -lt 2019101706 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro2e/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -lt 2020122218 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2/${yy}${mm}${dd}${hh}
  elif [ ${yy}${mm}${dd}${hh} -le 2021032506 ]; then
    directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16rt2n/${yy}${mm}${dd}${hh}
  else
    set +x
    echo NO DATA FOR ${yy}${mm}${dd}${hh}
    exit 3
  fi

  file=gdas_restarta.tar

  touch ./list.hires3
  htar -tvf  $directory/$file > ./list.hires1
  grep abias ./list.hires1 > ./list.hires2
  grep radstat ./list.hires1 >> ./list.hires2
  while read -r line
  do
    echo ${line##*' '} >> ./list.hires3
  done < "./list.hires2"

  htar -xvf $directory/$file -L ./list.hires3
  rc=$?
  [ $rc != 0 ] && exit $rc

fi # is this gdas or gfs CDUMP?

set +x
echo DATA PULL FOR v16 DONE

exit 0
