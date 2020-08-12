#!/bin/bash

#----------------------------------------------------------------------
# Retrieve gfs v16 data from hpss.
#
#----------------------------------------------------------------------

set -x

cd $EXTRACT_DIR

date10_m6=`$NDATE -6 $yy$mm$dd$hh`

echo $date10_m6
yy_m6=$(echo $date10_m6 | cut -c1-4)
mm_m6=$(echo $date10_m6 | cut -c5-6)
dd_m6=$(echo $date10_m6 | cut -c7-8)
hh_m6=$(echo $date10_m6 | cut -c9-10)

 directory=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_D/gfsv16/v16retro1e/${yy_m6}${mm_m6}${dd_m6}${hh_m6}
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


set +x
echo DATA PULL FOR v16 DONE

exit 0
