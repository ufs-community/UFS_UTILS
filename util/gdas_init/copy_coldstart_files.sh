#!/bin/bash

copy_data()
{

set -x

mkdir -p $SAVEDIR
cp gfs_ctrl.nc $SAVEDIR

for tile in 'tile1' 'tile2' 'tile3' 'tile4' 'tile5' 'tile6'
do
  cp out.atm.${tile}.nc  ${SAVEDIR}/gfs_data.${tile}.nc
  cp out.sfc.${tile}.nc  ${SAVEDIR}/sfc_data.${tile}.nc
done
}

set -x

MEMBER=$1
OUTDIR=$2
yy=$3
mm=$4
dd=$5
hh=$6
INPUT_DATA_DIR=$7

if [ ${MEMBER} == 'hires' ]; then
  MEMBER='gdas'
fi

echo 'in new script'
echo $MEMBER $OUTDIR

if [ ${MEMBER} == 'gdas' ] || [ ${MEMBER} == 'gfs' ]; then
  SUBDIR=$OUTDIR/${MEMBER}.${yy}${mm}${dd}/${hh}
  rm -fr $SUBDIR
  SAVEDIR=$SUBDIR/atmos/INPUT
  copy_data
  if [ ${MEMBER} == 'gdas' ]; then
    cp ${INPUT_DATA_DIR}/*abias* $SAVEDIR/..
    cp ${INPUT_DATA_DIR}/*radstat $SAVEDIR/..
  fi
  touch $SAVEDIR/../${MEMBER}.t${hh}z.loginc.txt
elif [ ${MEMBER} == 'enkf' ]; then  # v16 retro data only.
  MEMBER=1
  while [ $MEMBER -le 80 ]; do
    if [ $MEMBER -lt 10 ]; then
      MEMBER_CH="00${MEMBER}"
    else
      MEMBER_CH="0${MEMBER}"
    fi
    SUBDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/mem${MEMBER_CH}
    rm -fr $SUBDIR
    SAVEDIR=$SUBDIR/atmos/INPUT
    copy_data
    touch $SAVEDIR/../enkfgdas.t${hh}z.loginc.txt
    MEMBER=$(( $MEMBER + 1 ))
  done
else
  SUBDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/mem${MEMBER}
  rm -fr $SUBDIR
  SAVEDIR=$SUBDIR/atmos/INPUT
  copy_data
  touch $SAVEDIR/../enkfgdas.t${hh}z.loginc.txt
fi

exit 0
