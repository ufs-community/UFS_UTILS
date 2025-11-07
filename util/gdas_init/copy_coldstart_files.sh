#!/bin/bash

# Copy files from the working directory to the
# output directory.

copy_data()
{

  set -x

  MEM=$1

  SAVEDIR_MODEL_DATA=$SUBDIR/model/atmos/input
  mkdir -p $SAVEDIR_MODEL_DATA
  cp gfs_ctrl.nc $SAVEDIR_MODEL_DATA

  for tile in 'tile1' 'tile2' 'tile3' 'tile4' 'tile5' 'tile6'
  do
    cp out.atm.${tile}.nc ${SAVEDIR_MODEL_DATA}/gfs_data.${tile}.nc
    cp out.sfc.${tile}.nc ${SAVEDIR_MODEL_DATA}/sfc_data.${tile}.nc
  done

  if [ ${MEM} == 'gdas' ]; then
    SAVEDIR_ANALYSIS=$SUBDIR/analysis/atmos
    mkdir -p $SAVEDIR_ANALYSIS
    for abias_file in ${INPUT_DATA_DIR}/*abias*; do
      base_abias_name=$(basename ${abias_file})
      # Test for v17 style abias file names
      if [[ "${base_abias_name}" == *".txt" ]]; then
        cp ${INPUT_DATA_DIR}/${base_abias_name} $SAVEDIR_ANALYSIS/${base_abias_name}
      else
        cp ${INPUT_DATA_DIR}/${base_abias_name} $SAVEDIR_ANALYSIS/${base_abias_name}.txt
      fi
    done
    for radstat_file in ${INPUT_DATA_DIR}/*radstat; do
      base_radstat_name=$(basename ${radstat_file})
      if [[ "${base_radstat_name}" == *".tar" ]]; then
        cp ${INPUT_DATA_DIR}/${base_radstat_name} $SAVEDIR_ANALYSIS/${base_radstat_name}
      else
        cp ${INPUT_DATA_DIR}/${base_radstat_name} $SAVEDIR_ANALYSIS/${base_radstat_name}.tar
      fi
      group=$(stat -c %G ${INPUT_DATA_DIR}/${base_radstat_name})
      if [[ "${group}" == "rstprod" ]]; then
        chgrp rstprod $SAVEDIR_ANALYSIS/${base_radstat_name}.tar
        chmod 640 $SAVEDIR_ANALYSIS/${base_radstat_name}.tar
      fi
    done
  fi
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

set +x
echo 'COPY DATA TO OUTPUT DIRECTORY'
set -x

if [ ${MEMBER} == 'gdas' ] || [ ${MEMBER} == 'gfs' ]; then
  SUBDIR=$OUTDIR/${MEMBER}.${yy}${mm}${dd}/${hh}
  rm -fr $SUBDIR
  copy_data ${MEMBER}
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
    copy_data ${MEMBER}
    MEMBER=$(( $MEMBER + 1 ))
  done
else
  SUBDIR=$OUTDIR/enkfgdas.${yy}${mm}${dd}/${hh}/mem${MEMBER}
  rm -fr $SUBDIR
  copy_data ${MEMBER}
fi





#------------------------------------------------------------------------------------
# Make the README files with all relevant info to reproduce the outputs
#------------------------------------------------------------------------------------

cd $UFS_DIR

commit_string=$(git log -1 --oneline)
commit_num=$(echo $commit_string | cut -c1-7)

cd ${SAVEDIR_MODEL_DATA}

cat <<EOF > README.TXT
The following parameters were used
creation date=$(date +%Y-%m-%d)
commit_num=$commit_num
yy=$yy
mm=$mm
dd=$dd
hh=$hh
LEVS=$LEVS
CRES_HIRES=$CRES_HIRES
CRES_ENKF=$CRES_ENKF
gfs_ver=$gfs_ver
use_v16retro=$use_v16retro
OUTDIR=$OUTDIR
EXTRACT_DIR=$EXTRACT_DIR
FIX_ORO_INPUT=$FIX_ORO_INPUT

EOF






exit 0
