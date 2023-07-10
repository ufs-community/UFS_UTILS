#!/bin/bash

##
## Original code to merge by Shan Sun (shan.sun@noaa.gov )
## 
## Modified by Rahul 
## Modified by Sanath to be included into one grid generation process May 2023
## 
##



UFS_DIR=$home_dir

## get the required inputs from command line


nargv=$#

if [ $nargv -eq 2 ];  then  # we have the right number of variables

  res=$1 
  ocn=$2

else
  echo "Number of arguments must be 2 for fv3gfs_ocean_merge"
  echo "Usage fv3gfs_ocean_merge.sh res ocn"
  exit

fi

## change directories to where the exec resides
cd $UFS_DIR/exec


## create the input.nml file with the right info

cat << EOF > input.nml

&mask_nml
 ocean_mask_dir="$ocean_mask_dir/${ocn}/"
 ocnres="mx${ocn}"
 lake_mask_dir="$out_dir/"
 atmres="C${res}"
 out_dir="$TEMP_DIR/ocean_merged/C${res}.mx${ocn}/"
/  
EOF


if [[ ! -f  $UFS_DIR/exec/ocean_merge ]]; then

    error "ocean_merge exe file is not found in $UFS_DIR/exec/ . Try -b to build or -h for help."
else
    echo "ocean_merge  exe file is found in  $UFS_DIR/exec/"
fi

results_dir=$TEMP_DIR/ocean_merged/C${res}.mx${ocn}

mkdir -p $results_dir


./ocean_merge


rm -f input.nml

err=$?
  if [ $err != 0 ]; then
    echo error in fv3gfs_ocean_merge.sh
    exit $err
  fi




