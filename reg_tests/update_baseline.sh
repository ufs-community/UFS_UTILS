#!/bin/bash

set -x

HOMEreg=$1
test_name=$2
commit_num=$3

base_dir=$HOMEreg/baseline_data
base_dir_commit=${base_dir}/$test_name.$commit_num

chmod 755 $base_dir

if [ -d $base_dir_commit ];then
  chmod 777 $base_dir_commit
  if [ -d $base_dir_commit/sfc ]; then
    chmod 777 $base_dir_commit/sfc
  fi
  rm -fr $base_dir_commit
fi

mkdir -p $base_dir_commit

for files in *.nc snogrb_model seaice.5min.blend
do
  if [ -f $files ]; then
    cp $files $base_dir_commit
    chmod 444 $base_dir_commit/$files
  fi
done

# The grid_gen tests have a subdirectory for
# the surface climo fields.

if [ -d ./sfc ]; then
  mkdir -p $base_dir_commit/sfc
  cd sfc
  for files in *.nc
  do 
    if [ -f $files ]; then
      cp $files $base_dir_commit/sfc
      chmod 444 $base_dir_commit/sfc/$files
    fi
  done
  chmod 555 $base_dir_commit/sfc
fi

chmod 555 $base_dir_commit
rm -f $base_dir/$test_name
cd $base_dir
ln -fs $test_name.$commit_num $test_name

# move this to driver?
###chmod 555 $base_dir

exit
