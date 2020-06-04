#! /usr/bin/env bash
set -eux

target=${target:-"NULL"}

if [[ $target == "linux.gnu" || $target == "linux.intel" ]]; then
 unset -f module
else
 source ./sorc/machine-setup.sh > /dev/null 2>&1
fi

export MOD_PATH
source ./modulefiles/build.$target             > /dev/null 2>&1

#
# --- Build all programs.
#

rm -fr ./build
mkdir ./build
cd ./build

if [[ $target == "wcoss_cray" ]]; then
  cmake .. -DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON
else
  cmake .. -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON
fi

make -j 8 VERBOSE=1
make install

exit
