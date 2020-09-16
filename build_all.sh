#! /usr/bin/env bash
set -eux

target=${target:-"NULL"}

if [[ "$target" == "linux.gnu" || "$target" == "linux.intel" ]]; then
 unset -f module
else
 set +x
 source ./sorc/machine-setup.sh > /dev/null 2>&1
 set -x
fi

export MOD_PATH
set +x
module use ./modulefiles/build.$target #> /dev/null 2>&1
module list
set -x

#
# --- Build all programs.
#

rm -fr ./build
mkdir ./build
cd ./build

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON"

if [[ "$target" != "wcoss_cray" && "$target" != "odin" ]]; then
  CMAKE_FLAGS+=" -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc"
fi

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

exit
