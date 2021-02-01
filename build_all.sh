#! /usr/bin/env bash
set -eux

target=${target:-"NULL"}

export MOD_PATH

if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
 unset -f module
 set +x
 source ./modulefiles/build.$target > /dev/null 2>&1 
 set -x
else
 set +x
 source ./sorc/machine-setup.sh
 module use ./modulefiles
 module load build.$target > /dev/null 2>&1
 module list
 set -x
fi

compiler=${compiler:-"intel"}

# --- Build all programs.
#

rm -fr ./build
mkdir ./build
cd ./build

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON"

if [[ "$compiler" == "intel" ]]; then
  if [[ "$target" != "wcoss_cray" && "$target" != "odin" && "$target" != "gaea" ]]; then
    CMAKE_FLAGS+=" -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc"
  fi
fi

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

exit
