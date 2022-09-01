#! /usr/bin/env bash
#
# This build script is only used on officially supported machines.  All other
# users should set module files as needed, and build directly with CMake.
#
# George Gayno

set -eux

target=${target:-"NULL"}
compiler=${compiler:-"intel"}
export MOD_PATH

if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
 unset -f module
 set +x
 source ./modulefiles/build.$target > /dev/null
 set -x
else
 set +x
 source ./sorc/machine-setup.sh
 module use ./modulefiles
 module load build.$target.$compiler > /dev/null
 module list
 set -x
fi


# The unit test data download is part of the build system. Not all machines can
# access the EMC ftp site, so turn off the build (-DBUILD_TESTING=OFF) of the units tests accordingly.
# Those with access to the EMC ftp site are: Orion and Hera.

if [[ "$target" == "hera" || "$target" == "orion" || "$target" == "wcoss2" ]]; then
  CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DBUILD_TESTING=OFF"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DBUILD_TESTING=ON"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DENABLE_DOCS=ON -DBUILD_TESTING=ON"
else
  CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DBUILD_TESTING=OFF"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DENABLE_DOCS=ON -DBUILD_TESTING=OFF"
fi

rm -fr ./build
mkdir ./build && cd ./build

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

#make test
#ctest -I 4,5

exit
