#! /usr/bin/env bash
set -eux

target=${target:-"NULL"}

if [[ $target == "linux.gnu" || $target == "linux.intel" ]]; then
 unset -f module
else
 source ./machine-setup.sh > /dev/null 2>&1
fi

cwd=`pwd`

USE_PREINST_LIBS=${USE_PREINST_LIBS:-"true"}
if [ $USE_PREINST_LIBS = true ]; then
  export MOD_PATH
  source ../modulefiles/modules.$target             > /dev/null 2>&1
else
  export MOD_PATH=${cwd}/lib/modulefiles
  if [ $target = wcoss_cray ]; then
    source ../modulefiles/modules.${target}_userlib > /dev/null 2>&1
  else
    source ../modulefiles/modules.$target           > /dev/null 2>&1
  fi
fi

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  mkdir ../exec
fi

#
# --- Chgres part
#
rm -fr ../build
mkdir ../build
cd ../build
if [[ $target == "wcoss_cray" ]]; then
  cmake .. -DCMAKE_INSTALL_PREFIX=../
else
  cmake .. -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DCMAKE_INSTALL_PREFIX=../
fi

make -j 8 VERBOSE=1
make install
exit
