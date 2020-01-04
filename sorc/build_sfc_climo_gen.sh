#! /usr/bin/env bash
set -eux

if [[ $target == "linux.gnu" || $target == "linux.intel" ]]; then
 unset -f module
else
 source ./machine-setup.sh > /dev/null 2>&1
fi

cwd=`pwd`

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  mkdir ../exec
fi

USE_PREINST_LIBS=${USE_PREINST_LIBS:-"true"}
if [ $USE_PREINST_LIBS = true ]; then
  export MOD_PATH
  source ../modulefiles/modulefile.sfc_climo_gen.${target}             > /dev/null 2>&1
else
  export MOD_PATH=${cwd}/lib/modulefiles
  if [ $target = wcoss_cray ]; then
    source ../modulefiles/modulefile.sfc_climo_gen.${target}_userlib   > /dev/null 2>&1
  else
    source ../modulefiles/modulefile.sfc_climo_gen.${target}           > /dev/null 2>&1
  fi
fi

cd ./sfc_climo_gen.fd

make clean
make
make install

echo; echo DONE BUILDING sfc_climo_gen
