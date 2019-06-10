#!/bin/sh
set -eux

source ./machine-setup.sh > /dev/null 2>&1
cwd=`pwd`

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  mkdir ../exec
fi

USE_PREINST_LIBS=${USE_PREINST_LIBS:-"true"}
if [ $USE_PREINST_LIBS = true ]; then
  export MOD_PATH=/scratch3/NCEPDEV/nwprod/lib/modulefiles
  source ../modulefiles/modulefile.gettrk.$target             > /dev/null 2>&1
else
  export MOD_PATH=${cwd}/lib/modulefiles
  if [ $target = wcoss_cray ]; then
    source ../modulefiles/modulefile.gettrk.${target}_userlib > /dev/null 2>&1
  else
    source ../modulefiles/modulefile.gettrk.$target           > /dev/null 2>&1
  fi
fi

module list

cd gettrk.fd

make clean
make -f makefile
make install
make clean

cd ../

exit
