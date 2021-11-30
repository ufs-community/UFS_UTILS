#! /usr/bin/env bash
set -eux

source ./machine-setup.sh > /dev/null 2>&1
cwd=`pwd`

USE_PREINST_LIBS=${USE_PREINST_LIBS:-"true"}
if [ $USE_PREINST_LIBS = true ]; then
  module use ../modulefiles/fv3gfs
  module load global_cycle.$target             > /dev/null 2>&1
else
  export MOD_PATH=${cwd}/lib/modulefiles
  source ../modulefiles/fv3gfs/global_cycle.$target           > /dev/null 2>&1
fi
module list

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  mkdir ../exec
fi

cd ${cwd}/global_cycle.fd
./makefile.sh
