#! /usr/bin/env bash
set -eux

source ./machine-setup.sh > /dev/null 2>&1
cwd=`pwd`

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  mkdir ../exec
fi

module purge
USE_PREINST_LIBS=${USE_PREINST_LIBS:-"true"}
if [ $USE_PREINST_LIBS = true ]; then
  source ../modulefiles/modulefile.regional_grid.${target}             > /dev/null 2>&1
else
  export MOD_PATH=${cwd}/lib/modulefiles
  if [ $target = wcoss_cray ]; then
    source ../modulefiles/modulefile.regional_grid.${target}_userlib   > /dev/null 2>&1
  else
    source ../modulefiles/modulefile.regional_grid.${target}           > /dev/null 2>&1
  fi
fi

module list

export INCS="${NETCDF_INCLUDE}"
export LIBS="${NETCDF_LDFLAGS} ${HDF5_LDFLAGS}"

cd ./regional_grid.fd

make clean
make
make install

echo; echo DONE BUILDING regional_grid
