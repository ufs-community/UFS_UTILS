#! /usr/bin/env bash
#
# This build script is only used on officially supported machines.  All other
# users should set module files as needed, and build directly with CMake.
#
#
# George Gayno

set -eux

target=${target:-"NULL"}
compiler=${compiler:-"intel"}
PW_CSP=${PW_CSP:-} # TODO: This is an implementation from EPIC and consistent with the UFS WM build system.
export MOD_PATH

if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
 unset -f module
 set +x
 source ./modulefiles/build.$target > /dev/null
 set -x
#TODO: This will need to be revisited once the EPIC supported-stacks come online.
elif [[ "${PW_CSP}" == "aws" ]]; then
  module use /contrib/global-workflow/spack-stack/envs/ufsutils/install/modulefiles/Core
  module load stack-intel
  module load stack-intel-oneapi-mpi
  module use ./modulefiles
  module load build.noaacloud.intel
  module list
  set -x
else
 set +x
 source ./sorc/machine-setup.sh
 module use ./modulefiles
 module load build.$target.$compiler > /dev/null
 module list
 set -x
fi

# Ensure the submodules have been initialized.

if [[ ! -d ./ccpp-physics/physics ]]; then
  git submodule init
  git submodule update
fi

# The unit test data download is part of the build system. Not all machines can
# access the EMC ftp site, so turn off the build (-DBUILD_TESTING=OFF) of the units tests accordingly.
# Those with access to the EMC ftp site are: Orion and Hera.

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DCMAKE_INSTALL_BINDIR=exec -DBUILD_TESTING=OFF"

# Allow users of this script to provide CMake options e.g. -DGFS=ON|OFF to build GFS specific utilities only
CMAKE_OPTS=${CMAKE_OPTS:-}

rm -fr ./build
mkdir ./build && cd ./build

cmake .. ${CMAKE_FLAGS} ${CMAKE_OPTS}

make -j ${BUILD_JOBS:-8} VERBOSE=${BUILD_VERBOSE:-}
make install

#ctest
#ctest -I 4,5

exit
