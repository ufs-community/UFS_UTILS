#! /usr/bin/env bash
#
# This build script is only used on officially supported machines.  All other
# users should set module files as needed, and build directly with CMake.
#
#
# George Gayno

set -eux

# Get the root of the cloned ufs-utils directory
if [[ $(uname -s) == Darwin ]]; then
  readonly DIR_ROOT=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly DIR_ROOT=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi

source "${DIR_ROOT}/sorc/machine-setup.sh"

# User Options
target=${target:-"NULL"}
compiler=${compiler:-"intel"}

if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
  unset -f module
  set +x
  source "${DIR_ROOT}/modulefiles/build.${target}" > /dev/null
  set -x
else
  set +x
  module use "${DIR_ROOT}/modulefiles"
  module load "build.$target.$compiler" > /dev/null
  module list
  set -x
fi

# Ensure the submodules have been initialized.

if [[ ! -d "${DIR_ROOT}/ccpp-physics/physics" ]]; then
  cd "${DIR_ROOT}"
  git submodule init
  git submodule update
fi

# Collect BUILD Options
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=${BUILD_TYPE:-Release}"
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug"

# Install options; destination for built executables, libraries, CMake Package config
CMAKE_FLAGS+=" -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX:-${DIR_ROOT}} -DCMAKE_INSTALL_BINDIR=${INSTALL_BINDIR:-exec}"

# Testing options
# The unit test data download is part of the build system. Not all machines can
# access the EMC ftp site, so turn off the build (-DBUILD_TESTING=OFF) of the units tests accordingly.
# Those with access to the EMC ftp site are: Orion and Hera.
CMAKE_FLAGS+=" -DBUILD_TESTING=${BUILD_TESTING:-OFF}"

# Allow users of this script to provide CMake options e.g. -DGFS=ON|OFF to build GFS specific utilities only
CMAKE_OPTS=${CMAKE_OPTS:-}

# Re-use or create a new BUILD_DIR (Default: create new BUILD_DIR)
BUILD_DIR=${BUILD_DIR:-"${DIR_ROOT}/build"}
[[ ${BUILD_CLEAN:-"YES"} =~ [yYtT] ]] && rm -rf "$BUILD_DIR"
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"

cmake ${CMAKE_FLAGS} ${CMAKE_OPTS} "${DIR_ROOT}"

make -j "${BUILD_JOBS:-8}" VERBOSE="${BUILD_VERBOSE:-}"
make install

#ctest
#ctest -I 4,5

exit 0
