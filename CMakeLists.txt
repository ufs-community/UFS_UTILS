# This is the main CMake file for NCEPLIBS-ip.
#
# George Gayno
cmake_minimum_required(VERSION 3.15)

# Get the version from the VERSION file.
file(STRINGS "VERSION" pVersion)

project(
  ufs_util
  VERSION ${pVersion}
  LANGUAGES C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# User options.
option(OPENMP "use OpenMP threading" ON)
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)

# Utilities to be built (Default: ALL)
option(ICEBLEND "Enable building emcsfc_ice_blend.fd" ON)
option(SNOW2MDL "Enable building emcsfc_snow2mdl.fd" ON)
option(GCYCLE "Enable building global_cycle.fd" ON)
option(FRENCTOOLS "Enable building fre-nctools.fd" ON)
option(GRIDTOOLS "Enable building grid_tools.fd" ON)
option(CHGRES "Enable building chgres_cube.fd" ON)
option(CHGRES_ALL "Build chgres with all input data options." OFF)
option(OROG_MASK_TOOLS "Enable building orog_mask_tools.fd" ON)
# OROG_MASK_TOOLS must be ON for OROG_NETCDF_TOOLS to build.
option(OROG_NETCDF_TOOLS "Enable building orog_netcdf_tools.fd" OFF)
option(SFC_CLIMO_GEN "Enable building sfc_climo_gen.fd" ON)
option(VCOORD_GEN "Enable building vcoord_gen.fd" ON)
option(FVCOMTOOLS "Enable building fvcom_tools.fd" ON)
option(GBLEVENTS "Enable building gblevents.fd" ON)
option(CPLD_GRIDGEN "Enable building cpld_gridgen.fd" ON)
option(WEIGHT_GEN "Enable building weight_gen.fd" ON)
option(OCEAN_MERGE "Enable building ocean_merge.fd" ON)

# Option to build application specific utilities
option(GFS "Enable building GFS-only utilities" OFF)

# When building the GFS, the following need not be built
if(GFS)
  message(STATUS "Building utilities specific to the GFS")
  set(FRENCTOOLS OFF CACHE BOOL "Disable building fre-nctools.fd" FORCE)
  set(GRIDTOOLS OFF CACHE BOOL "Disable building grid_tools.fd" FORCE)
  set(OROG_MASK_TOOLS OFF CACHE BOOL "Disable building orog_mask_tools.fd" FORCE)
  set(SFC_CLIMO_GEN OFF CACHE BOOL "Disable building sfc_climo_gen.fd" FORCE)
  set(VCOORD_GEN OFF CACHE BOOL "Disable building vcoord_gen.fd" FORCE)
  set(FVCOMTOOLS OFF CACHE BOOL "Disable building fvcom_tools.fd" FORCE)
  set(GBLEVENTS OFF CACHE BOOL "Disable building gblevents.fd" FORCE)
  set(CPLD_GRIDGEN OFF CACHE BOOL "Disable building cpld_gridgen.fd" FORCE)
  set(WEIGHT_GEN OFF CACHE BOOL "Disable building weight_gen.fd" FORCE)
  set(OCEAN_MERGE OFF CACHE BOOL "Disable building ocean_merge.fd" FORCE)
endif()

SET(TEST_FILE_DIR "." CACHE STRING "Check this directory for test files before using FTP.")

# Set the build type.
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# Set compiler flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -ftrapuv")
  if(APPLE)
    # The linker on macOS does not include `common symbols` by default.
    # Passing the -c flag includes them and fixes an error with undefined symbols.
    set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
    set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fcheck=all")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
endif()

# Find packages.
find_package(NetCDF 4.3.3 REQUIRED C Fortran)
find_package(MPI REQUIRED C Fortran)
find_package(ESMF 8.0.0 REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_package(sfcio 1.4.0 REQUIRED)
find_package(w3emc 2.9.0 REQUIRED)
find_package(bacio 2.4.0 REQUIRED)
if(CHGRES_ALL OR GBLEVENTS)
  find_package(nemsio 2.5.0 REQUIRED)
endif()
find_package(sigio 2.3.0 REQUIRED)
find_package(ip 3.3.3 REQUIRED)
if(ip_VERSION LESS 5.0)
  find_package(sp 2.3.3 REQUIRED)
endif()
find_package(g2 3.4.3 REQUIRED)
find_package(sigio 2.3.0 REQUIRED)

# If doxygen documentation we enabled, build it. This must come before
# adding the source code directories; the main documentation build
# must happen before any of the utility document builds.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  set(abs_top_srcdir "${CMAKE_SOURCE_DIR}")
  add_subdirectory(docs)
endif()

add_subdirectory(sorc)

# Run unit tests.
include(CTest)
if(BUILD_TESTING)
  message(STATUS "Set BUILD_TESTING to YES and build unit testing package under tests")
  add_subdirectory(tests)
endif()
