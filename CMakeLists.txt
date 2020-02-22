cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(
  ufs_util
  VERSION ${pVersion}
  LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

option(OPENMP "use OpenMP threading" OFF)

find_package(
  HDF5
  COMPONENTS HL Fortran
  REQUIRED)
find_package(NetCDF REQUIRED)
find_package(MPI REQUIRED)
find_package(Jasper REQUIRED)
find_package(PNG REQUIRED)
find_package(ESMF MODULE REQUIRED)
find_package(ZLIB REQUIRED)
find_package(WGRIB2 REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT TARGET landsfcutil_4)
  find_package(landsfcutil REQUIRED)
endif()

if(NOT TARGET sfcio_4)
  find_package(sfcio REQUIRED)
endif()

if(NOT TARGET w3nco_d)
  find_package(w3nco REQUIRED)
endif()

if(NOT TARGET bacio_4)
  find_package(bacio REQUIRED)
endif()

if(NOT TARGET nemsio)
  find_package(nemsio REQUIRED)
endif()

if(NOT TARGET sigio_4)
  find_package(sigio REQUIRED)
endif()

if(NOT TARGET sp_d)
  find_package(sp REQUIRED)
endif()

add_subdirectory(sorc)
