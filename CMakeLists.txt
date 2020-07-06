cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(
  ufs_util
  VERSION ${pVersion}
  LANGUAGES C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|Clang|AppleClang)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_C_COMPILER_ID}")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-g -traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " ")
endif()

find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Jasper REQUIRED)
find_package(NetCDF REQUIRED C Fortran)
find_package(MPI REQUIRED )
find_package(ESMF MODULE REQUIRED)
find_package(WGRIB2 REQUIRED)

option(OPENMP "use OpenMP threading" ON)
if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT TARGET gfsio_4)
  find_package(gfsio REQUIRED)
endif()

if(NOT TARGET sfcio_4)
  find_package(sfcio REQUIRED)
endif()

if(NOT TARGET w3nco_d)
  find_package(w3nco REQUIRED)
endif()

if(NOT TARGET landsfcutil_d)
  find_package(landsfcutil REQUIRED)
endif()

if(NOT TARGET bacio_4)
  find_package(bacio REQUIRED)
endif()

if(NOT TARGET nemsio)
  find_package(nemsio REQUIRED)
endif()

if(NOT TARGET nemsiogfs)
  find_package(nemsiogfs REQUIRED)
endif()

if(NOT TARGET sigio_4)
  find_package(sigio REQUIRED)
endif()

if(NOT TARGET sp_d)
  find_package(sp REQUIRED)
endif()

if(NOT TARGET ip_d)
  find_package(ip REQUIRED)
endif()

if(NOT TARGET w3emc_d)
  find_package(w3emc REQUIRED)
endif()

if(NOT TARGET g2_d)
  find_package(g2 REQUIRED)
endif()

# EMC requires executables in ./exec
set(exec_dir bin)
if(EMC_EXEC_DIR)
  set(exec_dir exec)
endif()

add_subdirectory(sorc)
