cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(
  ufs_util
  VERSION ${pVersion}
  LANGUAGES Fortran)

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

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

option(OPENMP "use OpenMP threading" OFF)

if(NOT NETCDF)
  if(NOT DEFINED ENV{NETCDF})
    message(FATAL_ERROR "fatal error netcdf not set")
  else()
    list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF})
    set(NETCDF $ENV{NETCDF})
  endif()
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)

if(NOT ESMFMKFILE)
  if(NOT DEFINED ENV{ESMFMKFILE})
    message(FATAL_ERROR "fatal error esfm not set")
  else()
    list(APPEND CMAKE_PREFIX_PATH $ENV{ESMFMKFILE})
    set(ESMFMKFILE $ENV{ESMFMKFILE})
  endif()
endif()

find_package(WGRIB2 REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT TARGET sfcio_4)
  set(SFCIO_LIB4 $ENV{SFCIO_LIB4})
endif()

if(NOT TARGET w3nco_d)
  set(W3NCO_LIBd $ENV{W3NCO_LIBd})
endif()

if(NOT TARGET bacio_4)
  set(BACIO_LIB4 $ENV{BACIO_LIB4})
endif()

if(NOT TARGET nemsio)
  set(NEMSIO_LIB $ENV{NEMSIO_LIB})
endif()

if(NOT TARGET sigio_4)
  set(SIGIO_LIB4 $ENV{SIGIO_LIB4})
endif()

if(NOT TARGET sp_d)
  set(SP_LIBd $ENV{SP_LIBd})
endif()

add_subdirectory(sorc)

install(FILES parm/varmap_tables/GFSphys_var_map.txt DESTINATION share)
