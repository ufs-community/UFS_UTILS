cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(
  ufs_util
  VERSION ${pVersion}
  LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/")

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

option(OPENMP "use OpenMP threading" ON)

if(NOT NETCDF)
  if(NOT DEFINED ENV{NETCDF})
    message(FATAL_ERROR "fatal error netcdf not set")
  else()
    list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF})
    set(NETCDF $ENV{NETCDF})
  endif()
  if(DEFINED_ENV{NETCDF_FORTRAN})
    list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF_FORTRAN})
  endif()
endif()
find_package(NetCDF MODULE REQUIRED COMPONENTS C Fortran)
if(NOT TARGET w3nco_d)
  find_package(w3nco REQUIRED)
endif()
if(NOT TARGET sp_d)
  find_package(sp REQUIRED)
endif()
if(NOT TARGET ip_d)
  find_package(ip REQUIRED)
endif()
if(NOT TARGET bacio_4)
  find_package(bacio REQUIRED)
endif()

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()


set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

set(fortran_src
    mtnlm7_oclsm.f
    netcdf_io.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -traceback -r8  -convert big_endian -fp-model precise  -assume byterecl")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-form -ffree-line-length-0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name orog.x)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE ${NETCDF_INCLUDE_DIRS})
target_link_libraries(
  ${exe_name}
  bacio_4
  w3nco_d
  ip_d
  sp_d
  ${NETCDF_LIBRARIES})
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
