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
find_package(MPI REQUIRED )
find_package(ESMF MODULE REQUIRED)


set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

set(fortran_src
    driver.F90
    interp.F90
    model_grid.F90
    output.f90
    program_setup.f90
    search.f90
    source_grid.F90
    utils.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise -r8 -i4 -convert big_endian")
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

set(exe_name sfc_climo_gen.exe)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE ${ESMF_INC} ${NETCDF_INCLUDES})
target_link_libraries(
  ${exe_name}
  esmf
  MPI::MPI_Fortran
  ${NETCDF_LIBRARIES})

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
