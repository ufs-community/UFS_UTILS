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

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
  message(STATUS "Set BUILD_TESTING to YES and build unit testing package under tests")
  set(BUILD_TESTING "YES")
endif()

# Set compiler flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
endif()

# Find packages.
find_package(NetCDF REQUIRED C Fortran)
find_package(MPI REQUIRED)
find_package(ESMF MODULE REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_package(gfsio REQUIRED)
find_package(sfcio REQUIRED)
find_package(w3nco REQUIRED)
find_package(landsfcutil REQUIRED)
find_package(bacio REQUIRED)
find_package(nemsio REQUIRED)
find_package(nemsiogfs REQUIRED)
find_package(sigio REQUIRED)
find_package(sp REQUIRED)
find_package(ip REQUIRED)
find_package(w3emc REQUIRED)
find_package(g2 REQUIRED)
find_package(wgrib2 REQUIRED)

# EMC requires executables in ./exec
set(exec_dir bin)
if(EMC_EXEC_DIR)
  set(exec_dir exec)
endif()

add_subdirectory(sorc)

# Run unit tests.
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(unit-tests)
endif()

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)  
endif()
