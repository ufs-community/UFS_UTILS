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
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv")
  if(APPLE)
    # The linker on macOS does not include `common symbols` by default.
    # Passing the -c flag includes them and fixes an error with undefined symbols.
    set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
    set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#  set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall")
#  set(CMAKE_Fortran_FLAGS_DEBUG "-O1 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check -fno-omit-frame-pointer -fno-optimize-sibling-calls")
  set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
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
find_package(w3nco 2.4.0 REQUIRED)
find_package(bacio 2.4.0 REQUIRED)
find_package(nemsio 2.5.0 REQUIRED)
find_package(sigio 2.3.0 REQUIRED)
find_package(sp 2.3.3 REQUIRED)
find_package(ip 3.3.3 REQUIRED)
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
