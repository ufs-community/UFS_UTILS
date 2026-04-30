# FindPIO.cmake
#
# Copyright UCAR 2020
# Copyright NOAA/NWS/NCEP/EMC 2020
#
# Find PIO: A high-level Parallel I/O Library for structured grid applications
# https://github.com/NCAR/ParallelIO
#
# Components available for query:
#  C - Has C support
#  Fortran - Has Fortran support
#  SHARED - Has shared targets
#  STATIC - Has static targets
#
# Variables provided:
#  PIO_FOUND - True if PIO was found
#  PIO_C_FOUND - True if PIO C support was found
#  PIO_Fortran_FOUND - True if PIO Fortran support was found
#  PIO_VERSION - Version of installed PIO
#
# Targets provided:
#  PIO::PIO_C - C interface target aliased to SHARED|STATIC as requested or to shared libraries if available else static libraries
#  PIO::PIO_Fortran - Fortran interface target aliases to SHARED|STATIC as requested or to shared libraries if available else static libraries
#
# To control finding of this package, set PIO_ROOT environment variable to the full path to the prefix
# under which PIO was installed (e.g., /usr/local)

set( _search_components )
set( _search_library_type )
foreach( _comp ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS} )
    if( _comp MATCHES "^(STATIC|SHARED)$" )
        list( APPEND _search_library_type ${_comp} )
    else()
        list( APPEND _search_components ${_comp} )
    endif()
endforeach()
set( ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS ${_search_components} )

# If no COMPONENTS are requested, seach both C and Fortran
if( NOT _search_components )
    list( APPEND _search_components C Fortran )
endif()

# Ensure there is only one type of library being requested
if( _search_library_type )
    list( LENGTH _search_library_type _len)
    if( _len GREATER 1 )
        message(FATAL_ERROR "User requesting both STATIC and SHARED is not permissible")
    endif()
    unset(_len)
endif()

## Find libraries and paths, and determine found components
find_path(PIO_INCLUDE_DIR NAMES pio.h HINTS "${PIO_PREFIX}" PATH_SUFFIXES include include/pio)
if(PIO_INCLUDE_DIR)
    string(REGEX REPLACE "/include(/.+)?" "" PIO_PREFIX ${PIO_INCLUDE_DIR})
    set(PIO_PREFIX ${PIO_PREFIX} CACHE STRING "")
    find_path(PIO_MODULE_DIR NAMES pio.mod PATHS "${PIO_PREFIX}"
              PATH_SUFFIXES include include/pio lib/pio/module module module/pio NO_DEFAULT_PATH)
    if(APPLE)
      set(_SHARED_LIB_EXT dylib)
    else()
      set(_SHARED_LIB_EXT so)
    endif()
    find_library(PIO_C_STATIC_LIB libpioc.a PATHS "${PIO_PREFIX}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
    find_library(PIO_C_SHARED_LIB libpioc.${_SHARED_LIB_EXT} PATHS "${PIO_PREFIX}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
    find_library(PIO_Fortran_STATIC_LIB libpiof.a PATHS "${PIO_PREFIX}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
    find_library(PIO_Fortran_SHARED_LIB libpiof.${_SHARED_LIB_EXT} PATHS "${PIO_PREFIX}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
    unset(_SHARED_LIB_EXT)
    #Check for Fortran components
    if(PIO_MODULE_DIR)
        if(PIO_Fortran_STATIC_LIB)
            set(PIO_Fortran_STATIC_FOUND 1)
        endif()
        if(PIO_Fortran_SHARED_LIB)
            set(PIO_Fortran_SHARED_FOUND 1)
        endif()
    endif()

    #Set version using installed .settings file
    file(READ "${PIO_PREFIX}/lib/libpio.settings" ver)
    string(REGEX MATCH "PIO Version:[ \t]*([0-9]*).([0-9]*).([0-9]*)" _ ${ver})
    set(PIO_VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")

    #Check for C components
    if(PIO_C_STATIC_LIB)
        set(PIO_C_STATIC_FOUND 1)
    endif()
    if(PIO_C_SHARED_LIB)
        set(PIO_C_SHARED_FOUND 1)
    endif()
endif()
## Debugging output
message(DEBUG "[FindPIO] PIO_INCLUDE_DIR: ${PIO_INCLUDE_DIR}")
message(DEBUG "[FindPIO] PIO_PREFIX: ${PIO_PREFIX}")
message(DEBUG "[FindPIO] PIO_MODULE_DIR: ${PIO_MODULE_DIR}")
message(DEBUG "[FindPIO] PIO_C_STATIC_LIB: ${PIO_C_STATIC_LIB}")
message(DEBUG "[FindPIO] PIO_C_SHARED_LIB: ${PIO_C_SHARED_LIB}")
message(DEBUG "[FindPIO] PIO_C_SHARED_FOUND: ${PIO_C_SHARED_FOUND}")
message(DEBUG "[FindPIO] PIO_C_STATIC_FOUND: ${PIO_C_STATIC_FOUND}")
message(DEBUG "[FindPIO] PIO_Fortran_STATIC_LIB: ${PIO_Fortran_STATIC_LIB}")
message(DEBUG "[FindPIO] PIO_Fortran_SHARED_LIB: ${PIO_Fortran_SHARED_LIB}")
message(DEBUG "[FindPIO] PIO_Fortran_SHARED_FOUND: ${PIO_Fortran_SHARED_FOUND}")
message(DEBUG "[FindPIO] PIO_Fortran_STATIC_FOUND: ${PIO_Fortran_STATIC_FOUND}")

## Create targets
set(_new_components)
# PIO_C_STATIC imported interface target
if(PIO_C_STATIC_FOUND AND NOT TARGET PIO_C_STATIC)
    add_library(PIO_C_STATIC INTERFACE IMPORTED)
    set_target_properties(PIO_C_STATIC PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${PIO_INCLUDE_DIR}
                          INTERFACE_LINK_LIBRARIES ${PIO_C_STATIC_LIB}
                          IMPORTED_GLOBAL True )
    set(_new_components 1)
endif()
# PIO_C_SHARED imported interface target
if(PIO_C_SHARED_FOUND AND NOT TARGET PIO_C_SHARED)
    add_library(PIO_C_SHARED INTERFACE IMPORTED)
    set_target_properties(PIO_C_SHARED PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${PIO_INCLUDE_DIR}
                          INTERFACE_LINK_LIBRARIES ${PIO_C_SHARED_LIB}
                          IMPORTED_GLOBAL True )
    set(_new_components 1)
endif()
# PIO_Fortran_STATIC imported interface target
if(PIO_Fortran_STATIC_FOUND AND NOT TARGET PIO_Fortran_STATIC)
    add_library(PIO_Fortran_STATIC INTERFACE IMPORTED)
    set_target_properties(PIO_Fortran_STATIC PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${PIO_INCLUDE_DIR}
                          INTERFACE_LINK_LIBRARIES ${PIO_Fortran_STATIC_LIB}
                          IMPORTED_GLOBAL True )
    if(PIO_MODULE_DIR AND NOT PIO_MODULE_DIR STREQUAL PIO_INCLUDE_DIR )
        set_property(TARGET PIO_Fortran_STATIC APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PIO_MODULE_DIR})
    endif()
    set(_new_components 1)
    target_link_libraries(PIO_Fortran_STATIC INTERFACE PIO_C_STATIC)
endif()
# PIO_Fortran_SHARED imported interface target
if(PIO_Fortran_SHARED_FOUND AND NOT TARGET PIO_Fortran_SHARED)
    add_library(PIO_Fortran_SHARED INTERFACE IMPORTED)
    set_target_properties(PIO_Fortran_SHARED PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${PIO_INCLUDE_DIR}
                          INTERFACE_LINK_LIBRARIES ${PIO_Fortran_SHARED_LIB}
                          IMPORTED_GLOBAL True )
    if(PIO_MODULE_DIR AND NOT PIO_MODULE_DIR STREQUAL PIO_INCLUDE_DIR )
        set_property(TARGET PIO_Fortran_SHARED APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PIO_MODULE_DIR})
    endif()
    target_link_libraries(PIO_Fortran_SHARED INTERFACE PIO_C_SHARED)
    set(_new_components 1)
endif()

if( _search_library_type MATCHES "^(SHARED)$" )
  if( TARGET PIO_C_SHARED )
      add_library(PIO::PIO_C ALIAS PIO_C_SHARED)
      set(PIO_C_FOUND 1)
  endif()
  if( TARGET PIO_Fortran_SHARED )
      add_library(PIO::PIO_Fortran ALIAS PIO_Fortran_SHARED)
      set(PIO_Fortran_FOUND 1)
  endif()
elseif( _search_library_type MATCHES "^(STATIC)$" )
  if( TARGET PIO_C_STATIC )
      add_library(PIO::PIO_C ALIAS PIO_C_STATIC)
      set(PIO_C_FOUND 1)
  endif()
  if( TARGET PIO_Fortran_STATIC )
      add_library(PIO::PIO_Fortran ALIAS PIO_Fortran_STATIC)
      set(PIO_Fortran_FOUND 1)
  endif()
else()
  if( TARGET PIO_C_SHARED )
        add_library(PIO::PIO_C ALIAS PIO_C_SHARED)
        set(PIO_C_FOUND 1)
  elseif( TARGET PIO_C_STATIC )
        add_library(PIO::PIO_C ALIAS PIO_C_STATIC)
        set(PIO_C_FOUND 1)
  endif()
  if( TARGET PIO_Fortran_SHARED )
      add_library(PIO::PIO_Fortran ALIAS PIO_Fortran_SHARED)
      set(PIO_Fortran_FOUND 1)
  elseif( TARGET PIO_Fortran_STATIC )
      add_library(PIO::PIO_Fortran ALIAS PIO_Fortran_STATIC)
      set(PIO_Fortran_FOUND 1)
  endif()
endif()

## Check package has been found correctly
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  PIO
  REQUIRED_VARS
    PIO_PREFIX
    PIO_INCLUDE_DIR
  VERSION_VAR PIO_VERSION
  HANDLE_COMPONENTS
)
message(DEBUG "[FindPIO] PIO_FOUND: ${PIO_FOUND}")

## Print status
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND AND NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY AND _new_components)
    message( STATUS "Find${CMAKE_FIND_PACKAGE_NAME}:" )
    message( STATUS "  - ${CMAKE_FIND_PACKAGE_NAME}_PREFIX [${${CMAKE_FIND_PACKAGE_NAME}_PREFIX}]")
    message( STATUS "  - ${CMAKE_FIND_PACKAGE_NAME}_VERSION: [${${CMAKE_FIND_PACKAGE_NAME}_VERSION}]")
    set(_found_comps)
    foreach( _comp ${_search_components} )
        if( ${CMAKE_FIND_PACKAGE_NAME}_${_comp}_FOUND )
            list(APPEND _found_comps ${_comp})
        endif()
    endforeach()
    message( STATUS "  - ${CMAKE_FIND_PACKAGE_NAME} Components Found: ${_found_comps}")
    unset(_found_comps)
endif()
unset(_new_components)
unset(_search_components)
unset(_search_library_type)
unset(_library_type)
