# This module looks for environment variables detailing where SIGIO lib is
# If variables are not set, SIGIO will be built from external source
if(DEFINED ENV{SIGIO_LIB} )
  set(SIGIO_LIB $ENV{SIGIO_LIB} CACHE STRING "SIGIO Library Location" )
  set(SIGIO_INC $ENV{SIGIO_INC} CACHE STRING "SIGIO Include Location" )

  set(name "sigio")
  string(TOUPPER ${name} uppercase_name)


  set(lib_name ${name})

  if(EXISTS ${${uppercase_name}_LIB} )
    get_filename_component(lib_dir ${${uppercase_name}_LIB} DIRECTORY)
    find_library(sigio_path NAMES ${lib_name} PATHS ${lib_dir} NO_DEFAULT_PATH)
    if(NOT sigio_path)
      set(lib_name ${name})
      find_library(sigio_path NAMES ${lib_name} PATHS ${lib_dir} NO_DEFAULT_PATH)
    endif()
    add_library(${lib_name} STATIC IMPORTED)
    set_target_properties(${lib_name} PROPERTIES
      IMPORTED_LOCATION ${sigio_path}
      INTERFACE_INCLUDE_DIRECTORIES ${${uppercase_name}_INC${kind}})
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sigio
  REQUIRED_VARS sigio_path)
