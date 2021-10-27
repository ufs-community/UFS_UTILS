# This modulerlooks for environment variables detailing where SFCIO lib is
# If variables are not set, SFCIO will be built from external source
if(DEFINED ENV{SFCIO_LIB} )
  set(SFCIO_LIB $ENV{SFCIO_LIB} CACHE STRING "SFCIO Library Location" )
  set(SFCIO_INC $ENV{SFCIO_INC} CACHE STRING "SFCIO Include Location" )

  set(name "sfcio")
  string(TOUPPER ${name} uppercase_name)

  #string(REGEX MATCH "(v[0-9]+\\.[0-9]+\\.[0-9]+)" _ ${${uppercase_name}_LIB})
  #set(version ${CMAKE_MATCH_1})

  set(lib_name ${name})
  #set(versioned_lib_name ${name})
  if(EXISTS ${${uppercase_name}_LIB} )
    get_filename_component(lib_dir ${${uppercase_name}_LIB} DIRECTORY)
    find_library(sfcio_path NAMES ${lib_name} PATHS ${lib_dir} NO_DEFAULT_PATH)
    if(NOT sfcio_path)
      set(lib_name ${name}) #_${version}_${kind})
      find_library(sfcio_path NAMES ${lib_name} PATHS ${lib_dir} NO_DEFAULT_PATH)
    endif()
    add_library(${lib_name} STATIC IMPORTED)
    set_target_properties(${lib_name} PROPERTIES
      IMPORTED_LOCATION ${sfcio_path}
      INTERFACE_INCLUDE_DIRECTORIES ${${uppercase_name}_INC})
  endif()
      
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sfcio
  REQUIRED_VARS sfcio_path)
