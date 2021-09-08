# This module produces the target wgrib2::wgrib2

find_path(WGRIB2_INCLUDES wgrib2api.mod)

find_library(WGRIB2_LIBRARIES libwgrib2.a)

find_program(WGRIB2_EXE wgrib2)
execute_process(COMMAND ${WGRIB2_EXE} --version OUTPUT_VARIABLE version_str)
string(SUBSTRING ${version_str} 3 5 version)

mark_as_advanced(WGRIB2_INCLUDES WGRIB2_LIBRARIES)

add_library(wgrib2::wgrib2 UNKNOWN IMPORTED)

set_target_properties(wgrib2::wgrib2 PROPERTIES
  IMPORTED_LOCATION "${WGRIB2_LIBRARIES}"
  INTERFACE_INCLUDE_DIRECTORIES "${WGRIB2_INCLUDES}"
  INTERFACE_LINK_LIBRARIES "${WGRIB2_LIBRARIES}")

find_package_handle_standard_args(wgrib2
  REQUIRED_VARS WGRIB2_LIBRARIES WGRIB2_INCLUDES WGRIB2_EXE
  VERSION_VAR version
)
