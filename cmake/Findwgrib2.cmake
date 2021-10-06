# This module produces the target wgrib2::wgrib2

find_path(WGRIB2_INCLUDES wgrib2api.mod)
find_library(WGRIB2_LIB libwgrib2.a)
find_library(WGRIB2_API_LIB libwgrib2_api.a)

add_library(wgrib2::wgrib2 UNKNOWN IMPORTED)

# Library builds are different between CMake build and the make build.
# libwgrib2_api.a is only necessary in the CMake build and must come first when linking
if(WGRIB2_API_LIB)
  # CMake build. Need both.
  set(first_lib ${WGRIB2_API_LIB})
  set(second_lib ${WGRIB2_LIB})
else()
  # Makefile build. Only need libwgrib2.a
  set(first_lib ${WGRIB2_LIB})
  set(second_lib "")
endif()

set_target_properties(wgrib2::wgrib2 PROPERTIES
  IMPORTED_LOCATION "${first_lib}"
  INTERFACE_INCLUDE_DIRECTORIES "${WGRIB2_INCLUDES}"
  INTERFACE_LINK_LIBRARIES "${second_lib}")

set(WGRIB2_LIBRARIES "${first_lib}" "${second_lib}")

find_program(WGRIB2_EXE wgrib2)
execute_process(COMMAND ${WGRIB2_EXE} --version OUTPUT_VARIABLE version_str)

# Wgrib2 changed how it output --version from "v0.x.y.z" to "vx.y.z" starting in wgrib2 3.0
if(version_str MATCHES "^v0.*")
  string(SUBSTRING "${version_str}" 3 5 version)
else()
  string(SUBSTRING "${version_str}" 1 5 version)
endif()

find_package_handle_standard_args(wgrib2
  REQUIRED_VARS WGRIB2_LIBRARIES WGRIB2_INCLUDES WGRIB2_EXE
  VERSION_VAR version
  )
