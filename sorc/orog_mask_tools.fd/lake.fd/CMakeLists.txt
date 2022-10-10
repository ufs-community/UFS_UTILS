set(lib_src
    enclosure_cnvx.F90
    find_limit.F90)

set(exe_src lakefrac.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()

add_library(lakefrac_lib STATIC ${lib_src})
add_executable(lakefrac ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(lakefrac_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(lakefrac_lib INTERFACE ${mod_dir})

target_link_libraries(
  lakefrac_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(lakefrac PRIVATE lakefrac_lib)

install(TARGETS lakefrac)
