set(lib_src nb.F90)
set(exe_src inland.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()

add_library(inland_lib STATIC ${lib_src})
add_executable(inland ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(inland_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(inland_lib INTERFACE ${mod_dir})

target_link_libraries(
  inland_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(inland PRIVATE inland_lib)

install(TARGETS inland)
