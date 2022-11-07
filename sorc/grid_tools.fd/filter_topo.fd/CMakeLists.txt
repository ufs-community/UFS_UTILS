set(lib_src
    utils.F90)

set(exe_src
    filter_topo.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl -real-size 64 -fno-alias -stack-temps -safe-cray-ptr -ftz")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

set(exe_name filter_topo)

add_library(filter_topo_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(filter_topo_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(filter_topo_lib INTERFACE ${mod_dir})

target_link_libraries(
  filter_topo_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(${exe_name} PRIVATE filter_topo_lib)

install(TARGETS ${exe_name})
