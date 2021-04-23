set(lib_src
  bilinear_interp.c
  conserve_interp.c
  fregrid_util.c)

set(exe_src fregrid.c)

add_library(fregrid_lib STATIC ${lib_src})
add_executable(fregrid ${exe_src})


target_link_libraries(fregrid_lib
  PUBLIC
  m
  shared_lib
  NetCDF::NetCDF_C)

target_link_libraries(fregrid PRIVATE fregrid_lib)

install(TARGETS fregrid RUNTIME DESTINATION ${exec_dir})
