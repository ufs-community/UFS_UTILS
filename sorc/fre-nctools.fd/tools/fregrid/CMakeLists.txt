set(c_src
  bilinear_interp.c
  conserve_interp.c
  fregrid.c
  fregrid_util.c)

add_executable(fregrid ${c_src})

target_link_libraries(fregrid
  m
  shared_lib
  NetCDF::NetCDF_C)

install(TARGETS fregrid RUNTIME DESTINATION ${exec_dir})
