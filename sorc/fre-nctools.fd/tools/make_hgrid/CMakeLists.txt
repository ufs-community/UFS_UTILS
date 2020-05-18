set(c_src
  create_conformal_cubic_grid.c
  create_gnomonic_cubic_grid.c
  create_grid_from_file.c
  create_lonlat_grid.c
  make_hgrid.c)

add_executable(make_hgrid ${c_src})

target_link_libraries(make_hgrid
  m
  shared_lib
  NetCDF::NetCDF_C)

install(TARGETS make_hgrid RUNTIME DESTINATION ${exec_dir})
