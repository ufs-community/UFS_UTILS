set(lib_src
  create_conformal_cubic_grid.c
  create_gnomonic_cubic_grid.c
  create_grid_from_file.c
  create_lonlat_grid.c)

add_library(make_hgrid_lib STATIC ${lib_src})
add_executable(make_hgrid make_hgrid.c)

target_link_libraries(make_hgrid_lib
  PUBLIC
  m
  shared_lib
  NetCDF::NetCDF_C)

target_link_libraries(make_hgrid PRIVATE make_hgrid_lib)

install(TARGETS make_hgrid)
