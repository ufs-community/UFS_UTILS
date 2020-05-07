set(c_src
  create_conformal_cubic_grid.c
  create_gnomonic_cubic_grid.c
  create_grid_from_file.c
  create_lonlat_grid.c
  make_hgrid.c)

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
    
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " ")
endif()

add_executable(make_hgrid ${c_src})

target_link_libraries(make_hgrid
  m
  shared_lib
  NetCDF::NetCDF_C)

install(TARGETS make_hgrid RUNTIME DESTINATION ${exec_dir})
