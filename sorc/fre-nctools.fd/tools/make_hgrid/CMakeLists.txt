set(c_src
  create_conformal_cubic_grid.c
  create_gnomonic_cubic_grid.c
  create_grid_from_file.c
  create_lonlat_grid.c
  make_hgrid.c)

set ( PROJECT_LINK_LIBS libshared_all.a )
link_directories(${CMAKE_CURRENT_BINARY_DIR}/../../shared_lib)
include_directories(../../shared_lib)

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG
      "-O0"
  )
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " "
  )
endif()

set(exe_name make_hgrid)
add_executable(${exe_name} ${c_src})
target_include_directories(
  ${exe_name} PRIVATE)
target_link_libraries(
  ${exe_name}
  ${PROJECT_LINK_LIBS}
  NetCDF::NetCDF_C)

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
