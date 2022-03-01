set(c_srcs mppnccombine.c)

set(exe_name mppnccombine)
add_executable(${exe_name} ${c_srcs})

target_link_libraries(${exe_name}
  PUBLIC
  m
  shared_lib
  NetCDF::NetCDF_C)

install(TARGETS ${exe_name} DESTINATION ${exec_dir})
