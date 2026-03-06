set(lib_src
    version.c)

set(exe_src mppnccombine.c)

add_library(mppnccombine_lib STATIC ${lib_src})
add_executable(mppnccombine ${exe_src})

target_link_libraries(mppnccombine_lib
  PUBLIC
  m
  shared_lib
  NetCDF::NetCDF_C)

target_link_libraries(mppnccombine PRIVATE mppnccombine_lib)

install(TARGETS mppnccombine)
