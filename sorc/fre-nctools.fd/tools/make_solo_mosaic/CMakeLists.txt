set(lib_src
    make_solo_mosaic.c
    get_contact.c)

set(exe_src make_solo_mosaic.c)

add_library(make_solo_mosaic_lib STATIC ${lib_src})
add_executable(make_solo_mosaic ${exe_src})

target_link_libraries(make_solo_mosaic_lib
  PUBLIC
  m
  shared_lib
  NetCDF::NetCDF_C)

target_link_libraries(make_solo_mosaic PRIVATE make_solo_mosaic_lib)

install(TARGETS make_solo_mosaic)
