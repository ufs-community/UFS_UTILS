set(c_src
    make_solo_mosaic.c
    get_contact.c)

add_executable(make_solo_mosaic ${c_src})

target_link_libraries(make_solo_mosaic
  m
  shared_lib
  NetCDF::NetCDF_C)

install(TARGETS make_solo_mosaic RUNTIME DESTINATION ${exec_dir})
