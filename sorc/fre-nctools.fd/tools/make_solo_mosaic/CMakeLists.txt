set(c_src
    make_solo_mosaic.c
    get_contact.c)

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-g -traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " ")
endif()

add_executable(make_solo_mosaic ${c_src})

target_link_libraries(make_solo_mosaic
  m
  shared_all
  NetCDF::NetCDF_C)

install(TARGETS make_solo_mosaic RUNTIME DESTINATION ${exec_dir})
