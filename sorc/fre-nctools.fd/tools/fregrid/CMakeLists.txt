set(c_src
  bilinear_interp.c
  conserve_interp.c
  fregrid.c
  fregrid_util.c)

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
    
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " ")
endif()

add_executable(fregrid ${c_src})

target_link_libraries(fregrid
  shared_all
  NetCDF::NetCDF_C)

install(TARGETS fregrid RUNTIME DESTINATION exec)
