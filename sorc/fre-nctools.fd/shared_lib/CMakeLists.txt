set(c_src
    affinity.c
    create_xgrid.c
    gradient_c2l.c
    interp.c
    mosaic_util.c
    mpp.c
    mpp_domain.c
    mpp_io.c
    mpp_domain.c
    read_mosaic.c
    tool_util.c)

add_definitions(-Duse_netCDF)

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

set(lib_name shared_all)
add_library(${lib_name} STATIC ${c_src})
target_include_directories(
  ${lib_name} PRIVATE)
target_link_libraries(
  ${lib_name}
  NetCDF::NetCDF_C)

install(TARGETS ${lib_name} RUNTIME DESTINATION ./)
