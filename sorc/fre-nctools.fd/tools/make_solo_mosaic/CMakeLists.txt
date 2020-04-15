set(c_src
    make_solo_mosaic.c
    get_contact.c
    ../../shared/mosaic/mosaic_util.c
    ../../shared/mosaic/read_mosaic.c
    ../../shared/mosaic/interp.c
    ../shared/mpp_io.c
    ../shared/mpp_domain.c
    ../shared/mpp.c
    ../shared/tool_util.c 
    ../../shared/mosaic/create_xgrid.c
    )

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-traceback")
  set(CMAKE_C_FLAGS_RELEASE "-Duse_netCDF -O2 -I${CMAKE_CURRENT_SOURCE_DIR}/../shared -I${CMAKE_CURRENT_SOURCE_DIR}/../../shared/mosaic")
  set(CMAKE_C_FLAGS_DEBUG
      "-O0"
  )
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " "
  )
endif()

set(exe_name make_solo_mosaic)
add_executable(${exe_name} ${c_src})
target_include_directories(
  ${exe_name} PRIVATE)
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_C)

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
