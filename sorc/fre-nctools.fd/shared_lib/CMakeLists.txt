set(c_src
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

if(CMAKE_C_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_C_FLAGS "-traceback")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS " ")
  set(CMAKE_C_FLAGS_RELEASE " ")
  set(CMAKE_C_FLAGS_DEBUG " ")
endif()

add_library(shared_all STATIC ${c_src})
target_compile_definitions(shared_all PRIVATE use_netCDF)

target_include_directories(shared_all PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(shared_all NetCDF::NetCDF_C)

install(TARGETS shared_all RUNTIME DESTINATION lib)
