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

add_library(shared_lib STATIC ${c_src})
target_compile_definitions(shared_lib PRIVATE use_netCDF)

target_include_directories(shared_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(shared_lib NetCDF::NetCDF_C)
