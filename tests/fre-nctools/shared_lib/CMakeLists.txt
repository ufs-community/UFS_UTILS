# This is the cmake build file for the fre-nctools shared library in
# the UFS_UTILS project.
#
# Ed Hartnett, 2/17/21

include_directories(${CMAKE_SOURCE_DIR}/sorc/fre-nctools.fd/shared_lib)

add_executable(tst_create_xgrid tst_create_xgrid.c)
add_test(NAME fre-nctools-tst_create_xgrid COMMAND tst_create_xgrid)
target_link_libraries(tst_create_xgrid NetCDF::NetCDF_C shared_lib m)






