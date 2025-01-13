# This is the cmake build file.
#
# George Gayno, Ed Hartnett

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

# Copy necessary test files from the source data directory to the
# build test directory.

# The "read_global_mask" test needs this file. Note that the copy command changes the
# file name to "landcover.umd.30s.nc", which is the name expected by the "read_global_mask" routine.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/landcover.umd.lowres.nc ${CMAKE_CURRENT_BINARY_DIR}/landcover.umd.30s.nc)

# The "read_global_orog" test needs this file. Note that the copy command changes the
# file name to "topography.gmted2010.30s.nc", which is the name expected by the "read_global_orog" routine.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/topography.gmted2010.lowres.nc ${CMAKE_CURRENT_BINARY_DIR}/topography.gmted2010.30s.nc)

# The "read_mdl_dims" and "read_mdl_grid_file" tests expects this file.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/C12_grid.tile1.nc ${CMAKE_CURRENT_BINARY_DIR}/C12_grid.tile1.nc)

# The "read_mask" test expects this file.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/C48.mx500.tile1.nc ${CMAKE_CURRENT_BINARY_DIR}/C48.mx500.tile1.nc)

# The "qc_orog_by_ramp" test needs this file. Note that the copy command changes the
# file name to "topography.antarctica.ramp.30s.nc", which is the name expected by the "qc_orog_by_ramp" routine.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/topography.antarctica.ramp.lowres.nc ${CMAKE_CURRENT_BINARY_DIR}/topography.antarctica.ramp.30s.nc)

add_executable(ftst_ll2xyz ftst_ll2xyz.F90)
add_test(NAME orog-ftst_ll2xyz COMMAND ftst_ll2xyz)
target_link_libraries(ftst_ll2xyz orog_lib)

add_executable(ftst_minmax ftst_minmax.F90)
add_test(NAME orog-ftst_minmax COMMAND ftst_minmax)
target_link_libraries(ftst_minmax orog_lib)

add_executable(ftst_get_ll_angle ftst_get_ll_angle.F90)
add_test(NAME orog-ftst_get_ll_angle COMMAND ftst_get_ll_angle)
target_link_libraries(ftst_get_ll_angle orog_lib)

add_executable(ftst_get_index ftst_get_index.F90)
add_test(NAME orog-ftst_get_index COMMAND ftst_get_index)
target_link_libraries(ftst_get_index orog_lib)

add_executable(ftst_inside_polygon ftst_inside_polygon.F90)
add_test(NAME orog-ftst_inside_polygon COMMAND ftst_inside_polygon)
target_link_libraries(ftst_inside_polygon orog_lib)

add_executable(ftst_transpose ftst_transpose.F90)
add_test(NAME orog-ftst_transpose COMMAND ftst_transpose)
target_link_libraries(ftst_transpose orog_lib)

add_executable(ftst_find_poles ftst_find_poles.F90)
add_test(NAME orog-ftst_find_poles COMMAND ftst_find_poles)
target_link_libraries(ftst_find_poles orog_lib)

add_executable(ftst_find_nearest_pole_pts ftst_find_nearest_pole_pts.F90)
add_test(NAME orog-ftst_find_nearest_pole_pts COMMAND ftst_find_nearest_pole_pts)
target_link_libraries(ftst_find_nearest_pole_pts orog_lib)

add_executable(ftst_get_xnsum ftst_get_xnsum.F90)
add_test(NAME orog-ftst_get_xnsum COMMAND ftst_get_xnsum)
target_link_libraries(ftst_get_xnsum orog_lib)

add_executable(ftst_get_xnsum2 ftst_get_xnsum2.F90)
add_test(NAME orog-ftst_get_xnsum2 COMMAND ftst_get_xnsum2)
target_link_libraries(ftst_get_xnsum2 orog_lib)

add_executable(ftst_get_xnsum3 ftst_get_xnsum3.F90)
add_test(NAME orog-ftst_get_xnsum3 COMMAND ftst_get_xnsum3)
target_link_libraries(ftst_get_xnsum3 orog_lib)

add_executable(ftst_rm_isolated_pts ftst_rm_isolated_pts.F90)
add_test(NAME orog-ftst_rm_isolated_pts COMMAND ftst_rm_isolated_pts)
target_link_libraries(ftst_rm_isolated_pts orog_lib)

add_executable(ftst_read_global_mask ftst_read_global_mask.F90)
add_test(NAME orog-ftst_read_global_mask COMMAND ftst_read_global_mask)
target_link_libraries(ftst_read_global_mask orog_lib)

add_executable(ftst_read_global_orog ftst_read_global_orog.F90)
add_test(NAME orog-ftst_read_global_orog COMMAND ftst_read_global_orog)
target_link_libraries(ftst_read_global_orog orog_lib)

add_executable(ftst_read_mdl_dims ftst_read_mdl_dims.F90)
add_test(NAME orog-ftst_read_mdl_dims COMMAND ftst_read_mdl_dims)
target_link_libraries(ftst_read_mdl_dims orog_lib)

add_executable(ftst_read_mdl_grid_file ftst_read_mdl_grid_file.F90)
add_test(NAME orog-ftst_read_mdl_grid_file COMMAND ftst_read_mdl_grid_file)
target_link_libraries(ftst_read_mdl_grid_file orog_lib)

add_executable(ftst_read_mask ftst_read_mask.F90)
add_test(NAME orog-ftst_read_mask COMMAND ftst_read_mask)
target_link_libraries(ftst_read_mask orog_lib)

add_executable(ftst_qc_orog_by_ramp ftst_qc_orog_by_ramp.F90)
add_test(NAME orog-ftst_qc_orog_by_ramp COMMAND ftst_qc_orog_by_ramp)
target_link_libraries(ftst_qc_orog_by_ramp orog_lib)
