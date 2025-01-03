# This is the cmake build file.
#
# George Gayno, Ed Hartnett

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -warn unused")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

#include_directories(${PROJECT_SOURCE_DIR})

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
