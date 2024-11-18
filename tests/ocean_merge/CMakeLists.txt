# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

# Copy necessary test files from the source data directory to the
# build data directory.

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/input.nml ${CMAKE_CURRENT_BINARY_DIR}/input.nml)

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/oro.C48.tile1.nc ${CMAKE_CURRENT_BINARY_DIR}/oro.C48.tile1.nc)

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/C48.mx500.tile1.nc ${CMAKE_CURRENT_BINARY_DIR}/C48.mx500.tile1.nc)

add_executable(ftst_read_nml ftst_read_nml.F90)
target_link_libraries(ftst_read_nml om_lib)

add_executable(ftst_merge ftst_merge.F90)
target_link_libraries(ftst_merge om_lib)

add_executable(ftst_read_lake_mask ftst_read_lake_mask.F90)
target_link_libraries(ftst_read_lake_mask om_lib)

add_executable(ftst_read_ocean_frac ftst_read_ocean_frac.F90)
target_link_libraries(ftst_read_ocean_frac om_lib)

add_test(NAME ocean_merge-ftst_read_nml COMMAND ftst_read_nml)
add_test(NAME ocean_merge-ftst_merge COMMAND ftst_merge)
add_test(NAME ocean_merge-ftst_read_lake_mask COMMAND ftst_read_lake_mask)
add_test(NAME ocean_merge-ftst_read_ocean_frac COMMAND ftst_read_ocean_frac)
