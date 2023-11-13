# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett, Larissa Reames

# Include cmake to allow parallel I/O tests.
include (LibMPI)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

include_directories(${PROJECT_SOURCE_DIR})

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.001 ${CMAKE_CURRENT_BINARY_DIR}/xainc.001)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.002 ${CMAKE_CURRENT_BINARY_DIR}/xainc.002)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.003 ${CMAKE_CURRENT_BINARY_DIR}/xainc.003)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.004 ${CMAKE_CURRENT_BINARY_DIR}/xainc.004)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.005 ${CMAKE_CURRENT_BINARY_DIR}/xainc.005)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/xainc.006 ${CMAKE_CURRENT_BINARY_DIR}/xainc.006)

add_executable(ftst_land_increments ftst_land_increments.F90)
target_link_libraries(ftst_land_increments global_cycle_lib)
add_executable(ftst_read_increments ftst_read_increments.F90)
target_link_libraries(ftst_read_increments global_cycle_lib)

add_mpi_test(global_cycle-ftst_land_increments
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_land_increments
  NUMPROCS 1 TIMEOUT 60)
add_mpi_test(global_cycle-ftst_read_increments
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_increments
  NUMPROCS 6 TIMEOUT 60)

