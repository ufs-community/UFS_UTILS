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

# Copy necessary test files from the source data directory to the
# build data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/files.txt ${CMAKE_CURRENT_BINARY_DIR}/data/files.txt)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

add_executable(ftst_interp ftst_interp.F90)
target_link_libraries(ftst_interp sfc_climo_gen_lib)

# Cause test to be run with MPI.
add_mpi_test(sfc_climo_gen-ftst_interp
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_interp
  NUMPROCS 1
  TIMEOUT 60)
