# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett, Larissa Reames

# Include cmake to allow parallel I/O tests.
include (LibMPI)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
endif()

include_directories(${PROJECT_SOURCE_DIR})

# Copy necessary test files.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_gaussian_nemsio.nml ${CMAKE_CURRENT_BINARY_DIR}/fort.41)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_fv3_tiled.nml ${CMAKE_CURRENT_BINARY_DIR}/config_fv3_tiled.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_fv3_tiled_warm_restart.nml ${CMAKE_CURRENT_BINARY_DIR}/config_fv3_tiled_warm_restart.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_gaussian_nemsio.nml ${CMAKE_CURRENT_BINARY_DIR}/config_gaussian_nemsio.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_gfs_grib2.nml ${CMAKE_CURRENT_BINARY_DIR}/config_gfs_grib2.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config_spectral_sigio.nml ${CMAKE_CURRENT_BINARY_DIR}/config_spectral_sigio.nml)

add_executable(ftst_utils ftst_utils.F90)
add_test(NAME chgres_cube-ftst_utils COMMAND ftst_utils)
target_link_libraries(ftst_utils chgres_cube_lib)

add_executable(ftst_program_setup ftst_program_setup.F90)
target_link_libraries(ftst_program_setup chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_program_setup
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_program_setup
  NUMPROCS 4
  TIMEOUT 60)

add_executable(ftst_search_util ftst_search_util.F90)
target_link_libraries(ftst_search_util chgres_cube_lib)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(ftst_search_util OpenMP::OpenMP_Fortran)
endif()

# Cause test to be run with MPI. The search routine was
# written to work on a single task. So set NUMPROCS to 1.
add_mpi_test(chgres_cube-ftst_search_util
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_search_util
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_sfc_input_data ftst_sfc_input_data.F90)
target_link_libraries(ftst_sfc_input_data chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_sfc_input_data
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_sfc_input_data
  NUMPROCS 4
  TIMEOUT 60)
