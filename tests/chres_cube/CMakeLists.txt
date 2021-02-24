# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
endif()

include_directories(${PROJECT_SOURCE_DIR})

# Copy necessary test files.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/config.nml
    ${CMAKE_CURRENT_BINARY_DIR}/fort.41)

add_executable(ftst_utils "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/utils.f90" ftst_utils.F90)
add_test(NAME ftst_utils COMMAND ftst_utils)
target_link_libraries(
  ftst_utils
  nemsio::nemsio
  sfcio::sfcio
  sigio::sigio
  bacio::bacio_4
  sp::sp_d
  w3nco::w3nco_d
  esmf
  wgrib2::wgrib2_lib
  wgrib2::wgrib2_api
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(ftst_utils OpenMP::OpenMP_Fortran)
endif()

add_executable(ftst_program_setup "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/input_data.F90"
  "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/model_grid.F90"
  "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/program_setup.f90"
  "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/grib2_util.F90"
  "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/search_util.f90"
  "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/utils.f90"
  ftst_program_setup.F90)
add_test(NAME ftst_program_setup COMMAND ftst_program_setup)
target_link_libraries(
  ftst_program_setup
  nemsio::nemsio
  sfcio::sfcio
  sigio::sigio
  bacio::bacio_4
  sp::sp_d
  w3nco::w3nco_d
  esmf
  wgrib2::wgrib2_lib
  wgrib2::wgrib2_api
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(ftst_program_setup OpenMP::OpenMP_Fortran)
endif()



