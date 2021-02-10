set(fortran_src
    "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/utils.f90"
    to_upper_lower.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
endif()

include_directories(
      ${PROJECT_SOURCE_DIR}
)

set(exe_name to_upper_lower)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
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
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
