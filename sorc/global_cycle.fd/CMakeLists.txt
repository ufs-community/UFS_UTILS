set(fortran_src
    cycle.f90     
    machine.f90
    num_parthds.f90
    sfcsub.F
    read_write_data.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name global_cycle)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  w3nco_d
  bacio_4
  ip_d
  sp_d
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran
  OpenMP::OpenMP_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
