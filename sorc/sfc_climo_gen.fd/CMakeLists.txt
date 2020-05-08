set(fortran_src
    driver.F90
    interp.F90
    model_grid.F90
    output.f90
    program_setup.f90
    search.f90
    source_grid.F90
    utils.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name sfc_climo_gen)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
