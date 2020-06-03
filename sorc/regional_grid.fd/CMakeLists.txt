set(fortran_src
    gen_schmidt.f90
    hgrid_ak.f90
    pietc.f90
    pkind.f90
    pmat4.f90
    pmat5.f90
    pmat.f90
    psym2.f90
    regional_grid.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian -fno-range-check")
endif()

set(exe_name regional_grid)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
