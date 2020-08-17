set(fortran_src
    mtnlm7_oclsm.f
    netcdf_io.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian -fno-range-check")
endif()

set(exe_name orog)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  bacio_4
  w3nco_d
  ip_d
  sp_d
  NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
