set(fortran_src
    mtnlm7_oclsm.f
    netcdf_io.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -r8 -convert big_endian -assume byterecl")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian -fno-range-check")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name orog)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  bacio_4
  w3nco_d
  ip_d
  sp_d
  NetCDF::NetCDF_Fortran
  OpenMP::OpenMP_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
