set(fortran_src
    filter_topo.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -assume byterecl -i4 -real_size 64")
  set(CMAKE_Fortran_FLAGS_RELEASE "-fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -O2")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-0 -fdefault-real-8")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name filter_topo)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE)
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
