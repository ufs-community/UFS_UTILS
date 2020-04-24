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
  set(CMAKE_Fortran_FLAGS "-g -traceback -r8 -i4 -convert big_endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-form -ffree-line-length-0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name sfc_climo_gen)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE ${ESMF_INC})
target_link_libraries(
  ${exe_name}
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
