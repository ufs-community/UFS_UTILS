set(fortran_src
    atmosphere.F90
    chgres.F90
    grib2_util.F90
    input_data.F90
    model_grid.F90
    program_setup.f90
    search_util.f90
    static_data.F90
    surface.F90
    utils.f90
    write_data.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
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

set(exe_name chgres_cube.exe)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE ${ESMF_INC} ${WGRIB2_INCLUDES} ${NETCDF_INCLUDES})
target_link_libraries(
  ${exe_name}
  nemsio
  sfcio_4
  sigio_4
  bacio_4
  sp_d
  w3nco_d
  esmf
  MPI::MPI_Fortran
  ${WGRIB2_LIBRARIES}
  ${NETCDF_LIBRARIES})
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
