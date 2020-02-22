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
  set(CMAKE_Fortran_FLAGS "-g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -ffree-form -ffree-line-length-0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

set(exe_name chgres_cube.exe)
add_executable(${exe_name} ${fortran_src})
set_target_properties(${lib_name} PROPERTIES COMPILE_FLAGS
                                             "${OpenMP_Fortran_FLAGS}")
target_include_directories(
  ${exe_name}
  PRIVATE ${ESMF_INC} ${WGRIB2_INCLUDES} ${MPI_Fortran_INCLUDE_PATH}
          ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90})
target_link_libraries(
  ${exe_name}
  ${WGRIB2_LIBRARIES}
  ${NETCDF_LIBRARIES_F90}
  ${NETCDF_LIBRARIES}
  nemsio
  sfcio_4
  ${MPI_Fortran_LIBRARIES}
  sigio_4
  bacio_4
  sp_d
  w3nco_d
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${OpenMP_Fortran_FLAGS}
  esmf)

install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
