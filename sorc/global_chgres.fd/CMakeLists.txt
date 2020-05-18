set(fortran_src
    chgres.f90
    chgres_utils.f90
    funcphys.f90
    machine_8.f90
    nrlmsise00_sub.f90
    nsst_chgres.f90
    nstio_module.f90
    num_parthds.f90
    physcons.f90
    read_write.f90
    sfcsub.F
    surface_chgres.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name global_chgres)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  gfsio_4
  nemsiogfs
  nemsio
  sigio_4
  sfcio_4
  landsfcutil_d
  ip_d
  sp_d
  w3emc_d
  w3nco_d
  bacio_4
  NetCDF::NetCDF_Fortran
  OpenMP::OpenMP_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
