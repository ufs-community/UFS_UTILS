set(fortran_src
    driver.f
    grib_utils.f
    model_grid.f
    program_setup.f
    snow2mdl.f
    snowdat.f)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -FR -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name emcsfc_snow2mdl)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  g2::g2_d
  ip::ip_d
  sp::sp_d
  landsfcutil::landsfcutil_d
  bacio::bacio_4
  w3nco::w3nco_d)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
