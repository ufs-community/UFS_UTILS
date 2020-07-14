set(fortran_src
    global_equiv_resol.f90)

set(exe_name global_equiv_resol)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
