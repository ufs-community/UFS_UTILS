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

set(exe_name regional_grid)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
