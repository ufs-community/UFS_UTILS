set(fortran_src
    nemsio_read.f90)

set(exe_name nemsio_read)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  nemsio
  bacio_4
  w3nco_d)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
