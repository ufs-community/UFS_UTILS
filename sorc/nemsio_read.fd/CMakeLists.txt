set(fortran_src
    nemsio_read.f90)

set(exe_name nemsio_read)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  nemsio::nemsio)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
