set(fortran_src
    mkgfsnemsioctl.f90)

set(exe_name mkgfsnemsioctl)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  nemsio::nemsio)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
