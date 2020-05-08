set(fortran_src
    emcsfc_ice_blend.f90)

set(exe_name emcsfc_ice_blend)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  bacio_4
  g2_4
  w3nco_4
  ${JASPER_LIBRARIES}
  PNG::PNG
  ZLIB::ZLIB)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
