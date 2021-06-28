# This is the CMake build file for the emcsfc_ice_blend utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts

set(fortran_src
    emcsfc_ice_blend.f90)

set(exe_name emcsfc_ice_blend)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  bacio::bacio_4
  g2::g2_4
  w3nco::w3nco_4)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
