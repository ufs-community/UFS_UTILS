# This is the CMake build file for the emcsfc_snow2mdl utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts

set(fortran_src
    driver.F90
    grib_utils.F90
    model_grid.F90
    program_setup.F90
    snow2mdl.F90
    snowdat.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name emcsfc_snow2mdl)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  g2::g2_d
  ip::ip_d
  sp::sp_d
  bacio::bacio_4
  w3nco::w3nco_d)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
