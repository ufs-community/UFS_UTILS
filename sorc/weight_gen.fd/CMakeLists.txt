list(APPEND fortran_src
  scrip.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i4 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(exe_name scrip.x)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  sp::sp_d
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
