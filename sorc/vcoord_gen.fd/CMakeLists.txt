set(fortran_src
    driver.f90
    matrix_utils.f90
    vcoord_gen.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

set(exe_name vcoord_gen)
add_executable(${exe_name} ${fortran_src})

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
