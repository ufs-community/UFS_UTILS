set(fortran_src
    enclosure_cnvx.F90  
    find_limit.F90  
    lakefrac.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()

set(exe_name lakefrac)

add_executable(${exe_name} ${fortran_src})

target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
