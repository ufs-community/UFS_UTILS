set(lib_src nb.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()

add_library(inland_lib ${lib_src})
add_executable(inland inland.F90)

target_link_libraries(
  inland_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(inland PRIVATE inland_lib)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
