set(fortran_src1
    enclosure_cnvx.F90  
    find_limit.F90  
    lakefrac.F90)

set(fortran_src2
    nb.F90  
    inland.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()

set(exe_name1 lakefrac)
set(exe_name2 inland)

add_executable(${exe_name1} ${fortran_src1})
add_executable(${exe_name2} ${fortran_src2})

target_link_libraries(
  ${exe_name1}
  NetCDF::NetCDF_Fortran)

target_link_libraries(
  ${exe_name2}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name1} ${exe_name2} RUNTIME DESTINATION ${exec_dir})
