list(APPEND fortran_src
  scrip.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i4 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

set(exe_name weight_gen)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

if(sp_FOUND)
  target_link_libraries(${exe_name} sp::sp_d)
else()
  target_link_libraries(${exe_name} ip::ip_d)
endif()

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()
