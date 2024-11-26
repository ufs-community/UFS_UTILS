set(lib_src read_write.F90 merge.F90 utils.F90 namelist.F90)
set(exe_src merge_lake_ocnmsk.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i4 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(exe_name ocean_merge)

add_library(om_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

target_link_libraries(
  om_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(${exe_name} PRIVATE om_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()
