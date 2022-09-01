set(lib_src
    module_gsl_oro_data_lg_scale.f90
    module_gsl_oro_data_sm_scale.f90)

set(exe_src gsl_oro_data.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fno-range-check")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
endif()

add_library(orog_gsl_lib STATIC ${lib_src})
set(exe_name orog_gsl)
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(orog_gsl_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(orog_gsl_lib INTERFACE ${mod_dir})

target_link_libraries(
  orog_gsl_lib
  PUBLIC
  NetCDF::NetCDF_Fortran)

target_link_libraries(${exe_name} PRIVATE orog_gsl_lib)

install(TARGETS ${exe_name})
