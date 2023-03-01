# This is the CMake build file for the emcsfc_snow2mdl utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    grib_utils.F90
    model_grid.F90
    program_setup.F90
    snow2mdl.F90
    snowdat.F90)

set(exe_src driver.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian")
endif()
if(ip_VERSION GREATER_EQUAL 4.0.0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIP_V4")
endif()

set(exe_name emcsfc_snow2mdl)

add_library(snow2mdl_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(snow2mdl_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(snow2mdl_lib INTERFACE ${mod_dir})

target_link_libraries(snow2mdl_lib
  PUBLIC
  g2::g2_d
  ip::ip_d
  sp::sp_d
  bacio::bacio_4
  w3nco::w3nco_d)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE snow2mdl_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
