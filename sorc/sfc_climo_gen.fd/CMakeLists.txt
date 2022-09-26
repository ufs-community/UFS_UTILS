# This is the CMake build file for the sfc_climo_gen utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    interp.F90
    model_grid.F90
    output.f90
    program_setup.f90
    search.f90
    source_grid.F90
    utils.f90)

set(exe_src driver.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

set(exe_name sfc_climo_gen)

add_library(sfc_climo_gen_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(sfc_climo_gen_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(sfc_climo_gen_lib INTERFACE ${mod_dir})

target_link_libraries(
  sfc_climo_gen_lib
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

target_link_libraries(${exe_name} PRIVATE sfc_climo_gen_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
