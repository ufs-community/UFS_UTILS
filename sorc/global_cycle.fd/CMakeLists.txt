# This is the CMake file for the global_cycle utility in the UFS_UTILS
# project.
#
# George Gayno, Kyle Gerheiser

set(lib_src
    machine.f90
    num_parthds.f90
    ../../ccpp-physics/physics/sfcsub.F
    read_write_data.f90
    utils.F90
    land_increments.f90)

set(exe_src cycle.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(exe_name global_cycle)

add_library(global_cycle_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

add_dependencies(${exe_name} noah_lsm_routines)

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(global_cycle_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(global_cycle_lib INTERFACE ${mod_dir})

target_link_libraries(
  global_cycle_lib
  PUBLIC
  w3nco::w3nco_d
  bacio::bacio_4
  ip::ip_d
  sp::sp_d
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran
  noah_lsm_routines)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(global_cycle_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE global_cycle_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
