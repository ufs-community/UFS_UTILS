# This is the CMake file for the regrid_sfc utility in the UFS_UTILS
# project.
#
# George Gayno, Clara Draper

set(lib_src
   grids_IO.f90
   utils.f90)

set(exe_src  regridStates.f90)

# coppied from global_cycle 
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffree-line-length-0")
endif()

set(exe_name regridStates.x)

add_library(regrid_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(regrid_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(regrid_lib INTERFACE ${mod_dir})
target_link_libraries(
  regrid_lib
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(regrid_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE regrid_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
