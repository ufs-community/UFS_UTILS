# This is the CMake file for the fvcom_tools utility in the UFS_UTILS
# project.
#
# David Wright

set(fortran_src
   kinds.f90
   module_ncio.f90
   module_nwp_base.f90
   module_nwp.f90)
#   process_FVCOM.f90)

set(exe_src process_FVCOM.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

set(exe_name fvcom_to_FV3)

add_library(fvcom_tools_lib STATIC ${fortran_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(fvcom_tools_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(fvcom_tools_lib INTERFACE ${mod_dir})

target_link_libraries(
#  ${exe_name}
  fvcom_tools_lib
  PUBLIC
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

target_link_libraries(${exe_name} PRIVATE fvcom_tools_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()

