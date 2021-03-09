# This is the CMake file for the fvcom_tools utility in the UFS_UTILS
# project.
#
# David Wright

set(fortran_src
   kinds.f90
   module_ncio.f90
   module_nwp_base.f90
   module_nwp.f90
   process_FVCOM.f90)


if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
endif()

set(exe_name fvcom_to_FV3)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
