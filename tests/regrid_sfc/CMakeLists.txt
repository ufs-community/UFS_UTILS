# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

# Copy necessary test files from the source data directory to the
# build data directory.

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/regrid.nml ${CMAKE_CURRENT_BINARY_DIR}/regrid.nml)

add_executable(ftst_read_namelist ftst_read_namelist.F90)
target_link_libraries(ftst_read_namelist regrid_lib)

add_test(NAME regrid_sfc-ftst_read_namelist COMMAND ftst_read_namelist)
