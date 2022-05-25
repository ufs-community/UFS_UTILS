# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

# Copy necessary test files from the source data directory to the
# build data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/data/config.nml ${CMAKE_CURRENT_BINARY_DIR}/fort.41)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

add_executable(ftst_program_setup_snow ftst_program_setup.F90)
target_link_libraries(ftst_program_setup_snow snow2mdl_lib)
add_test(NAME emcsfc_snow2mdl-ftst_program_setup COMMAND ftst_program_setup_snow)
