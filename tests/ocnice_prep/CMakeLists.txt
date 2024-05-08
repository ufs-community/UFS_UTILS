# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Ed Hartnett

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

include_directories(${PROJECT_SOURCE_DIR})

# Copy necessary test files from the source data directory to the
# build data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

# Copy necessary test files from the source data directory to the
# build data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/invalid.model.nml ${CMAKE_CURRENT_BINARY_DIR}/data/invalid.model.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/invalid.srcdim.nml ${CMAKE_CURRENT_BINARY_DIR}/data/invalid.srcdim.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/invalid.dstdim.nml ${CMAKE_CURRENT_BINARY_DIR}/data/invalid.dstdim.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ice.badvecpairs.csv ${CMAKE_CURRENT_BINARY_DIR}/data/ice.badvecpairs.csv)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ocean.badvecpairs.csv ${CMAKE_CURRENT_BINARY_DIR}/data/ocean.badvecpairs.csv)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ice.badvecgrid.csv ${CMAKE_CURRENT_BINARY_DIR}/data/ice.badvecgrid.csv)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ocean.badvecgrid.csv ${CMAKE_CURRENT_BINARY_DIR}/data/ocean.badvecgrid.csv)

add_executable(ocnice_ftst_program_setup ftst_program_setup.F90)
add_test(NAME ocnice-ftst_program_setup COMMAND ocnice_ftst_program_setup)
target_link_libraries(ocnice_ftst_program_setup ocnice_prep_lib)
