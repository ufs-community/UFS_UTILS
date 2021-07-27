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
    ${CMAKE_CURRENT_SOURCE_DIR}/data/input.nml ${CMAKE_CURRENT_BINARY_DIR}/input.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

add_executable(ftst_read_filter_topo_nml ftst_readnml.F90)
add_test(NAME filter_topo-ftst_read_namelist COMMAND ftst_read_filter_topo_nml)
target_link_libraries(ftst_read_filter_topo_nml filter_topo_lib)
