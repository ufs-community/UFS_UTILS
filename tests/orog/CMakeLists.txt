# This is the cmake build file.
#
# George Gayno, Ed Hartnett

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

include_directories(${PROJECT_SOURCE_DIR})

add_executable(ftst_ll2xyz ftst_ll2xyz.F90)
add_test(NAME orog-ftst_ll2xyz COMMAND ftst_ll2xyz)
target_link_libraries(ftst_ll2xyz orog_lib)

add_executable(ftst_minmax ftst_minmax.F90)
add_test(NAME orog-ftst_minmax COMMAND ftst_minmax)
target_link_libraries(ftst_minmax orog_lib)
