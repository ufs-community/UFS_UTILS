# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

add_executable(ftst_enclosure_cnvx ftst_enclosure_cnvx.F90)
target_link_libraries(ftst_enclosure_cnvx lakefrac_lib)

add_executable(ftst_middle ftst_middle.F90)
target_link_libraries(ftst_middle lakefrac_lib)

add_test(NAME lake-ftst_enclosure_cnvx COMMAND ftst_enclosure_cnvx)
add_test(NAME lake-ftst_middle COMMAND ftst_middle)
