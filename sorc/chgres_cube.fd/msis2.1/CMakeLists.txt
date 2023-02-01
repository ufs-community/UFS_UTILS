include(FetchContent)
include(CheckFortranCompilerFlag)

add_library(msis2 msis_utils.F90 msis_constants.F90 msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90 msis_gtd8d.F90)
set_target_properties(msis2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(msis2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # msis_calc:bspline has argument mismatch on nodes variable
  target_compile_options(msis2 PRIVATE -std=legacy)
endif()

if(BUILD_TESTING)
  add_executable(msis2test msis2.1_test.F90)
  target_link_libraries(msis2test PRIVATE msis2)

  add_test(NAME MSIS2
    COMMAND $<TARGET_FILE:msis2test>
    WORKING_DIRECTORY ${msis2proj_SOURCE_DIR})
endif()
