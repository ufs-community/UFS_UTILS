# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett

set(fortran_src
    "${CMAKE_SOURCE_DIR}/sorc/chgres_cube.fd/utils.f90"
    to_upper_lower.F90)

include_directories(
      ${PROJECT_SOURCE_DIR}
)

add_test(NAME to_upper_lower COMMAND to_upper_lower)


