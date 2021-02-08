# This is the cmake file to build unit tests for the UFS_UTILS
# project.
#
# Ed Hartnett, 2/8/21

add_executable(ftst_input_data ftst_input_data.F90)
add_test(ftst_input_data ftst_input_data)
