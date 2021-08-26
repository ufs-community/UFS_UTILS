HOW TO CREATE A UNIT TEST.

Unit tests should test only small parts of a program. For example,
a single routine or function.

Source code for the test shall be placed in a program
specific directory under ./tests.

The name of the unit test source file shall begin with
"ftst_" (for Fortran) or "tst_" (for C/C++).

A unit test should have the following components.
For an example, see ./filter_topo/ftst_readnml.F90

1) A prolog describing the test and its author.

2) The test should begin with a print to standard
output, such as "Starting test of ..."

3) A call to the function or routine to be tested.

4) Output from the function or routine should be 
compared to expected values. If wrong values are
returned, the test should stop with a bad status.

5) If the test is successful, that message
should be printed to standard output, i.e.,
print*,"SUCCESS!"

How to compile under CMake. The above test is
compiled and run with these statements in the
CMakeLists.txt file. See ./filter_topo/CMakeLists.txt.

This names the unit test executable:
  add_executable(ftst_read_filter_topo_nml ftst_readnml.F90)

This maps in the program to the unit test. Other than the driver,
all ufs_utils program modules are compiled in a library.
  target_link_libraries(ftst_read_filter_topo_nml filter_topo_lib)

This adds the test so it can be run under CMake.
  add_test(NAME filter_topo-ftst_read_namelist COMMAND ftst_read_filter_topo_nml)

The above is for a serial test. To run an mpi-based test, use this command:
  add_mpi_test(filter_topo-ftst_read_namelist
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_filter_topo_nml
    NUMPROCS 4 TIMEOUT 60)
