## HOW TO CREATE A UNIT TEST.

Unit tests should test only small parts of a program. For example,
a single routine or function.

Source code for a test shall be placed in a program
specific directory under ./tests.

The name of the unit test source file shall begin with
"ftst_" (for Fortran) or "tst_" (for C/C++).

A unit test shall have the following components.
For an example, see ./filter_topo/ftst_readnml.F90

- A prolog describing the test and its author.
- The test should begin with a print to standard
output, such as "Starting test of ..."
- A call to the function or routine to be tested.
- Output from the function or routine should be 
compared to expected values. If wrong values are
returned, the test should stop with a bad status.
- If the test is successful, that message
should be printed to standard output, i.e.,
print*,"SUCCESS!"

### HOW TO COMPILE THE TEST UNDER CMAKE. 

The above test is compiled and run with these
statements in the CMakeLists.txt file. For more
details, see ./filter_topo/CMakeLists.txt.

This statement names the unit test executable:
```
add_executable(ftst_read_filter_topo_nml ftst_readnml.F90)
```

This maps in the program library to the unit test. Other than the driver,
all ufs_utils program modules are compiled in a library.
  target_link_libraries(ftst_read_filter_topo_nml filter_topo_lib)

This adds the test so it can be run under CMake. Since this
test is part of the filter_topo program, note that the name 
begins with "filter_topo-". 

  add_test(NAME filter_topo-ftst_read_namelist COMMAND ftst_read_filter_topo_nml)

The above is for a serial test. To run an mpi-based test, use this command:
  add_mpi_test(filter_topo-ftst_read_namelist
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_filter_topo_nml
    NUMPROCS 4 TIMEOUT 60)

### TESTS WITH INPUT DATA:

Input data such as small text files may be placed in a ./data
sub-directory. The CMakeLists.txt file should then be updated
to copy this data to the run directory. For example, if
you test needs an input namelist (called input.nml), use the
following command:

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/input.nml ${CMAKE_CURRENT_BINARY_DIR}/input.nml)

Do not add large, binary input data to the ./data sub-directory. If
your test requires these data, contact the repository managers
for assistance.

### RUNNING THE TESTS:

Tests will automatically be run by Github actions upon each commit.
To run the tests locally, do a 'make test' after the 'make install'
command. For the parallel tests to run locally, update the
machine-specific "mpi_exec" script under UFS_UTILS/cmake with your
run account and queue.
