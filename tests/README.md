## HOW TO CREATE A UNIT TEST.

Unit tests should test only small parts of a program. For example,
a single routine or function, or multiple closely-linked routines.

Source code for a test shall be placed in a program
specific directory under ./tests.

The name of the unit test source file shall begin with
"ftst_" (for Fortran) or "tst_" (for C/C++).

A unit test shall have the following components.
For an example, see [ftst_readnml.F90](filter_topo/ftst_readnml.F90).

- A prologue describing the test and its author.
- The test should begin with a print to standard
output, such as "Starting test of ..."
- A call to the function or routine to be tested.
- Output from the function or routine must be 
compared to expected values. If wrong values are
returned, the test must stop with a non-zero exit code.
- If the test is successful, the test must return an
exit code of 0. And this message should be printed to
standard output: print*,"SUCCESS!"

**The test must never rely on human inspection of results - it should
be fully automated.**

### HOW TO COMPILE THE TEST UNDER CMAKE. 

The above test is compiled and run with these
statements in the CMakeLists.txt file. For more
details, see [CMakeLists.txt](filter_topo/CMakeLists.txt).

This statement names the unit test executable:
```
add_executable(ftst_read_filter_topo_nml ftst_readnml.F90)
```

This maps the program library to the unit test. Other than the driver,
all ufs_utils program modules are compiled in a library.
```
target_link_libraries(ftst_read_filter_topo_nml filter_topo_lib)
```

This adds the test so it can be run under CMake. Since this
test is part of the filter_topo program, note that the name 
begins with "filter_topo-". 
```
add_test(NAME filter_topo-ftst_read_namelist COMMAND ftst_read_filter_topo_nml)
```

The above is for a serial test. To run an mpi-based test, use this command:
```
add_mpi_test(filter_topo-ftst_read_namelist
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_filter_topo_nml
  NUMPROCS 4 TIMEOUT 60)
```

### TESTS WITH INPUT DATA:

Input data such as small text files may be placed in a ./data
sub-directory. The CMakeLists.txt file should then be updated
to copy this data to the run directory. For example, if
your test needs an input namelist (called input.nml), use the
following command:

```
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/input.nml ${CMAKE_CURRENT_BINARY_DIR}/input.nml)
```

**Do not** add large, binary input data to the ./data sub-directory. If
your test requires these data, contact the repository managers
for assistance.

### RUNNING THE TESTS:

Tests will automatically be run by Github actions upon each commit.

To run the tests locally, invoke one of the following after the `make install`:
- `make test` (Standard output for all tests sent to ./build/Testing/Temporary)
- `ctest --verbose` (Standard output for all tests sent to the screen)
- `ctest -R test_name` (Standard output for one test sent to the screen)

For the parallel tests to run locally, update the machine-specific
"mpi_exec" script under [cmake](../cmake) with your run account and queue.

### EXAMPLE UNIT TEST:

An example unit test, ftst_example.F90, exists in the tests/chgres_cube directory
and is currently compiled and run with existing tests

This simple test checks whether the routine rh2spfh in chgres_cube.fd/grib2_util.F90
is working correctly, and contains detailed comments explaining what each section 
of the test does. It also prompts the user to create additional lines of code to 
test one more subroutine from grib2_util.F90. Use this to get started understanding
the unit test framework.

### QUESTIONS

Please contact the repository managers: https://github.com/ufs-community/UFS_UTILS/wiki
