# This file contains CMake code related to MPI. (Taken from the
# ParallelIO project.)

# Jim Edwards
include (CMakeParseArguments)

# Find Valgrind to perform memory leak check
if (PIO_VALGRIND_CHECK)
  find_program (VALGRIND_COMMAND NAMES valgrind)
  if (VALGRIND_COMMAND)
    set (VALGRIND_COMMAND_OPTIONS --leak-check=full --show-reachable=yes)
  else ()
    message (WARNING "Valgrind not found: memory leak check could not be performed")
    set (VALGRIND_COMMAND "")
  endif ()
endif ()

#
# - Functions for parallel testing with CTest
#

#==============================================================================
# - Get the machine platform-specific
#
# Syntax:  platform_name (RETURN_VARIABLE)
#
function (platform_name RETURN_VARIABLE)

  # Determine platform name from site name...
  site_name (SITENAME)
  
  # ursa
  if (SITENAME MATCHES "^ufe01" OR
      SITENAME MATCHES "^ufe02" OR
      SITENAME MATCHES "^ufe03" OR
      SITENAME MATCHES "^ufe04")

    set (${RETURN_VARIABLE} "ursa" PARENT_SCOPE)

  # wcoss2
  elseif (SITENAME MATCHES "^along01" OR
      SITENAME MATCHES "^alogin02" OR
      SITENAME MATCHES "^clogin" OR
      SITENAME MATCHES "^dlogin")

    set (${RETURN_VARIABLE} "wcoss2" PARENT_SCOPE)

  # gaea c5
  elseif (SITENAME MATCHES "^gaea51" OR
      SITENAME MATCHES "^gaea52" OR
      SITENAME MATCHES "^gaea53" OR
      SITENAME MATCHES "^gaea54" OR
      SITENAME MATCHES "^gaea55" OR
      SITENAME MATCHES "^gaea56" OR
      SITENAME MATCHES "^gaea57" OR
      SITENAME MATCHES "^gaea58" OR
      SITENAME MATCHES "^gaea51.ncrc.gov" OR
      SITENAME MATCHES "^gaea52.ncrc.gov" OR
      SITENAME MATCHES "^gaea52.ncrc.gov" OR
      SITENAME MATCHES "^gaea54.ncrc.gov" OR
      SITENAME MATCHES "^gaea55.ncrc.gov" OR
      SITENAME MATCHES "^gaea56.ncrc.gov" OR
      SITENAME MATCHES "^gaea57.ncrc.gov" OR
      SITENAME MATCHES "^gaea58.ncrc.gov")

    set (${RETURN_VARIABLE} "gaeac5" PARENT_SCOPE)
  
  # gaea c6
  elseif (SITENAME MATCHES "^gaea61" OR
      SITENAME MATCHES "^gaea62" OR
      SITENAME MATCHES "^gaea63" OR
      SITENAME MATCHES "^gaea64" OR
      SITENAME MATCHES "^gaea65" OR
      SITENAME MATCHES "^gaea66" OR
      SITENAME MATCHES "^gaea67" OR
      SITENAME MATCHES "^gaea68" OR
      SITENAME MATCHES "^gaea61.ncrc.gov" OR
      SITENAME MATCHES "^gaea62.ncrc.gov" OR
      SITENAME MATCHES "^gaea62.ncrc.gov" OR
      SITENAME MATCHES "^gaea64.ncrc.gov" OR
      SITENAME MATCHES "^gaea65.ncrc.gov" OR
      SITENAME MATCHES "^gaea66.ncrc.gov" OR
      SITENAME MATCHES "^gaea67.ncrc.gov" OR
      SITENAME MATCHES "^gaea68.ncrc.gov")

    set (${RETURN_VARIABLE} "gaeac6" PARENT_SCOPE)

  elseif (SITENAME MATCHES "^Orion-login-1.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-2.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-3.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-4.HPC.MsState.Edu")

    set (${RETURN_VARIABLE} "orion" PARENT_SCOPE)

  elseif (SITENAME MATCHES "^Hercules-login-1.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Hercules-login-2.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Hercules-login-3.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Hercules-login-4.HPC.MsState.Edu" OR
      SITENAME MATCHES "^hercules-login-1.hpc.msstate.edu" OR
      SITENAME MATCHES "^hercules-login-2.hpc.msstate.edu" OR
      SITENAME MATCHES "^hercules-login-3.hpc.msstate.edu" OR
      SITENAME MATCHES "^hercules-login-4.hps.msstate.edu")

    set (${RETURN_VARIABLE} "hercules" PARENT_SCOPE)

  elseif (SITENAME MATCHES "^cheyenne1.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne1.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne2.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne3.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne4.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne5.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne6.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne1.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne2.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne3.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne4.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne5.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^cheyenne6.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin1.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin2.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin3.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin4.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin5.ib0.cheyenne.ucar.edu" OR
      SITENAME MATCHES "^chadmin6.ib0.cheyenne.ucar.edu")

    set (${RETURN_VARIABLE} "cheyenne" PARENT_SCOPE)
  elseif (SITENAME MATCHES "^login1.stampede2.tacc.utexas.edu" OR
      SITENAME MATCHES "^login2.stampede2.tacc.utexas.edu" OR
      SITENAME MATCHES "^login3.stampede2.tacc.utexas.edu" OR
      SITENAME MATCHES "^login4.stampede2.tacc.utexas.edu")
    

    set (${RETURN_VARIABLE} "stampede" PARENT_SCOPE)

  else ()

    set (${RETURN_VARIABLE} "unknown" PARENT_SCOPE)

  endif ()
endfunction ()

#==============================================================================
# - Add a new parallel test
#
# Syntax:  add_mpi_test (<TESTNAME>
#                        EXECUTABLE <command>
#                        ARGUMENTS <arg1> <arg2> ...
#                        NUMPROCS <num_procs>
#                        TIMEOUT <timeout>)
function (add_mpi_test TESTNAME)

  # Parse the input arguments
  set (options)
  set (oneValueArgs NUMPROCS TIMEOUT EXECUTABLE)
  set (multiValueArgs ARGUMENTS)
  cmake_parse_arguments (${TESTNAME} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Store parsed arguments for convenience
  set (exec_file ${${TESTNAME}_EXECUTABLE})
  set (exec_args ${${TESTNAME}_ARGUMENTS})
  set (num_procs ${${TESTNAME}_NUMPROCS})
  set (timeout ${${TESTNAME}_TIMEOUT})

  # Get the platform name
  platform_name (PLATFORM)

  get_property(WITH_MPIEXEC GLOBAL PROPERTY WITH_MPIEXEC)
  if (WITH_MPIEXEC)
    set(MPIEXEC "${WITH_MPIEXEC}")
  endif ()

  # Default ("unknown" platform) execution
  if (PLATFORM STREQUAL "unknown")

    # Run tests directly from the command line
    set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${num_procs}
      ${MPIEXEC_PREFLAGS} ${VALGRIND_COMMAND} ${VALGRIND_COMMAND_OPTIONS} ${exec_file}
      ${MPIEXEC_POSTFLAGS} ${exec_args})

  else ()

    # Run tests from the platform-specific executable
    set (EXE_CMD ${CMAKE_SOURCE_DIR}/cmake/mpiexec.${PLATFORM}
      ${num_procs} ${VALGRIND_COMMAND} ${VALGRIND_COMMAND_OPTIONS} ${exec_file} ${exec_args})

  endif ()

  # Add the test to CTest
  add_test(NAME ${TESTNAME} COMMAND ${EXE_CMD})

  # Adjust the test timeout
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${timeout})

endfunction()
