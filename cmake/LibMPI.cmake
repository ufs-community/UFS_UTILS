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
  
  # hera
  if (SITENAME MATCHES "^hfe01" OR
      SITENAME MATCHES "^hfe02" OR
      SITENAME MATCHES "^hfe03" OR
      SITENAME MATCHES "^hfe04" OR
      SITENAME MATCHES "^hfe05" OR
      SITENAME MATCHES "^hfe06" OR
      SITENAME MATCHES "^hfe07" OR
      SITENAME MATCHES "^hfe08" OR
      SITENAME MATCHES "^hfe09" OR
      SITENAME MATCHES "^hfe10" OR
      SITENAME MATCHES "^hfe11" OR
      SITENAME MATCHES "^hfe12")

    set (${RETURN_VARIABLE} "hera" PARENT_SCOPE)

  # wcoss_cray (Luna)
  elseif (SITENAME MATCHES "^llogin1" OR
      SITENAME MATCHES "^llogin2" OR
      SITENAME MATCHES "^llogin3")

    set (${RETURN_VARIABLE} "wcoss_cray" PARENT_SCOPE)
    
  # wcoss_cray (Surge)
  elseif (SITENAME MATCHES "^slogin1" OR
      SITENAME MATCHES "^slogin2" OR
      SITENAME MATCHES "^slogin3")

    set (${RETURN_VARIABLE} "wcoss_cray" PARENT_SCOPE)
    
  # wcoss_dell_p3 (Venus)
  elseif (SITENAME MATCHES "^v71a1.ncep.noaa.gov" OR
      SITENAME MATCHES "^v71a2.ncep.noaa.gov" OR
      SITENAME MATCHES "^v71a3.ncep.noaa.gov" OR
      SITENAME MATCHES "^v72a1.ncep.noaa.gov" OR
      SITENAME MATCHES "^v72a2.ncep.noaa.gov" OR
      SITENAME MATCHES "^v72a3.ncep.noaa.gov")

    set (${RETURN_VARIABLE} "wcoss_dell_p3" PARENT_SCOPE)
  
  # wcoss_dell_p3 (Mars)
  elseif (SITENAME MATCHES "^m71a1.ncep.noaa.gov" OR
      SITENAME MATCHES "^m71a2.ncep.noaa.gov" OR
      SITENAME MATCHES "^m71a3.ncep.noaa.gov" OR
      SITENAME MATCHES "^m72a1.ncep.noaa.gov" OR
      SITENAME MATCHES "^m72a2.ncep.noaa.gov" OR
      SITENAME MATCHES "^m72a3.ncep.noaa.gov")

    set (${RETURN_VARIABLE} "wcoss_dell_p3" PARENT_SCOPE)

  # wcoss2
  elseif (SITENAME MATCHES "^along01" OR
      SITENAME MATCHES "^alogin02")

    set (${RETURN_VARIABLE} "wcoss2" PARENT_SCOPE)

  # gaea
  elseif (SITENAME MATCHES "^gaea9" OR
      SITENAME MATCHES "^gaea10" OR
      SITENAME MATCHES "^gaea11" OR
      SITENAME MATCHES "^gaea12" OR
      SITENAME MATCHES "^gaea13" OR
      SITENAME MATCHES "^gaea14" OR
      SITENAME MATCHES "^gaea15" OR
      SITENAME MATCHES "^gaea16" OR
      SITENAME MATCHES "^gaea9.ncrc.gov" OR
      SITENAME MATCHES "^gaea10.ncrc.gov" OR
      SITENAME MATCHES "^gaea11.ncrc.gov" OR
      SITENAME MATCHES "^gaea12.ncrc.gov" OR
      SITENAME MATCHES "^gaea13.ncrc.gov" OR
      SITENAME MATCHES "^gaea14.ncrc.gov" OR
      SITENAME MATCHES "^gaea15.ncrc.gov" OR
      SITENAME MATCHES "^gaea16.ncrc.gov")

    set (${RETURN_VARIABLE} "gaea" PARENT_SCOPE)
    
  # jet 
  elseif (SITENAME MATCHES "^fe1" OR
      SITENAME MATCHES "^fe2" OR
      SITENAME MATCHES "^fe3" OR
      SITENAME MATCHES "^fe4" OR
      SITENAME MATCHES "^fe5" OR
      SITENAME MATCHES "^fe6" OR
      SITENAME MATCHES "^fe7" OR
      SITENAME MATCHES "^fe8" OR
      SITENAME MATCHES "^tfe1" OR
      SITENAME MATCHES "^tfe2")

    set (${RETURN_VARIABLE} "jet" PARENT_SCOPE)

  elseif (SITENAME MATCHES "^Orion-login-1.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-2.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-3.HPC.MsState.Edu" OR
      SITENAME MATCHES "^Orion-login-4.HPC.MsState.Edu")

    set (${RETURN_VARIABLE} "orion" PARENT_SCOPE)

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

  elseif (SITENAME MATCHES "^s4-submit.ssec.wisc.edu")

    set (${RETURN_VARIABLE} "s4" PARENT_SCOPE)

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
