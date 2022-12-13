# This is the CMake file for the tests directory of the UFS_UTILS
# project.
#
# Ed Hartnett 2/11/21

# Some test files are large and are kept on the NOAA EMC FTP
# site. This function is used to download such test data. It takes two
# arguments, the URL and the file to be downloaded.
function(PULL_DATA THE_URL THE_FILE)
  # If the TEST_FILE_DIR was specified, look for our test data files
  # there before FTPing them. Developers can keep all test files on
  # their machines, and save the time of downloading them every time.
  message(STATUS "TEST_FILE_DIR ${TEST_FILE_DIR} THE_FILE ${THE_FILE}")  
  if(NOT ${TEST_FILE_DIR} STREQUAL ".")
    if (EXISTS ${TEST_FILE_DIR}/${THE_FILE})
      message(STATUS "Copying file ${TEST_FILE_DIR}/${THE_FILE} to test data directory.")
      FILE(COPY ${TEST_FILE_DIR}/${THE_FILE}
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
    endif()
  endif()
  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/data/${THE_FILE}")
    message(STATUS "Downloading file ${CMAKE_CURRENT_BINARY_DIR}/data/${THE_FILE}.")
    file(DOWNLOAD
      ${THE_URL}/${THE_FILE}
      ${CMAKE_CURRENT_BINARY_DIR}/data/${THE_FILE}
      SHOW_PROGRESS
      STATUS status
      INACTIVITY_TIMEOUT 30
    )
    list(GET status 0 status_num)
    if(NOT status_num EQUAL 0 OR NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/data/${THE_FILE}")
      message(FATAL_ERROR "Could not download ${THE_FILE}")
    endif()
  endif()
endfunction()

# Add the test subdirecotries.
# fvcom test only works for Intel. Comment out for now.
#add_subdirectory(fvcom_tools)
add_subdirectory(filter_topo)
add_subdirectory(chgres_cube)
add_subdirectory(fre-nctools)
add_subdirectory(global_cycle)
add_subdirectory(sfc_climo_gen)
add_subdirectory(cpld_gridgen)
add_subdirectory(emcsfc_snow2mdl)
