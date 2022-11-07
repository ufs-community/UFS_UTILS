# This is the cmake build file for the fvcom_tools test directory of the
# UFS_UTILS project.
#
# David Wright

set(FVCOM_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/ufs_utils/unit_tests/fvcom_tools")

set(FVCOM_FILE "fvcom_unittest.nc")
set(SFCDATA_FILE "sfcdata_unittest.nc")

foreach(THE_FILE IN LISTS FVCOM_FILE SFCDATA_FILE)
  PULL_DATA(${FVCOM_URL} ${THE_FILE})
endforeach()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

include_directories(${PROJECT_SOURCE_DIR})

# Copy necessary test files from the source data directory to the
# build data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

add_executable(ftst_readfvcomnetcdf ftst_readfvcomnetcdf.F90)
add_test(NAME fvcom_tools-ftst_readfvcomnetcdf COMMAND ftst_readfvcomnetcdf)
target_link_libraries(ftst_readfvcomnetcdf fvcom_tools_lib)
