# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Ed Hartnett

set(CPLDGRIDGEN_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/ufs_utils/unit_tests/cpld_gridgen")

set(GRID_FILE "ocean_hgrid.nc")

PULL_DATA(${CPLDGRIDGEN_URL} ${GRID_FILE})

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

add_executable(ftst_find_angq ftst_find_angq.F90)
add_test(NAME cpld_gridgen-ftst_find_angq COMMAND ftst_find_angq)
target_link_libraries(ftst_find_angq cpld_gridgen_lib)
