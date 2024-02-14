# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett, Larissa Reames

set(CYCLE_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/ufs_utils/unit_tests/global_cycle")

set(FILE1 "soil_xainc.001")
set(FILE2 "soil_xainc.002")
set(FILE3 "soil_xainc.003")
set(FILE4 "soil_xainc.004")
set(FILE5 "soil_xainc.005")
set(FILE6 "soil_xainc.006")

foreach(THE_FILE IN LISTS FILE1 FILE2 FILE3 FILE4 FILE5 FILE6)
  PULL_DATA(${CYCLE_URL} ${THE_FILE})
endforeach()

# Include cmake to allow parallel I/O tests.
include (LibMPI)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

include_directories(${PROJECT_SOURCE_DIR})

add_executable(ftst_land_increments ftst_land_increments.F90)
target_link_libraries(ftst_land_increments global_cycle_lib)
add_executable(ftst_read_increments ftst_read_increments.F90)
target_link_libraries(ftst_read_increments global_cycle_lib)

add_mpi_test(global_cycle-ftst_land_increments
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_land_increments
  NUMPROCS 1 TIMEOUT 60)
add_mpi_test(global_cycle-ftst_read_increments
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_increments
  NUMPROCS 6 TIMEOUT 60)
