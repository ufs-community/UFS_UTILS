# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett, Larissa Reames

set(SFCGEN_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/ufs_utils/unit_tests/sfc_climo_gen")
set(MOSAIC_FILE "C450_mosaic.nc")
set(GRID_FILE "C450_grid.tile7.nc")
set(ORO_FILE "C450_oro_data.tile7.nc")

foreach(THE_FILE IN LISTS MOSAIC_FILE GRID_FILE ORO_FILE)
  PULL_DATA(${SFCGEN_URL} ${THE_FILE}) 
endforeach()

set(SFCGEN2_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/GFS/fix/fix_sfc_climo")
set(VEG_FILE "vegetation_type.viirs.igbp.0.1.nc")

PULL_DATA(${SFCGEN2_URL} ${VEG_FILE})

# Include cmake to allow parallel I/O tests.
include (LibMPI)

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

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/fort.41 ${CMAKE_CURRENT_BINARY_DIR}/fort.41)

add_executable(ftst_interp ftst_interp.F90)
target_link_libraries(ftst_interp sfc_climo_gen_lib)

add_executable(ftst_pgm_setup ftst_program_setup.F90)
target_link_libraries(ftst_pgm_setup sfc_climo_gen_lib)

add_mpi_test(sfc_climo_gen-ftst_pgm_setup
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_pgm_setup
  NUMPROCS 1
  TIMEOUT 60)

add_mpi_test(sfc_climo_gen-ftst_interp
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_interp
  NUMPROCS 1
  TIMEOUT 60)
