# This is the cmake build file for the tests directory of the
# UFS_UTILS project.
#
# George Gayno, Lin Gan, Ed Hartnett, Larissa Reames

set(CHGRES_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/UFS/ufs_utils/unit_tests/chgres_cube")

set(V16SFC_FILE "gfs.v16.sfc.history.nc")
set(V16SFC_FILE2 "gfs.v16.sfc2.history.nc")
set(V16ATM_FILE "gfs.v16.atm.history.nc")
set(V14SFC_FILE "gfs.v14.sfc.history.nemsio")
set(V15SFC_FILE "gfs.v15.sfc.history.nemsio")
set(GFS_GRIB_FILE "gfs.t00z.pgrb2.0p50.f000")

foreach(THE_FILE IN LISTS V16SFC_FILE V16SFC_FILE2 V16ATM_FILE V14SFC_FILE V15SFC_FILE GFS_GRIB_FILE) 
  PULL_DATA(${CHGRES_URL} ${THE_FILE})
endforeach()

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
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_fv3_tiled.nml ${CMAKE_CURRENT_BINARY_DIR}/data/config_fv3_tiled.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_fv3_tiled_warm_restart.nml ${CMAKE_CURRENT_BINARY_DIR}/data/config_fv3_tiled_warm_restart.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_gaussian_nemsio.nml ${CMAKE_CURRENT_BINARY_DIR}/data/config_gaussian_nemsio.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_gfs_grib2.nml ${CMAKE_CURRENT_BINARY_DIR}/data/config_gfs_grib2.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/global_hyblev.l28.txt ${CMAKE_CURRENT_BINARY_DIR}/data/global_hyblev.l28.txt)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_spectral_sigio.nml ${CMAKE_CURRENT_BINARY_DIR}/data/config_spectral_sigio.nml)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/GFSphys_varmap.txt ${CMAKE_CURRENT_BINARY_DIR}/data/GFSphys_varmap.txt)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/parm/msis_lib/msis21.parm ${CMAKE_CURRENT_BINARY_DIR}/data/msis21.parm)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/msis2.1_test_ref_dp.txt ${CMAKE_CURRENT_BINARY_DIR}/data/msis2.1_test_ref_dp.txt)
# This one does not end up in the data directory.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/data/config_gaussian_nemsio.nml ${CMAKE_CURRENT_BINARY_DIR}/fort.41)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/LSanSuppress.supp ${CMAKE_CURRENT_BINARY_DIR}/LSanSuppress.supp)

add_executable(ftst_utils ftst_utils.F90)
add_test(NAME chgres_cube-ftst_utils COMMAND ftst_utils)
target_link_libraries(ftst_utils chgres_cube_lib)

add_executable(ftst_read_atm_grib2 ftst_read_atm_grib2.F90)
target_link_libraries(ftst_read_atm_grib2 chgres_cube_lib)
add_mpi_test(chgres_cube-ftst_read_atm_grib2
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_atm_grib2
  NUMPROCS 1 TIMEOUT 60)

add_executable(ftst_read_sfc_grib2 ftst_read_sfc_grib2.F90)
target_link_libraries(ftst_read_sfc_grib2 chgres_cube_lib)
add_mpi_test(chgres_cube-ftst_read_sfc_grib2
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_sfc_grib2
  NUMPROCS 1 TIMEOUT 60)

add_executable(ftst_program_setup ftst_program_setup.F90)
target_link_libraries(ftst_program_setup chgres_cube_lib)
add_mpi_test(chgres_cube-ftst_program_setup
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_program_setup
  NUMPROCS 4 TIMEOUT 60)

add_executable(ftst_program_setup_varmaps ftst_program_setup_varmaps.F90)
target_link_libraries(ftst_program_setup_varmaps chgres_cube_lib)
add_mpi_test(chgres_cube-ftst_program_setup_varmaps
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_program_setup_varmaps
  NUMPROCS 4 TIMEOUT 60)

add_executable(ftst_search_util ftst_search_util.F90)
target_link_libraries(ftst_search_util chgres_cube_lib)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(ftst_search_util OpenMP::OpenMP_Fortran)
endif()

# Cause test to be run with MPI. The search routine was
# written to work on a single task. So set NUMPROCS to 1.
add_mpi_test(chgres_cube-ftst_search_util
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_search_util
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_sfc_input_data ftst_sfc_input_data.F90)
target_link_libraries(ftst_sfc_input_data chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_sfc_input_data
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_sfc_input_data
  NUMPROCS 4
  TIMEOUT 60)

add_executable(ftst_dint2p ftst_dint2p.F90)
target_link_libraries(ftst_dint2p
  chgres_cube_lib)

add_test(NAME chgres_cube-ftst_dint2p COMMAND ftst_dint2p)

add_executable(ftst_rh2spfh_gfs ftst_rh2spfh_gfs.F90)
target_link_libraries(ftst_rh2spfh_gfs
  chgres_cube_lib)

add_test(NAME chgres_cube-ftst_rh2spfh_gfs COMMAND ftst_rh2spfh_gfs)


add_executable(ftst_quicksort ftst_quicksort.F90)
target_link_libraries(ftst_quicksort
  chgres_cube_lib)

add_test(NAME chgres_cube-ftst_quicksort COMMAND ftst_quicksort)

add_executable(ftst_convert_winds ftst_convert_winds.F90)
target_link_libraries(ftst_convert_winds chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_convert_winds
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_convert_winds
  NUMPROCS 3
  TIMEOUT 60)

add_executable(ftst_read_sfc_gfs_nemsio ftst_read_sfc_gfs_nemsio.F90)
target_link_libraries(ftst_read_sfc_gfs_nemsio chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_sfc_gfs_nemsio
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_sfc_gfs_nemsio
  NUMPROCS 1
  TIMEOUT 60)

## Comment out this unit test until ESMF memory leaks are solved
## add_executable(ftst_surface_interp ftst_surface_interp.F90)
## target_link_libraries(ftst_surface_interp chgres_cube_lib)
##
## Cause test to be run with MPI.
## add_mpi_test(chgres_cube-ftst_surface_interp
##  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_surface_interp
##  NUMPROCS 1
##  TIMEOUT 60)

add_executable(ftst_read_sfc_nemsio ftst_read_sfc_nemsio.F90)
target_link_libraries(ftst_read_sfc_nemsio chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_sfc_nemsio
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_sfc_nemsio
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_read_sfc_netcdf ftst_read_sfc_netcdf.F90)
target_link_libraries(ftst_read_sfc_netcdf chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_sfc_netcdf
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_sfc_netcdf
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_read_nst_netcdf ftst_read_nst_netcdf.F90)
target_link_libraries(ftst_read_nst_netcdf chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_nst_netcdf
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_nst_netcdf
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_read_nst_nemsio ftst_read_nst_nemsio.F90)
target_link_libraries(ftst_read_nst_nemsio chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_nst_nemsio
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_nst_nemsio
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_read_atm_gaussian_netcdf ftst_read_atm_gaussian_netcdf.F90)
target_link_libraries(ftst_read_atm_gaussian_netcdf chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_read_atm_gaussian_netcdf
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_read_atm_gaussian_netcdf
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_surface_nst_landfill ftst_surface_nst_landfill.F90)
target_link_libraries(ftst_surface_nst_landfill chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_surface_nst_landfill
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_surface_nst_landfill
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_surface_regrid_many ftst_surface_regrid_many.F90)
target_link_libraries(ftst_surface_regrid_many chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_surface_regrid_many
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_surface_regrid_many
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_surface_search_many ftst_surface_search_many.F90)
target_link_libraries(ftst_surface_search_many chgres_cube_lib)

# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_surface_search_many
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_surface_search_many
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_read_vcoord ftst_read_vcoord.F90)
target_link_libraries(ftst_read_vcoord chgres_cube_lib)
add_test(NAME chgres_cube-ftst_read_vcoord COMMAND ftst_read_vcoord)

add_executable(ftst_msis_lib ftst_msis2.1_lib.F90)
target_link_libraries(ftst_msis_lib msis2)
# Cause test to be run with MPI.
add_mpi_test(chgres_cube-ftst_msis_lib
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ftst_msis_lib
  NUMPROCS 1
  TIMEOUT 60)

add_executable(ftst_example ftst_example.F90)
target_link_libraries(ftst_example chgres_cube_lib)
add_test(NAME chgres_cube-ftst_example COMMAND ftst_example)
