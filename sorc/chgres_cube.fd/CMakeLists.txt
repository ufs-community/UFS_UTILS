# This is the CMake build file for the chgres_cube utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    atmosphere.F90
    grib2_util.F90
    input_data.F90
    model_grid.F90
    program_setup.F90
    search_util.F90
    static_data.F90
    surface.F90
    thompson_mp_climo_data.F90
    utils.F90
    write_data.F90)

set(exe_src chgres.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
  
  # Turn on this argument mismatch flag for gfortran10.
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(exe_name chgres_cube)

add_library(chgres_cube_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(chgres_cube_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(chgres_cube_lib INTERFACE ${mod_dir})

target_link_libraries(
  chgres_cube_lib
  PUBLIC
  nemsio::nemsio
  sfcio::sfcio
  sigio::sigio
  bacio::bacio_4
  sp::sp_d
  w3nco::w3nco_d
  esmf
  wgrib2::wgrib2_lib
  wgrib2::wgrib2_api
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE chgres_cube_lib)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()

