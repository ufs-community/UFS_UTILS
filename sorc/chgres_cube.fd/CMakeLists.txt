# This is the CMake build file for the chgres_cube utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser, Larissa Reames

set(lib_src
    atmosphere.F90
    atmosphere_target_data.F90
    grib2_util.F90
    atm_input_data.F90
    sfc_input_data.F90
    nst_input_data.F90
    model_grid.F90
    program_setup.F90
    search_util.F90
    static_data.F90
    surface_target_data.F90
    surface.F90
    thompson_mp_climo_data.F90
    wam_climo_data.f90
    utils.F90
    write_data.F90)

set(exe_src chgres.F90)

add_subdirectory(msis2.1.fd)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")

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
  g2::g2_d
  nemsio::nemsio
  sfcio::sfcio
  sigio::sigio
  bacio::bacio_4
  ip::ip_d
  sp::sp_d
  w3emc::w3emc_d
  esmf
  msis2
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE chgres_cube_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()

