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
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

if(CHGRES_ALL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCHGRES_ALL")
endif()

set(exe_name chgres_cube)

add_library(chgres_cube_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(chgres_cube_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(chgres_cube_lib INTERFACE ${mod_dir})

if(CHGRES_ALL)
  target_link_libraries(
    chgres_cube_lib 
    PUBLIC
    nemsio::nemsio
    sfcio::sfcio
    sigio::sigio)
endif()

target_link_libraries(
  chgres_cube_lib
  PUBLIC
  g2::g2_d
  bacio::bacio_4
  ip::ip_d
  w3emc::w3emc_d
  esmf
  msis2
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(sp_FOUND)
  target_link_libraries(chgres_cube_lib PUBLIC sp::sp_d)
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE chgres_cube_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()

