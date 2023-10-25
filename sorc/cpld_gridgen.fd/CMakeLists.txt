# This is the CMake build file for the chgres_cube utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    angles.F90
    charstrings.F90
    cicegrid.F90
    debugprint.F90
    gengrid_kinds.F90
    grdvars.F90
    inputnml.F90
    mapped_mask.F90
    postwgts.F90
    scripgrid.F90
    topoedits.F90
    tripolegrid.F90
    vartypedefs.F90
    vertices.F90)

set(exe_src gen_fixgrid.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
  
  # Turn on this argument mismatch flag for gfortran10.
# if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
# endif()
endif()

set(exe_name cpld_gridgen)

add_library(cpld_gridgen_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(cpld_gridgen_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(cpld_gridgen_lib INTERFACE ${mod_dir})

target_link_libraries(
  cpld_gridgen_lib
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE cpld_gridgen_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()
