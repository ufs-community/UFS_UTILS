# This is the CMake build file for the chgres_cube utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    arrays_mod.F90
    init_mod.F90
    ocncalc_mod.F90
    restarts_mod.F90
    utils_mod.F90
    utils_esmf_mod.F90)

set(exe_src ocniceprep.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fdefault-real-8")
endif()

set(exe_name oiprep)

add_library(ocnice_prep_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(ocnice_prep_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(ocnice_prep_lib INTERFACE ${mod_dir})

target_link_libraries(
  ocnice_prep_lib
  PUBLIC
  esmf
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(${exe_name} PRIVATE ocnice_prep_lib)

install(TARGETS ${exe_name})

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()
