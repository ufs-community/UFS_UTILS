# This is the CMake build file for the vcoord_gen utility in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts, Kyle Gerheiser

set(lib_src
    driver.f90
    matrix_utils.f90
    vcoord_gen.f90)

set(exe_src driver.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

set(exe_name vcoord_gen)

add_library(vcoord_gen_lib STATIC ${lib_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(vcoord_gen_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(vcoord_gen_lib INTERFACE ${mod_dir})

target_link_libraries(${exe_name} vcoord_gen_lib)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})


# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
