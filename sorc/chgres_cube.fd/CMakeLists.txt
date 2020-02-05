
file(GLOB fortran_src
  ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/*.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/*.F)

if(IntelComp)
  set(fortran_flags "-g" "-O3" "-axCORE-AVX2")
elseif(GNUComp)
  set(fortran_flags "-ffree-form" "-ffree-line-length-0" "-g" "-O3")
else()
  message("unknown compiler!")
  exit()
endif()

set(exe_name chgres_cube.exe)
add_executable(${exe_name} ${fortran_src})
target_compile_options(${exe_name} PRIVATE "${fortran_flags}")

target_include_directories(${exe_name} PRIVATE ${ESMF_INC}
  ${NETCDF_INCLUDES}
  ${MPI_Fortran_INCLUDE_PATH}
  ${NETCDF_INCLUDES_F90})


set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include_4")
set_target_properties(${exe_name} PROPERTIES Fortran_MODULE_DIRECTORY ${module_dir})

target_include_directories(${exe_name} PUBLIC
  $<BUILD_INTERFACE:${module_dir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include_4>)

target_include_directories(${exe_name}
  PRIVATE
  ${WGRIB2_INCLUDES}
  ${NETCDF_INCLUDES}
  ${NETCDF_INCLUDES_F90}
  ${MPI_Fortran_INCLUDE_PATH})

target_link_libraries(${exe_name}
  ${WGRIB2_LIBRARIES}
  ${NETCDF_LIBRARIES_F90} ${NETCDF_LIBRARIES}
  nemsio
  sfcio_4 ${MPI_Fortran_LIBRARIES}
  sigio_4
  bacio_4
  sp_d
  w3nco_d
  ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_LIBRARIES}
  ${OpenMP_Fortran_FLAGS}
  esmf) 

install(TARGETS ${exe_name}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib) 
install(DIRECTORY ${module_dir} DESTINATION include_4)

