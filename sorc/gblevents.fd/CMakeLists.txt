
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

add_library(gblevents STATIC gblevents.f)

target_link_libraries(gblevents PRIVATE nemsio::nemsio sigio::sigio
  NetCDF::NetCDF_Fortran)

install(
  TARGETS gblevents
  EXPORT gblevents-exports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/gblevents)

export(EXPORT gblevents-exports
  NAMESPACE gblevents::
  FILE gblevents-targets.cmake)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/gblevents-package-config.cmake ${CMAKE_BINARY_DIR}/gblevents-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(FILES ${CMAKE_BINARY_DIR}/gblevents-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/gblevents-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(EXPORT gblevents-exports
  NAMESPACE gblevents::
  FILE gblevents-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
