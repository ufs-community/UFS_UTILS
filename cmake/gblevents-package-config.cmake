@PACKAGE_INIT@

# Include targets file.  This will create IMPORTED target @PROJECT_NAME@
include("${CMAKE_CURRENT_LIST_DIR}/gblevents-targets.cmake")
include(CMakeFindDependencyMacro)

find_dependency(nemsio CONFIG)
find_dependency(sigio CONFIG)
find_dependency(NetCDF COMPONENTS Fortran)

get_target_property(gblevents_BUILD_TYPES gblevents::gblevents IMPORTED_CONFIGURATIONS)

check_required_components("gblevents")

get_target_property(location gblevents::gblevents LOCATION)
message(STATUS "Found gblevents: ${location} (found version \"@PROJECT_VERSION@\")")
