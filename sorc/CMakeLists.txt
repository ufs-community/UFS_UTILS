# This is the cmake build file for the sorc directory in the
# UFS_UTILS project.
#
# George Gayno
add_subdirectory(emcsfc_ice_blend.fd)
add_subdirectory(emcsfc_snow2mdl.fd)
if (OpenMP_FOUND)
  add_subdirectory(global_cycle.fd)
  add_subdirectory(lsm_routines.fd)
else()
  message(STATUS "OpenMP is required for global_cycle.fd and was NOT found, skipping ...")
endif()
add_subdirectory(fre-nctools.fd)
add_subdirectory(grid_tools.fd)
add_subdirectory(chgres_cube.fd)
add_subdirectory(chgres_tools.fd)
add_subdirectory(chgres_tools2.fd)
add_subdirectory(orog_mask_tools.fd)
add_subdirectory(sfc_climo_gen.fd)
add_subdirectory(vcoord_gen.fd)
add_subdirectory(fvcom_tools.fd)
add_subdirectory(gblevents.fd)
add_subdirectory(cpld_gridgen.fd)
add_subdirectory(weight_gen.fd)
