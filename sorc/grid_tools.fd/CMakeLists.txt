# This is the CMake file for the grid_tools utilites of UFS_UTILS.
#
# George Gayno
add_subdirectory(shave.fd)
add_subdirectory(filter_topo.fd)
add_subdirectory(regional_esg_grid.fd)
add_subdirectory(global_equiv_resol.fd)

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
