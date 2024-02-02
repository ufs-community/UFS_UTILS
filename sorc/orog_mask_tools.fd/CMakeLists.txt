# This is the CMake build file for the orog_mask_tools in the
# UFS_UTILS package.
#
# George Gayno, Mark Potts

add_subdirectory(orog.fd)
add_subdirectory(orog_gsl.fd)
add_subdirectory(lake.fd)
add_subdirectory(inland.fd)
add_subdirectory(mask.fd)
add_subdirectory(topo.fd)
add_subdirectory(ramp.fd)

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  add_subdirectory(docs)  
endif()
