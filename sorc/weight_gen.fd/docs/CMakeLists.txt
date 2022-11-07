# This is the CMake file for building the docs directory of UFS_UTILS
# utility weight_gen
#
# Create doxyfile.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# Build documentation with target all.
add_custom_target(weight_gen_doc ALL
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API Documentation with Doxygen" VERBATIM)

# Ensure top-level docs have been generated.
add_dependencies(weight_gen_doc doc)
