# This is the CMake file for building this docs directory.
#
# Ed Hartnett 3/9/21

# Create doxyfile.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

# Build documentation with target all.
add_custom_target(ocean_merge_doc ALL
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API Documentation with Doxygen" VERBATIM)

# Ensure top-level docs have been generated.
add_dependencies(ocean_merge_doc doc)
