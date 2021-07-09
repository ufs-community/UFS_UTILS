# This is the cmake build file for the lsm_routines directory in the
# UFS_UTILS project.
#
# Clara Draper

add_subdirectory(noah.fd)

if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()
