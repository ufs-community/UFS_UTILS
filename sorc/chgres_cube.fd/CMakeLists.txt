cmake_minimum_required(VERSION 2.6)
set(EXENAME  chgres_cube.exe )
file(GLOB EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.F90 ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.F)
#file(GLOB EXE_C_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
if(IntelComp)
    message("setting intel flags")
    set(CMAKE_Fortran_4_FLAGS " -g -O3 -axCORE-AVX2" CACHE INTERNAL "" )
else()
    message("unkown compiler!")
    set(CMAKE_Fortran_4_FLAGS "-ffree-form -ffree-line-length-0 -g -O3 " CACHE INTERNAL "" )
endif()


#set_source_files_properties( ${EXE_C_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS} )
#add_definitions( ${CMAKE_Fortran_4_FLAGS} )
add_executable(${EXENAME} ${EXE_SRC} ${EXE_C_SRC} )
set_target_properties( ${EXENAME} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_4_FLAGS} )
set_target_properties( ${EXENAME} PROPERTIES LINK_FLAGS ${CMAKE_Fortran_4_FLAGS} )
message("sfcio lib is ${SFCIO_LIB} ${W3NCO_LIBd}")
include_directories( ${G2_INC4} ${ESMF_INC} ${NEMSIO_INC} ${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_PATH} ${SFCIO_INC} ${NETCDF_INCLUDES_F90})
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include_4")
set_target_properties(${EXENAME} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include_4")
target_include_directories(${EXENAME} PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include_4")
target_link_libraries(${EXENAME}  ${NETCDF_LIBRARIES_F90} ${NETCDF_LIBRARIES} ${NEMSIO_LIB} ${ESMF_LIB} ${SFCIO_LIB} ${MPI_Fortran_LIBRARIES} ${SIGIO_LIB} ${BACIO_LIB4} ${SP_LIBd} ${W3NCO_LIBd} ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_LIBRARIES} ${OpenMP_Fortran_FLAGS} )

file(MAKE_DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

install(TARGETS ${EXENAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
install(DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY}/ DESTINATION include_4 )

