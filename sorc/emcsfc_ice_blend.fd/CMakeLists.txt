set(fortran_src
    emcsfc_ice_blend.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -i4")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O0")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-0")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name emcsfc_ice_blend)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE)
target_link_libraries(
  ${exe_name}
  bacio_4
  g2_4
  w3nco_4
  ${JASPER_LIBRARIES}
  PNG::PNG
  ZLIB::ZLIB)

install(TARGETS ${exe_name} RUNTIME DESTINATION exec)
