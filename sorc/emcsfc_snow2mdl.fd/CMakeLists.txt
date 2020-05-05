set(fortran_src
    driver.f
    grib_utils.f
    model_grid.f
    program_setup.f
    snow2mdl.f
    snowdat.f)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -r8 -i4 -FR -convert big_endian -assume byterecl")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O0")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-form -ffree-line-length-0 -fdefault-real-8 -fconvert=big-endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check"
  )
endif()

set(exe_name emcsfc_snow2mdl)
add_executable(${exe_name} ${fortran_src})
target_include_directories(
  ${exe_name} PRIVATE)
target_link_libraries(
  ${exe_name}
  g2_d
  ip_d
  sp_d
  landsfcutil_d
  bacio_4
  w3nco_d
  ${JASPER_LIBRARIES}
  PNG::PNG
  ZLIB::ZLIB
  OpenMP::OpenMP_Fortran)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
