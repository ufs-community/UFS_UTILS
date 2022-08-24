# This is the CMake file for the regional_esf_grid utility in
# UFS_UTILS.
#
# George Gayno
set(fortran_src
    pesg.f90
    pfun.f90
    pietc.f90
    pietc_s.f90
    pkind.f90
    pmat.f90
    pmat2.f90
    pmat4.f90
    pmat5.f90
    psym2.f90
    regional_esg_grid.f90)

set(exe_name regional_esg_grid)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(
  ${exe_name}
  NetCDF::NetCDF_Fortran)

install(TARGETS ${exe_name})
