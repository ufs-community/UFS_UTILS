!> @file
!! @brief unit test for input parameters
!! @author Denise.Worthen@noaa.gov
!!
!! Test that readnml and readcsv will correctly fail if ocniceprep.nml, ocean.csv
!! or ice.csv contain non-valid values
!!
!! @author Denise.Worthen@noaa.gov
program ftst_program_setup

  use init_mod, only: readnml, readcsv

  implicit none

  character(len=120) :: errmsg
  integer :: rc

  !call readnml('data/input.nml',errmsg, rc)
  !print *,trim(errmsg),rc

  call readnml('data/invalid.model.nml', errmsg, rc)
  print *,trim(errmsg),rc

  !call readnml('data/invalid.srcdim.nml',errmsg, rc)
  !print *,trim(errmsg),rc

  !call readnml('data/invalid.dstdim.nml',errmsg, rc)
  !print *,trim(errmsg),rc
end program ftst_program_setup
