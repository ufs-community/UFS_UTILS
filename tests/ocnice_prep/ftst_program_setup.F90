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

  !integer, dimension(2) :: invalidsrcdims = (/2880, 2160/)
  !integer, dimension(2) :: invaliddstdims = (/1440, 1080/)

  !character(len=10) ::  invalidnml  = 'input.nml '
  !character(len=10) :: invalidfile  =  'atm      '
  !character(len=10) :: invalidangle = 'angchk   '
  !character(len=4), dimension(2) :: invaliduname = (/'vo  ', 'vvel'/)
  !character(len=4), dimension(2) :: invalidvname = (/'uo  ', 'uvel'/)

  !call readnml('data/input.nml',errmsg, rc)
  !print *,trim(errmsg),rc

  call readnml('data/invalid.model.nml', errmsg, rc)
  print *,trim(errmsg),rc

  !call readnml('data/invalid.srcdim.nml',errmsg, rc)
  !print *,trim(errmsg),rc

  !call readnml('data/invalid.dstdim.nml',errmsg, rc)
  !print *,trim(errmsg),rc
end program ftst_program_setup
