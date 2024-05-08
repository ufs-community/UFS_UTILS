!> @file
!! @brief unit test for input parameters
!! @author Denise.Worthen@noaa.gov
!!
!! Test that readnml and readcsv will correctly fail if ocniceprep.nml, ocean.csv
!! or ice.csv contain non-valid values. The test passes if the correct number
!! of invalid conditions are indentified
!!
!! @author Denise.Worthen@noaa.gov
program ftst_program_setup

  use init_mod, only: readnml, readcsv

  implicit none

  character(len=120) :: testpath
  character(len=120) :: errmsg
  integer :: rc
  integer :: passed, ntests

  ntests = 4
  passed = 0
  testpath = 'tests/ocnice_prep/'

  call readnml(trim(testpath)//'data/input.nml',errmsg, rc)
  print *,trim(errmsg),rc
  if (rc .eq. 1)passed = passed+1

  call readnml(trim(testpath)//'data/invalid.model.nml', errmsg, rc)
  print *,trim(errmsg),rc
  if (rc .eq. 1)passed = passed+1

  call readnml(trim(testpath)//'data/invalid.srcdim.nml',errmsg, rc)
  print *,trim(errmsg),rc
  if (rc .eq. 1)passed = passed+1

  call readnml(trim(testpath)//'data/invalid.dstdim.nml',errmsg, rc)
  print *,trim(errmsg),rc
  if (rc .eq. 1)passed = passed+1

  if (passed .eq. ntests) then
     print '(2(a10,i6),a20)','SUCCESS! ',passed ,' tests out of ',ntests,' were identified '
  else
     print '(2(a10,i6),a20)','FAIL! only ',passed,' out of ',ntests,' were identified'
     stop 1
  endif

end program ftst_program_setup
