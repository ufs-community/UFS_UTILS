!> @file
!! @brief unit test for QC of source CICE restarts
!! @author Denise.Worthen@noaa.gov
!!
!! Test that land values and non-physical values will be found and
!! removed and that good values are not removed
!! @author Denise.Worthen@noaa.gov
program ftst_qcice

  use utils_mod , only : zero_out_land_ice, zero_out_phantom_ice

  implicit none

  integer, parameter :: ncat=3, nxy = 5

  integer      :: kmt(nxy)
  real(kind=8) :: field1(ncat,nxy), field2(ncat,nxy)

  integer :: passed, ntests_land_ice, ntests_phantom_ice
  integer :: expected_cnt, icnt

  ntests_phantom_ice = 0
  ntests_land_ice = 0
  passed = 0

  ! ------------------------------------
  ! Tests for zero_out_land_ice
  ! ------------------------------------

  ! test 1: land only, values everywhere
  ntests_land_ice = ntests_land_ice + 1
  kmt = 0
  field1(:,:) = 1.0e-3
  expected_cnt = nxy
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test2: ocean only, values everywhere
  ntests_land_ice = ntests_land_ice + 1
  kmt = 1
  field1(:,:) = 1.0e-3
  expected_cnt = 0
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test3: mix land and ocean, values everywhere
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 1.0e-3
  expected_cnt = 3
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test4: mix land and ocean, values nowhere
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 0.0
  expected_cnt = 0
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test5: mix land and ocean, values at one category/location only at land
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 0.0
  field1(ncat,5) = 1.0e-3
  expected_cnt = 1
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test6: mix land and ocean, values at all categories/single location only at land
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 0.0
  field1(:,1) = 1.0e-3
  expected_cnt = 1
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test7: mix land and ocean, values at one category/location only at ocean
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 0.0
  field1(ncat,3) = 1.0e-3
  expected_cnt = 0
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  ! test8: mix land and ocean, values at all categories/single location only at ocean
  ntests_land_ice = ntests_land_ice + 1
  kmt = (/0,1,1,0,0/)
  field1(:,:) = 0.0
  field1(:,3) = 1.0e-3
  expected_cnt = 0
  call zero_out_land_ice(kmt, field1, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,expected_cnt,icnt,passed

  if (passed .eq. ntests_land_ice) then
     print '(2(a36,i6),a20)','land-ice QC test: SUCCESS! ',passed ,' tests out of ',ntests_land_ice
  else
     print '(2(a36,i6),a20)','land-ice QC test: FAIL! only ',passed,' out of ',ntests_land_ice
     print '(a)',' Phantom ice QC tests will not run '
     stop 1
  endif

  passed = 0
  ! ------------------------------------
  ! Tests for zero_out_phantom_ice
  ! ------------------------------------

  ! test1: all zeros, field1 and 2
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,:) = 0.0
  field2(:,:) = 0.0
  expected_cnt = 0
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test1: ',expected_cnt,icnt,passed

  ! test2: all values, field1 and 2
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,:) = 2.0e-1
  field2(:,:) = 1.0e-3
  expected_cnt = 0
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test2: ',expected_cnt,icnt,passed

  ! test3: field1 contains values, field2 all zero
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,:) = 1.e-3
  field2(:,:) = 0.0
  ! the expected_cnt = 0 because the QC does not check for this case
  expected_cnt = 0
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test3: ',expected_cnt,icnt,passed

  ! test4: field1 is zero, field2 contains values
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,:) = 0.0
  field2(:,:) = 1.e-3
  expected_cnt = nxy
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test4: ',expected_cnt,icnt,passed

  ! test5: field1 has values single cat & location, field2 is zero
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(2,1) = 1.e-3
  field2(:,:) = 0.0
  ! the expected_cnt = 0 because the QC does not check for this case
  expected_cnt = 0
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test5: ',expected_cnt,icnt,passed

  ! test6: field2 has values single cat & location, field1 is zero
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,:) = 0.0
  field2(ncat,3) = 1.e-3
  expected_cnt = 1
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test6: ',expected_cnt,icnt,passed

  ! test7: field1 has values all cat & 1 location, field2 is zero
  ntests_phantom_ice = ntests_phantom_ice + 1
  field1(:,3) = 1.e-3
  field2(:,:) = 0.0
  ! the expected_cnt = 0 because the QC does not check for this case
  expected_cnt = 0
  call zero_out_phantom_ice(field1, field2, icnt)
  if (icnt == expected_cnt) passed = passed + 1
  !print *,'test7: ',expected_cnt,icnt,passed

  if (passed .eq. ntests_phantom_ice) then
     print '(2(a36,i6),a20)','phantom-ice QC test: SUCCESS! ',passed ,' tests out of ',ntests_phantom_ice
  else
     print '(2(a36,i6),a20)','phantom-ice QC test: FAIL! only ',passed,' out of ',ntests_phantom_ice
     stop 1
  endif

end program ftst_qcice
