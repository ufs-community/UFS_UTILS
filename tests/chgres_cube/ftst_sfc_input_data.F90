 program test_sfc_input_data

! Test sfc_input_data. Tests include:
!
! subroutine soilt_check using a sample 3x3 matrix
!
! subroutine cnwat_check using a sample 3x3 matrix
!
! author: Larissa Reames (larissa.reames@noaa.gov)

 use mpi
 use esmf
 use utilities
 use model_grid

 implicit none

 integer :: ierr
 
 integer(esmf_kind_i4), allocatable :: mask(:,:)
 real(esmf_kind_r8), allocatable    :: soilt_bad(:,:,:), &
                                       soilt_updated(:,:,:), &
                                       soilt_correct(:,:,:)
 real(esmf_kind_r8), allocatable    :: skint(:,:)
 real(esmf_kind_r8), allocatable    :: cnwat_bad(:,:), &
                                       cnwat_updated(:,:), &
                                       cnwat_correct(:,:)

 call mpi_init(ierr)

!--------------------------------------------------------!
!------------------- CHECK_SOILT TEST -------------------!
!--------------------------------------------------------!

 lsoil_input = 2
 i_input = 3
 j_input = 3

 allocate(mask(i_input,j_input))
 allocate(skint(i_input,j_input))
 allocate(soilt_bad(i_input,j_input,lsoil_input))
 allocate(soilt_updated(i_input,j_input,lsoil_input))
 allocate(soilt_correct(i_input,j_input,lsoil_input))

!--------------------------------------------------------
! These variables are used to test all three functions of
! the check_soilt routine (i.e., replace water point 
! temperature and excessive land point soil temperature 
! with skin temperature and store a default ice column
! temperature in the soil temperature array at sea ice
! points since this array isn't available in grib2 files)
!--------------------------------------------------------

! Definition of the mask. The '0' indicates open water, '1' 
! any land point, and '2' sea ice.

 mask = reshape((/0, 1, 0,  0, 0, 1, 0, 0, 2/),(/3,3/))

! Soil temperature array with some values that are all incorrect
! based on landmask type. This will be the input soil array

 soilt_bad(1:i_input,1:j_input,1) = reshape((/0., 9999999.9, 0., &
                                              0., 0., 295.0, &
                                              0., 25.0, 25.0 /),(/3,3/))
 soilt_bad(1:i_input,1:j_input,2) = reshape((/0., 9999999.9, 0., &
                                              0., 0., 295.0, &
                                              0., 25.0, 25.0 /),(/3,3/))

! Subjective, reasonable skin temperature array.

 skint = reshape((/285.0, 295.0, 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 260.0/), &
                 (/3,3/))

!--------------------------------------------------------
! This is the corrected soil temperature array that 
! should be passed back from the check_soilt routine
!--------------------------------------------------------

 soilt_correct(1:i_input,1:j_input,1) = reshape( (/285.0, 295.0, 280.0, &
                                                   281.0, 282.0, 295.0, &
                                                   284.0, 285.0, 265.0/), (/3,3/))
 soilt_correct(1:i_input,1:j_input,2) = reshape( (/285.0, 295.0, 280.0, &
                                                   281.0, 282.0, 295.0, &
                                                   284.0, 285.0, 265.0/), (/3,3/))


 print*,"Starting test of check_soilt subroutine."

 soilt_updated = soilt_bad
 call check_soilt(soilt_updated,mask,skint)

 if (any(soilt_updated /= soilt_correct)) then
   print*,'SOILT TEST FAILED '
   print*,'SOILT ARRAY SHOULD BE:', soilt_correct
   print*,'SOILT ARRAY FROM TEST:', soilt_updated
   stop 2
 else
   print*, 'SOILT TEST SUCCEEDED. CONTINUING'
 endif
 deallocate(mask,skint,soilt_bad,soilt_updated,soilt_correct)

!--------------------------------------------------------!
!------------------- CHECK_CNWAT TEST -------------------!
!--------------------------------------------------------!

 i_input = 3
 j_input = 3

 allocate(cnwat_bad(i_input,j_input))
 allocate(cnwat_updated(i_input,j_input))
 allocate(cnwat_correct(i_input,j_input))

! Canopy water content array with some values that are too high. 
! This will be the input array for check_cnwat
 cnwat_bad = reshape((/999.0, 0.3, 0.4, 0.1, 0.3, 999.0, 999.0, 0.3, 0.15/), &
                 (/3,3/))

! Corrected canopy water content array that should be passed
! back from check_cnwat.
 cnwat_correct = reshape((/0.0, 0.3, 0.4, 0.1, 0.3, 0.0, 0.0, 0.3, 0.15/), &
                 (/3,3/))

 print*,"Starting test of check_cnwat subroutine."

 cnwat_updated = cnwat_bad
 call check_cnwat(cnwat_updated)

 if (any(cnwat_updated /= cnwat_correct)) then
   print*,'CNWAT TEST FAILED '
   print*,'CNWAT ARRAY SHOULD BE:', cnwat_correct
   print*,'CNWAT ARRAY FROM TEST:', cnwat_updated
   stop 2
 else
   print*, 'CNWAT TEST SUCCEEDED.'
 endif
 deallocate(cnwat_bad,cnwat_updated,cnwat_correct)

 call mpi_finalize(ierr)

 print*,"SUCCESS!"


 end program test_sfc_input_data
