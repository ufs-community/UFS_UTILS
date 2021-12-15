! Unit test for rh2spfh in grib2_utils to be used as an example for users
! learning how to write unit tests. Users are prompted to add an additional
! test for convert_omega
! Larissa Reames OU/CIMMS/NOAA/NSSL/FRDD
 
program ftst_example

 use esmf
 use model_grid, only : i_input, j_input 
 use grib2_util, only : rh2spfh
 
 implicit none

 real(esmf_kind_r4), allocatable :: rh_spfh(:,:), &
                                    rh_orig(:,:), &
                                    spfh_returned(:,:), &
                                    spfh_correct(:,:)
 real(esmf_kind_r8),allocatable  :: t(:,:)
 real(esmf_kind_r8) :: p
 real,parameter          :: EPS = 1.0E-6

 i_input = 2
 j_input = 2
 allocate(rh_spfh(i_input,j_input))
 allocate(rh_orig(i_input,j_input))
 allocate(spfh_returned(i_input,j_input))
 allocate(spfh_correct(i_input,j_input))
 allocate(t(i_input,j_input))

 ! -------------------------------------------------------------------------
 ! Set constants/values to be passed to the unit test. In this case it's a
 ! set of single values, but it could be more complicated, like an 
 ! n-dimensional array or ESMF objects.
 ! -------------------------------------------------------------------------

 rh_spfh(:,:) = 50.0  ! Relative humidity (%)
 p = 100000.0    ! Pressure (Pa)
 t(:,:) = 300.0       ! Temperature (K)
 spfh_correct(:,:) = 10.978297E-3 ! Correct specific humidity value (kg/kg)

 print*, "Starting Unit Testing rh2spfh."

 !-------------------------------------------------------------------------
 ! Execute testing below by running target function rh2spfh and providing
 ! values set above
 !-------------------------------------------------------------------------

 rh_orig = rh_spfh          ! Save the original RH value for posterity
 call rh2spfh(rh_spfh,p,t) 
 spfh_returned = rh_spfh    ! Rename the returned value for clarity

 !-------------------------------------------------------------------------
 ! Check the returned value against what we know to be the correct answer. 
 ! If the correct result is returned (within a certain small tolerance), 
 ! then the test passes and the called routine is working as expected. If the 
 ! incorrect value is passed back, the test fails and an error is returned.
 !-------------------------------------------------------------------------

 if ( any(abs(spfh_returned - spfh_correct) .gt. EPS)) then 
   stop 1
   print*, "RH2SPFH TEST FAILED."
 endif

 !-------------------------------------------------------------------------
 ! If you are trying to debug a test failure, code like the commented
 ! section below might prove useful.
 !-------------------------------------------------------------------------

 ! if ( any(abs(spfh_returned - spfh_correct) .lt. EPS)) then
 !   print*, "RH2SPFH TEST PASSED. SUCCESS!"
 ! else
 !   print*, "RH2SPFH TEST FAILED."
 !   print*, "TEST RETURNED VALUE OF ", spfh_returned
 !   print*, "RETURNED VALUE EXPECT TO BE ", spfh_correct
 !   stop 1
 ! endif

 !-------------------------------------------------------------------------
 ! Make sure to deallocate any and all allocatable arrays that you use
 !-------------------------------------------------------------------------

 deallocate(rh_spfh,spfh_correct,rh_orig,spfh_returned,t)

 ! -------------------------------------------------------------------------
 ! You can test multiple subroutines (units) in any test file. This would
 ! be a good place to test the other subroutine in grib2_util, 
 ! convert_omega. Make note of the difference in variable size for this 
 ! routine. You don't have to pass an array of size you'd normally 
 ! encounter for these variables (like 200x200x64), just choose a small
 ! size with the proper dimensionality, say 3x3x2, and fill it with values
 ! typical of the various arrays. You can check the returned array element-
 ! by-element, or use the any() command to check all elements at once. Make
 ! sure to provide a helpful failure message indicating where the failure 
 ! occured and what the returned/expected values were at that location. Also,
 ! don't forget to deallocate any allocatable arrays as this will cause a
 ! failure when compiling the test with gnu compilers and address sanitizers.
 ! -------------------------------------------------------------------------

end program ftst_example
