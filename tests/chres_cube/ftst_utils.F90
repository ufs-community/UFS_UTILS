! Unit test for to_upper() and to_lower() functions under UFS_UTILS
! package, chres_cube utility.
!
! Lin Gan NCEP/EMC
 
program ftst_utils

 
 implicit none

 logical                          :: match_result

 character(len=12)                :: to_upper
 character(len=12)                :: test_input_char_1, test_input_char_2, u_st_base, l_st_base

 u_st_base="STAGGERLOCCE"
 l_st_base="staggerlocce"
 test_input_char_1="sTAGGErLOCCE"
 test_input_char_2="staGGErLOCCE"

 print*, "Starting Unit Testing to_upper_lower."
 print*, "testing to_lower and to_upper..."

!-------------------------------------------------------------------------
! Execute testing below by running target function with testing string 
!   When match_result set to be T - compare to base line is identical
!   When match_result set to be F - compare to base line is NOT identical
!-------------------------------------------------------------------------

 call to_lower(test_input_char_1) 
 match_result = test_input_char_1 == l_st_base
 if (.not.match_result) then
   stop 1
 endif

 if (to_upper(test_input_char_2) /= u_st_base) then
   stop 2
 endif

!-------------------------------------------------------------------------
! Display final result
!-------------------------------------------------------------------------

 print*, "OK"
 print*, "SUCCESS!"

end program ftst_utils
