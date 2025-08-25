 program ftst_middle
! 
! Unit test for subroutine middle.
! Computes the 'middle' of a great circle
! given the two end points.
!
! Test two points against expected values.
!
! Author G. Gayno
!
 implicit none

 real*8 :: p1(2), p2(2), p(2)
 real*8, parameter :: d2r = acos(-1.0)/180.0D0
 real*8, parameter :: epsilon = 0.001D0

 print*,"Start test of middle"

 print*,"Test point 1"

 p1(1) = 0.0D0*d2r ; p1(2) = 359.0D0*d2r ! Cross greenwich at
 p2(1) = 0.0D0*d2r ; p2(2) = 1.0D0*d2r   ! the equator.
 
 call middle (p1,p2,p)

 p = p/d2r

 if (abs(p(1)-0.0D0) > epsilon) stop 2   ! Should be equator/greenwich
 if (abs(p(2)-0.0D0) > epsilon) stop 4

 print*,"Test point 2"

 p1(1) = 89.0D0*d2r ; p1(2) = 270.0D0*d2r ! Cross the north pole.
 p2(1) = 89.0D0*d2r ; p2(2) = 90.0D0*d2r
 
 call middle (p1,p2,p)

 p = p/d2r

 if (abs(p(1)-90.0D0) > epsilon) stop 6  ! Should be NP/180.0
 if (abs(p(2)-180.0D0) > epsilon) stop 8

 print*,"OK"
 print*,"SUCCESS"

 end program ftst_middle
