 program enclosure
!
! Unit test for function enclosure_cnvx.
!
! Test the function using five test points
! and check against expected values.
!
! Authors N. Wang and G. Gayno
!
 implicit none

 real*8 :: v(2,4)
 real*8 :: p(2)
 real*8, parameter :: d2r = acos(-1.0)/180.0D0

 integer :: co_gc

 logical :: enclosure_cnvx, inside

 print*,"Start test of enclosure_cnvx"

! Lat/lon of vertices of the polygon.

 v(1,1) = 10.0D0*d2r; v(2,1) = 20.0D0*d2r
 v(1,2) = 15.0D0*d2r; v(2,2) = 30.0D0*d2r
 v(1,3) = 17.7D0*d2r; v(2,3) = 25.0D0*d2r
 v(1,4) = 20.0D0*d2r; v(2,4) = 20.0D0*d2r

! Points to test.

 print*,"Test point 1"
 p(1) = 17.7D0*d2r; p(2) = 25.000000001D0*d2r
 inside = enclosure_cnvx(v,4,p,co_gc)
 print*,'Point 1 ',inside,co_gc
 if (inside .or. co_gc /= 0) stop 2

 print*,"Test point 2"
 p(1) = 15.0D0*d2r; p(2) = 30.00000001D0*d2r
 inside = enclosure_cnvx(v,4,p,co_gc)
 print*,'Point 2 ',inside,co_gc
 if (inside .or. co_gc /= 0) stop 4

 print*,"Test point 3"
 p(1) = 20.00000000D0*d2r; p(2) = 20.0D0*d2r
 inside = enclosure_cnvx(v,4,p,co_gc)
 print*,'Point 3 ',inside,co_gc
 if (.not.inside .or. co_gc /= 4) stop 6

 print*,"Test point 4"
 p(1) = 9.999999999D0*d2r; p(2) = 20.0D0*d2r
 inside = enclosure_cnvx(v,4,p,co_gc)
 print*,'Point 4 ',inside,co_gc
 if (inside .or. co_gc /= 4) stop 8

 print*,"Test point 5"
 p(1) = 10.00000000*d2r; p(2) = 20.000000001D0*d2r
 inside = enclosure_cnvx(v,4,p,co_gc)
 print*,'Point 5 ',inside,co_gc
 if (inside .or. co_gc /= 0) stop 10

 print*,"OK"
 print*,"SUCCESS"

 end program enclosure
