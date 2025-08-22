 program ftst_middle

 implicit none

 real*8 :: p1(2), p2(2), p(2)
 real*8, parameter :: d2r = acos(-1.0)/180.0D0

 print*,"Start test of middle"

 p1(1) = 0.0D0*d2r ; p1(2) = 90.0D0*d2r
 p2(1) = 1.0D0*d2r ; p2(2) = 90.0D0*d2r
 
 
 call middle (p1,p2,p)

 print*,'got here ',(p/d2r)

 print*,"OK"
 print*,"SUCCESS"

 end program ftst_middle
