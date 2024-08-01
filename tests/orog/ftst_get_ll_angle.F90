 program get_ll_angle

! Unit test for function get_lat_angle.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : get_lat_angle

 implicit none

 real            :: dlat, dy
 real, parameter :: EPSILON=0.001

! dy is the approximate distance in meters of one 
! degree of latitude.

 dy = 111139.0

 dlat = get_lat_angle(dy)

! Is dlat approximately one degree?

 if (abs(dlat - 1.0) > EPSILON) stop 2

 print*,"OK"

 print*,"SUCCESS"

 end program get_ll_angle
