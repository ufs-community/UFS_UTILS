 program get_ll_angle

! Unit test for functions get_lat_angle and
! get_lon_angle.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : get_lat_angle, get_lon_angle

 implicit none

 real            :: dlat, dlon, dy, lat
 real, parameter :: EPSILON=0.001

 print*,'Test get_lat_angle'

! dy is the approximate distance in meters of one 
! degree of latitude (or longitude at the equator).

 dy = 111139.0

 dlat = get_lat_angle(dy)

! Is dlat approximately one degree?

 if (abs(dlat - 1.0) > EPSILON) stop 2

 print*,'Test get_lon_angle'

! Test equator point. Should be about 1-degree.

 lat = 0.0
 dlon = get_lon_angle(dy,lat)
 if (abs(dlon - 1.0) > EPSILON) stop 3

! Test point at 60S. Should be about 2-degrees.

 lat = -60.0
 dlon = get_lon_angle(dy,lat)
 if (abs(dlon - 2.0) > EPSILON) stop 4

! Test both poles. To prevent a divide by zero,
! the function has special logic at the poles.
! The result is about 176 degrees.
 
 lat = -90.0
 dlon = get_lon_angle(dy,lat)
 if (abs(dlon - 176.254) > EPSILON) stop 5

 lat =  90.0
 dlon = get_lon_angle(dy,lat)
 if (abs(dlon - 176.254) > EPSILON) stop 6

 print*,"OK"

 print*,"SUCCESS"

 end program get_ll_angle
