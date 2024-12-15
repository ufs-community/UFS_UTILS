 program ll2xyz

! Unit test for routine latlon2xyz, which converts
! lat/lon to x/y/z coordinates.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : latlon2xyz

 implicit none

 integer, parameter :: siz = 6

 real, parameter    :: d2r = 3.14159265358979/180.
 real, parameter    :: EPSILON=0.0001

 integer    :: j

 real       :: lon(siz), lat(siz), x(siz), y(siz), z(siz)
 real       :: expected_x_component(siz)
 real       :: expected_y_component(siz)
 real       :: expected_z_component(siz)

! These are the expected x/y/z components returned from
! latlon2xyz for our test points.

 data expected_x_component/1.0, 0.0, -1.0, &
                           0.0, 0.0, 0.7071068/

 data expected_y_component/0.0, 1.0,  0.0, &
                           -1.0, 0.0, 0.0/

 data expected_z_component/0.0, 0.0, 0.0, &
                           0.0, 1.0, -0.7071068/

 print*,"Starting test of latlon2xyz."

! Test point 1 - the equator/greenwich.

 lat(1) = 0.0
 lon(1) = 0.0

! Test point 2 - the equator/90E

 lat(2) = 0.0
 lon(2) = 90.0

! Test point 3 - the equator/dateline

 lat(3) = 0.0
 lon(3) = 180.0

! Test point 4 - the equator/90W

 lat(4) = 0.0
 lon(4) = 270.0

! Test point 5 - the north pole/greenwich

 lat(5) = 90.0
 lon(5) = 0.0

! Test point 6 - 45S/greenwich

 lat(6) = -45.0
 lon(6) = 0.0

 lat = lat * d2r
 lon = lon * d2r

! Call the routine to unit test.

 call latlon2xyz(siz,lon,lat,x,y,z)

! Check results.

 do j = 1, siz
   if (abs(x(j) - expected_x_component(j)) > EPSILON) stop 2
   if (abs(y(j) - expected_y_component(j)) > EPSILON) stop 3
   if (abs(z(j) - expected_z_component(j)) > EPSILON) stop 4
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program ll2xyz
