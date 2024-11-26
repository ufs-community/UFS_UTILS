 program inside_polygon

 use orog_utils, only    : inside_a_polygon

 implicit none

 integer, parameter    :: npts=4

 real, parameter           :: D2R = 3.14159265358979/180.
 logical    :: inside

 real       :: lon1, lat1
 real       :: lon2(npts), lat2(npts)

! Test to trip the first 'if' range check

 print*, "Test point 1"

 lon1 = 90.0 * D2R
 lat1 = 0.0 * D2R

 lon2(1) = 94.0 * D2R
 lat2(1) = -1.0 * D2R
 lon2(2) = 94.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 95.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 95.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 2

! Test to trip the second 'if' range check

 print*, "Test point 2"

 lon1 = 90.0 * D2R
 lat1 = 0.0 * D2R

 lon2(1) = 84.0 * D2R
 lat2(1) = -1.0 * D2R
 lon2(2) = 84.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 85.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 85.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 4

 print*,"OK"
 print*,"SUCCSSS"

 end program inside_polygon
