 program inside_polygon

! Unit test for function inside_a_polygon, which determines
! whether a test point is located within an area specified
! by a polygon.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only    : inside_a_polygon

 implicit none

 integer, parameter   :: npts=4

 real, parameter      :: D2R = 3.14159265358979/180.

 logical              :: inside

 real                 :: lon1, lat1
 real                 :: lon2(npts), lat2(npts)

 print*,"Starting test of inside_a_polygon."

! The first part of the function checks if the test point is well outside
! the neighborhood of the polygon - i.e., a gross check. There 
! are six separate checks. The first six tests are designed
! so that the polygon is far enough from the test point that
! each check is tripped.

! Test to trip the first 'if' gross check. Is point well outside the
! neighborhood in the plus 'x' direction?

 print*, "Test point 1"

 lon1 = 90.0 * D2R      ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 94.0 * D2R   ! Polygon.
 lat2(1) = -1.0 * D2R
 lon2(2) = 94.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 95.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 95.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 2  ! Test point should be outside polygon.

! Test to trip the second 'if' gross check. Is point well outside the
! neighborhood in the minus 'x' direction?

 print*, "Test point 2"

 lon1 = 90.0 * D2R     ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 84.0 * D2R  ! Polygon.
 lat2(1) = -1.0 * D2R
 lon2(2) = 84.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 85.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 85.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 4 ! Test point should be outside polygon.

! Test to trip the third 'if' gross check. Is point well outside the
! neighborhood in the plus 'y' direction?

 print*, "Test point 3"

 lon1 = 0.0 * D2R       ! Test point.
 lat1 = 0.0 * D2R       

 lon2(1) = 354.0 * D2R  ! Polygon.
 lat2(1) = -1.0 * D2R
 lon2(2) = 354.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 355.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 355.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 6 ! Test point should be outside polygon.

! Test to trip the fourth 'if' gross check. Is point well outside the
! neighborhood in the minus 'y' direction?

 print*, "Test point 4"

 lon1 = 0.0 * D2R       ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 4.0 * D2R    ! Polygon.
 lat2(1) = -1.0 * D2R
 lon2(2) = 4.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 5.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 5.0 * D2R
 lat2(4) = -1.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 8 ! Test point should be outside polygon.

! Test to trip the fifth 'if' gross check. Is point well outside the
! neighborhood in the plus 'z' direction?

 print*, "Test point 5"

 lon1 = 0.0 * D2R       ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 359.0 * D2R  ! Polygon.
 lat2(1) = -6.0 * D2R
 lon2(2) = 359.0 * D2R
 lat2(2) =  -5.0 * D2R
 lon2(3) = 1.0 * D2R
 lat2(3) =  -5.0 * D2R
 lon2(4) =  1.0 * D2R
 lat2(4) =  -6.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 10 ! Test point should be outside polygon.

! Test to trip the sixth 'if' gross check. Is point well outside the
! neighborhood in the minus 'z' direction?

 print*, "Test point 6"

 lon1 = 0.0 * D2R       ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 359.0 * D2R  ! Polygon.
 lat2(1) = 5.0 * D2R
 lon2(2) = 359.0 * D2R
 lat2(2) =  6.0 * D2R
 lon2(3) = 1.0 * D2R
 lat2(3) =  6.0 * D2R
 lon2(4) =  1.0 * D2R
 lat2(4) =  5.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 12 ! Test point should be outside polygon.

! Test the case when the test point coincides with a corner point
! of the polygon. Result should be 'inside'.

 print*, "Test point 7"

 lon1 = 90.0 * D2R      ! Test point.
 lat1 = 0.0 * D2R

 lon2(1) = 90.0 * D2R   ! Polygon.
 lat2(1) =  0.0 * D2R
 lon2(2) = 90.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 91.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 91.0 * D2R
 lat2(4) =  0.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (.not.inside) stop 14 ! Test point should be inside polygon.

! Test the case when the test point is inside the polygon.
! Here, the test point is exactly in the middle of a square
! polygon. That means the computed angle to each corner
! is 90 degrees.

 print*, "Test point 8"

 lon1 = 90.5 * D2R      ! Test point.
 lat1 = 0.5 * D2R

 lon2(1) = 90.0 * D2R   ! Polygon.
 lat2(1) =  0.0 * D2R
 lon2(2) = 90.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 91.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 91.0 * D2R
 lat2(4) =  0.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (.not.inside) stop 16  ! Test point should be inside polygon.

! Test the case when the test point just outside the polygon.

 print*, "Test point 9"

 lon1 = 89.5 * D2R      ! Test point.
 lat1 = 0.5 * D2R

 lon2(1) = 90.0 * D2R   ! Polygon.
 lat2(1) =  0.0 * D2R
 lon2(2) = 90.0 * D2R
 lat2(2) =  1.0 * D2R
 lon2(3) = 91.0 * D2R
 lat2(3) =  1.0 * D2R
 lon2(4) = 91.0 * D2R
 lat2(4) =  0.0 * D2R

 inside=inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 if (inside) stop 18 ! Test point should be outside polygon.

 print*,"OK"
 print*,"SUCCESS"

 end program inside_polygon
