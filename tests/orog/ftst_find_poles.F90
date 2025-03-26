 program ftst_find_poles

! Unit test for routine find_poles.
!
! Checks if a model tile contains a 'pole'
! point based on the latitude of each point.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only    : find_poles

 implicit none

 integer, parameter     :: nx=9
 integer, parameter     :: ny=9

 integer                :: i_north_pole, j_north_pole
 integer                :: i_south_pole, j_south_pole

 real                   :: geolat(nx+1,ny+1)

 print*,"Starting test of find_poles."

! Point 1 - North Pole at (4,5)

 print*,"North pole test."

 geolat      = 89.0
 geolat(4,5) = 89.95  ! North pole

 call find_poles(geolat, nx, ny, i_north_pole, j_north_pole, &
                 i_south_pole, j_south_pole)

 if (i_north_pole /= 4) stop 2
 if (j_north_pole /= 5) stop 4
 if (i_south_pole /= 0) stop 6
 if (j_south_pole /= 0) stop 8

! Point 2 - South Pole at (2,8)

 print*,"South pole test."

 geolat      = -89.0
 geolat(2,8) = -89.95  ! South pole

 call find_poles(geolat, nx, ny, i_north_pole, j_north_pole, &
                 i_south_pole, j_south_pole)

 if (i_north_pole /= 0) stop 12
 if (j_north_pole /= 0) stop 14
 if (i_south_pole /= 2) stop 16
 if (j_south_pole /= 8) stop 18

! Point 3 - No pole points.

 print*,"No pole test."

 geolat      = -89.3

 call find_poles(geolat, nx, ny, i_north_pole, j_north_pole, &
                 i_south_pole, j_south_pole)

 if (i_north_pole /= 0) stop 22
 if (j_north_pole /= 0) stop 24
 if (i_south_pole /= 0) stop 26
 if (j_south_pole /= 0) stop 28

 print*,"OK"

 print*,"SUCCESS"

 end program ftst_find_poles
