 program check_get_xnsum

! Unit test for function get_xnsum, which computes the
! the number of high-resolution orography data points
! within a model grid box that are above the average
! terrain height.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : get_xnsum

 implicit none

 integer, parameter :: imn=360  ! i-dimension of high-res grid
 integer, parameter :: jmn=181  ! j-dimension of high-res grid

 integer            :: j
 integer            :: zavg(imn,jmn)  ! High-res grid terrain
 integer            :: zslm(imn,jmn)  ! High-res grid land mask

 real               :: delxn=360.0/float(imn) ! Resolution of high-res grid in degrees.
 real               :: glat(jmn)              ! Latitude of each high-res grid row in degrees.
 real               :: lon1,lat1,lon2,lat2
 real               :: xnsum

 print*,"Begin test of function get_xnsum"

! The high-res grid is a global one-degree lat/lon grid. Point (1,1)
! is the south pole/greenwich.

 do j = 1, jmn
   glat(j) = -90.0 + float(j-1) * delxn
 enddo

! Bounds of model grid box - straddles equator/greenwich.
! There are approximately 16 high-res points located within
! the model box.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

! First test point. The high-res grid is all ocean. Since all points in
! the model grid box will be sea level, there will be no points
! higher than the average of 0 meters.

 print*,"Test point 1."

 zslm = 0    ! high-res grid all water
 zavg = -999 ! high-res grid all sea level

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 0) stop 2

! Second test point. All high-res points are land with an elevation
! of 50 meters. There will be no points higher than the average.

 print*,"Test point 2."

 zslm = 1   ! high-res grid all land
 zavg = 50  ! constant elevation of 50 meters.

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 0) stop 4

! Third test point. All high-res points are land. One point is
! 100 meters, the rest are 50 meters. Therefore, there should
! be one point higher than the average.

 print*,"Test point 3."

 zslm = 1           ! all land
 zavg = 50          ! initialize elevation to 50 meters.
 zavg(359,91) = 100 ! one point within the model box is 100 meters.

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 1) stop 6

! Fourth test point. All high-res points are water. One point is
! 20 meters, which represents an inland lake. The other points
! are ocean. Therefore, there should be one point higher than the average.

 print*,"Test point 4."

 zslm = 0         ! all water
 zavg = -999      ! initialize to sea level.
 zavg(2,93) = 20  ! this represents an inland lake above sea level.

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 1) stop 8

! Fifth test point. Half of the high-res points within the grid
! box are land and half are water. Four high-res points are
! above the average.

 print*,"Test point 5."

 zslm = 0            ! Initialize to all water
 zavg = -999         ! Initialize to sea level.
 zslm(1:2,90:93) = 1 ! Half the points in the grid box are land.
 zavg(1,90:93) = 25  ! Set the elevation at the
 zavg(2,90) = 100    ! land points so that
 zavg(2,91) = 110    ! four points are above the average.
 zavg(2,92) = 107
 zavg(2,93) = 207

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 4) stop 10

 print*,"OK"

 print*,"SUCCESS"

 end program check_get_xnsum
