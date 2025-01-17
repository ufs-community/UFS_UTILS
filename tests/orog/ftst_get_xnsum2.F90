 program check_get_xnsum2

! Unit test for routine get_xnsum2, which counts the
! number of high-resolution orography points within a
! model grid box, and counts the number of high-resolution
! points higher than a critical height. The critical
! height is a function of the standard deviation of height.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : get_xnsum2

 implicit none

 integer, parameter :: imn=360  ! i-dimension of high-res grid
 integer, parameter :: jmn=181  ! j-dimension of high-res grid
 
 integer            :: j
 integer            :: zavg(imn,jmn)  ! High-res grid terrain

 real               :: delxn=360.0/float(imn) ! Resolution of high-res grid in degrees.
 real               :: glat(jmn)              ! Latitude of each high-res grid row in degrees.
 real               :: lon1,lat1,lon2,lat2,hc
 real               :: xnsum1,xnsum2

! Variables holding the expected test results.

 real, parameter    :: EPSILON=0.001
 real               :: expected_hc
 integer            :: expected_xnsum1, expected_xnsum2

 data expected_hc /677.2/   ! Expected critical height in meters.
 data expected_xnsum1 /8/   ! Expected number of high-res pts in model
                            ! grid box that are above the critical height.
 data expected_xnsum2 /16/  ! Expected total number of high-res pts in model grid box.

 print*,"Begin test of routine get_xnsum2."

! The high-res grid is a global one-degree lat/lon grid. Point (1,1)
! is the south pole/greenwich.

 do j = 1, jmn
   glat(j) = -90.0 + float(j-1) * delxn
 enddo

! Bounds of model grid box - straddles equator/greenwich.
! There are 16 high-res points located within the model
! box. The i/j range of these 16 points is i=359,360,1,2
! and j=0,91,92,93.

 lon1 = -2.5  ! in degrees.
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

! Initialize high-res orography. Half of the points
! within the model grid box are at sea level and the
! other half at 1000 meters.

 zavg = -999  ! Flag value for sea level.
 zavg(359:360,90:93) = 1000

 call get_xnsum2(lon1,lat1,lon2,lat2,imn,jmn, &
                 glat,zavg,delxn,xnsum1,xnsum2,hc)

 if (nint(xnsum1) /= expected_xnsum1) stop 2
 if (nint(xnsum2) /= expected_xnsum2) stop 4
 if (abs(hc-expected_hc) > EPSILON) stop 6

 print*,"OK"

 print*,"SUCCESS"

 end program check_get_xnsum2
