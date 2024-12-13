 program check_get_xnsum

 use orog_utils, only : get_xnsum

 implicit none

 integer, parameter :: imn=360
 integer, parameter :: jmn=181

 integer            :: j
 integer            :: zavg(imn,jmn)
 integer            :: zslm(imn,jmn)

 real               :: delxn=360.0/float(imn)
 real               :: glat(jmn)
 real               :: lon1,lat1,lon2,lat2
 real               :: xnsum

 print*,"Begin test."

! Set up a 'high-res' 1-degree grid.

 do j = 1, jmn
   glat(j) = -90.0 + float(j-1) * delxn
   print*,'j lat ',j,glat(j)
 enddo

! First test point. The high-res grid is all ocean. Since all points in
! the model grid box will be sea level, there will be no points
! higher than the average of 0 meters.

 print*,"Test point 1."

 zslm = 0    ! all water
 zavg = -999 ! all sea level

! Bounds of model grid box - straddles greenwich.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 0) stop 2

 print*,"Test point 2."

 zslm = 1    ! all land
 zavg = 50  ! constant elevation of 50 meters.

! Bounds of model grid box - straddles greenwich.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 0) stop 4

 print*,"Test point 3."

 zslm = 1    ! all land
 zavg = 50   ! constant elevation of 50 meters.
 zavg(359,91) = 100

! Bounds of model grid box - straddles greenwich.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 1) stop 6

 print*,"Test point 4."

 zslm = 0    ! all water
 zavg = -999   ! set to sea level.
 zavg(2,93) = 20 ! this represents an inland lake above sea level.

! Bounds of model grid box - straddles greenwich.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 1) stop 8

 print*,"Test point 5."

 zslm = 0      ! all water
 zavg = -999   ! set to sea level.
 zslm(1:2,90:93) = 1 ! half points in grid box land
 zavg(1,90:93) = 25 ! land points above sea level
 zavg(2,90) = 100 ! land points above sea level
 zavg(2,91) = 110 ! land points above sea level
 zavg(2,92) = 107 ! land points above sea level
 zavg(2,93) = 207 ! land points above sea level

! Bounds of model grid box - straddles greenwich.

 lon1 = -2.5
 lon2 = 2.5
 lat1 = -1.5
 lat2 = 1.5

 xnsum = get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat, zavg, zslm, delxn)

 if (nint(xnsum) /= 4) stop 10

 print*,"OK"

 print*,"SUCCESS"

 end program check_get_xnsum
