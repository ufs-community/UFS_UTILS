!> @file
!! @brief Utilities for orog code.
!! @author George Gayno NOAA/EMC

!> Module containing utilites used by the orog program.
!!
!! @author George Gayno NOAA/EMC
 module orog_utils

 implicit none

 private

 real, parameter    :: earth_radius = 6371200. ! meters
 real, parameter    :: pi=3.1415926535897931
 real, parameter    :: rad2deg = 180./3.14159265358979
 real, parameter    :: deg2rad = 3.14159265358979/180.

 public :: latlon2xyz
 public :: minmax
 public :: get_lat_angle
 public :: get_lon_angle
 public :: timef
 public :: transpose_mask
 public :: spherical_angle

 contains

!> Print out the maximum and minimum values of
!! an array and optionally pass back the i/j
!! location of the maximum.
!!
!! @param[in] im The 'i' dimension of the array.
!! @param[in] jm The 'j' dimension of the array.
!! @param[in] a The array to check.
!! @param[in] title Name of the data to be checked.
!! @param[out] imax The 'i' location of the maximum.
!! @param[out] jmax The 'j' location of the maximum.
!!
!! @author Jordan Alpert NOAA/EMC
 subroutine minmax(im,jm,a,title,imax,jmax)

 implicit none

 character(len=8), intent(in)   :: title
 
 integer, intent(in)            :: im, jm
 integer, intent(out), optional :: imax, jmax

 real,    intent(in)            :: a(im,jm)

 integer                        :: i, j

 real                           :: rmin,rmax

 rmin=huge(a)
 rmax=-rmin

 if (present(imax) .and. present(jmax)) then
   imax=0
   jmax=0
 endif

 do j=1,jm
   do i=1,im
     if(a(i,j) >= rmax) then
       rmax=a(i,j)
       if (present(imax) .and. present(jmax)) then
         imax = i
         jmax = j
       endif
     endif
     if(a(i,j) <= rmin)rmin=a(i,j)
   enddo
 enddo

 write(6,150) title,rmin,rmax
150 format(' - ',a8,' MIN=',e13.4,2x,'MAX=',e13.4)

 end subroutine minmax

!> Convert from latitude and longitude to x,y,z coordinates.
!!
!! @param[in] siz Number of points to convert.
!! @param[in] lon Longitude (radians) of points to convert.
!! @param[in] lat Latitude (radians) of points to convert.
!! @param[out] x 'x' Coordinate of the converted points.
!! @param[out] y 'y' Coordinate of the converted points.
!! @param[out] z 'z' Coordinate of the converted points.
!!
!! @author GFDL programmer
 subroutine latlon2xyz(siz,lon, lat, x, y, z)

 implicit none

 integer, intent(in) :: siz
 real, intent(in)    :: lon(siz), lat(siz)
 real, intent(out)   :: x(siz), y(siz), z(siz)

 integer             :: n

 do n = 1, siz
   x(n) = cos(lat(n))*cos(lon(n))
   y(n) = cos(lat(n))*sin(lon(n))
   z(n) = sin(lat(n))
 enddo

 end subroutine latlon2xyz

!> Convert the 'y' direction distance of a cubed-sphere grid
!! point to the corresponding distance in latitude.
!!
!! @param[in] dy Distance along the 'y' direction of a cubed-sphere
!! point in meters.
!! @return get_lat_angle Corresponding latitudinal distance in degrees.
!!
!! @author GFDL programmer

 function get_lat_angle(dy)

 implicit none

 real, intent(in)   :: dy

 real               :: get_lat_angle

 get_lat_angle = dy/earth_radius*rad2deg

 end function get_lat_angle

!> Convert the 'x' direction distance of a cubed-sphere grid
!! point to the corresponding distance in longitude.
!!
!! @param[in] dx Distance along the 'x' direction of a
!! cubed-sphere grid point in meters.
!! @param[in] lat_in Latitude of the cubed-sphere point in
!! degrees.
!! @return get_lon_angle Corresponding distance in longitude
!! in degrees.
!!
!! @author GFDL programmer

 function get_lon_angle(dx,lat_in)

 implicit none

 real, intent(in)     :: dx, lat_in

 real                 :: get_lon_angle, lat

 lat = lat_in
 if (lat > 89.5) lat = 89.5
 if (lat < -89.5) lat = -89.5

 get_lon_angle = 2*asin( sin(dx/earth_radius*0.5)/cos(lat*deg2rad) )*rad2deg

 end function get_lon_angle

!> Transpose the global landmask by flipping
!! the poles and moving the starting longitude to
!! Greenwich.
!!
!! @param[in] imn i-dimension of landmask data.
!! @param[in] jmn j-dimension of landmask data.
!! @param[inout] mask The global landmask data.
!! @author G. Gayno

 subroutine transpose_mask(imn, jmn, mask)

 implicit none

 integer, intent(in)       :: imn, jmn
 integer(1), intent(inout) :: mask(imn,jmn)

 integer    :: i, j, it, jt
 integer(1) :: isave

! Transpose from S to N to the NCEP standard N to S.

 do j=1,jmn/2
 do I=1,imn
   jt=jmn - j + 1
   isave = mask(I,j)
   mask(I,j)=mask(I,jt)
   mask(I,jt) = isave
 enddo
 enddo

! Data begins at dateline. NCEP standard is Greenwich.

 do j=1,jmn
 do I=1,imn/2
   it=imn/2 + i
   isave = mask(i,J)
   mask(i,J)=mask(it,J)
   mask(it,J) = isave
 enddo
 enddo

 end subroutine transpose_mask

!> Compute spherical angle.
!!
!! @param[in] v1 Vector 1.
!! @param[in] v2 Vector 2.
!! @param[in] v3 Vector 3.
!! @return spherical_angle Spherical Angle.
!! @author GFDL programmer

 function spherical_angle(v1, v2, v3)

 implicit none

 real              :: spherical_angle

 real, parameter   :: EPSLN30 = 1.e-30

 real, intent(in)  :: v1(3), v2(3), v3(3)
 
 real              :: px, py, pz, qx, qy, qz, ddd

! vector product between v1 and v2

 px = v1(2)*v2(3) - v1(3)*v2(2)
 py = v1(3)*v2(1) - v1(1)*v2(3)
 pz = v1(1)*v2(2) - v1(2)*v2(1)

! vector product between v1 and v3

 qx = v1(2)*v3(3) - v1(3)*v3(2);
 qy = v1(3)*v3(1) - v1(1)*v3(3);
 qz = v1(1)*v3(2) - v1(2)*v3(1);

 ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
 if ( ddd <= 0.0 ) then
   spherical_angle = 0.
 else
   ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
   if( abs(ddd-1) < EPSLN30 ) ddd = 1;
   if( abs(ddd+1) < EPSLN30 ) ddd = -1;
   if ( ddd>1. .or. ddd<-1. ) then
    !FIX to correctly handle co-linear points (angle near pi or 0) */
     if (ddd < 0.) then
       spherical_angle = PI
     else
       spherical_angle = 0.
     endif
   else
     spherical_angle = acos( ddd )
   endif
 endif

 end function spherical_angle

!> Get the date/time from the system clock.
!!
!! @return timef
!! @author Mark Iredell

 real function timef()

 implicit none

 character(8)          :: date
 character(10)         :: time
 character(5)          :: zone
 integer,dimension(8)  :: values
 integer               :: total
 real                  :: elapsed

 call date_and_time(date,time,zone,values)
 total=(3600*values(5)) + (60*values(6))+values(7)
 elapsed=float(total) + (1.0e-3*float(values(8)))
 timef=elapsed

 end function timef

 end module orog_utils
