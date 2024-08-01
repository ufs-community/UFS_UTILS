!> @file
!! @brief Utilities for orog code.
!! @author George Gayno NOAA/EMC

!> Module containing utilites used by the orog program.
!!
!! @author George Gayno NOAA/EMC
 module orog_utils

 implicit none

 private

 public :: latlon2xyz
 public :: minmax
 public :: get_lat_angle

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
 real, parameter    :: earth_radius = 6371200 ! meters
 real, parameter    :: rad2deg = 180./3.14159265358979

 get_lat_angle = dy/earth_radius*rad2deg

 end function get_lat_angle

 end module orog_utils
