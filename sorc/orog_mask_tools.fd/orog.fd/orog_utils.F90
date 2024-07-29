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

 contains

!> Convert from latitude and longitude to x,y,z coordinates.
!!
!! @param[in] siz Number of points to convert.
!! @param[in] lon Longitude (radians) of points to convert.
!! @param[in] lat Latitude (radians) of points to convert.
!! @param[out] x 'x' Coordinate of the converted points.
!! @param[out] y 'y' Coordinate of the converted points.
!! @param[out] z 'z' Coordinate of the converted points.
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

 end module orog_utils
