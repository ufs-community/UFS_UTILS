!> @file
!! @brief Utilities for orog code.
!! @author George Gayno NOAA/EMC

!> Module containing utilites used by the orog program.
!!
!! @author George Gayno NOAA/EMC
 module orog_utils

 implicit none

 private

 real, parameter    :: earth_radius = 6371200.          !< earth radius in meters.
 real, parameter    :: pi=3.1415926535897931            !< pi.
 real, parameter    :: rad2deg = 180./3.14159265358979  !< radians per degrees.
 real, parameter    :: deg2rad = 3.14159265358979/180.  !< degrees per radians.

 public :: find_nearest_pole_points
 public :: find_poles
 public :: get_index
 public :: get_lat_angle
 public :: get_lon_angle
 public :: get_xnsum
 public :: get_xnsum2
 public :: get_xnsum3
 public :: inside_a_polygon
 public :: latlon2xyz
 public :: minmax
 public :: remove_isolated_pts
 public :: timef
 public :: transpose_orog
 public :: transpose_mask

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

!> Transpose the global orography data by flipping
!! the poles and moving the starting longitude to
!! Greenwich.
!!
!! @param[in] imn i-dimension of orography data.
!! @param[in] jmn j-dimension of orography data.
!! @param[inout] glob The global orography data.
!! @author G. Gayno

 subroutine transpose_orog(imn, jmn, glob)

 implicit none

 integer, intent(in)       :: imn, jmn
 integer(2), intent(inout) :: glob(imn,jmn)

 integer    :: i, j, it, jt
 integer(2) :: i2save

! Transpose from S to N to the NCEP standard N to S.

 do j=1,jmn/2
 do I=1,imn
   jt=jmn - j + 1
   i2save = glob(I,j)
   glob(I,j)=glob(I,jt)
   glob(I,jt) = i2save
 enddo
 enddo

! Data begins at dateline. NCEP standard is Greenwich.

 do j=1,jmn
 do I=1,imn/2
   it=imn/2 + i
   i2save = glob(i,J)
   glob(i,J)=glob(it,J)
   glob(it,J) = i2save
 enddo
 enddo

 end subroutine transpose_orog

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

!> Check if a point is inside a polygon.
!!
!! @param[in] lon1 Longitude of the point to check.
!! @param[in] lat1 Latitude of the point to check.
!! @param[in] npts Number of polygon vertices.
!! @param[in] lon2 Longitude of the polygon vertices.
!! @param[in] lat2 Latitude of the polygon vertices.
!! @return inside_a_polygon When true, point is within
!! the polygon.
!! @author GFDL programmer

 function inside_a_polygon(lon1, lat1, npts, lon2, lat2)

 implicit none

 logical inside_a_polygon

 real, parameter          :: EPSLN10 = 1.e-10
 real, parameter          :: EPSLN8 = 1.e-8
 real, parameter          :: RANGE_CHECK_CRITERIA=0.05

 integer, intent(in)      :: npts

 real, intent(in)         :: lon1, lat1
 real, intent(in)         :: lon2(npts), lat2(npts)

 integer                  :: i, ip1

 real                     :: anglesum, angle
 real                     :: x2(npts), y2(npts), z2(npts)
 real                     :: lon1_1d(1), lat1_1d(1)
 real                     :: x1(1), y1(1), z1(1)
 real                     :: pnt0(3),pnt1(3),pnt2(3)
 real                     :: max_x2,min_x2,max_y2,min_y2,max_z2,min_z2

! first convert to cartesian grid.

 call latlon2xyz(npts,lon2, lat2, x2, y2, z2);
 lon1_1d(1) = lon1
 lat1_1d(1) = lat1
 call latlon2xyz(1,lon1_1d, lat1_1d, x1, y1, z1);
 inside_a_polygon = .false.
 max_x2 = maxval(x2)
 if( x1(1) > max_x2+RANGE_CHECK_CRITERIA ) return
 min_x2 = minval(x2)
 if( x1(1)+RANGE_CHECK_CRITERIA < min_x2 ) return
 max_y2 = maxval(y2)
 if( y1(1) > max_y2+RANGE_CHECK_CRITERIA ) return
 min_y2 = minval(y2)
 if( y1(1)+RANGE_CHECK_CRITERIA < min_y2 ) return
 max_z2 = maxval(z2)
 if( z1(1) > max_z2+RANGE_CHECK_CRITERIA ) return
 min_z2 = minval(z2)
 if( z1(1)+RANGE_CHECK_CRITERIA < min_z2 ) return

 pnt0(1) = x1(1)
 pnt0(2) = y1(1)
 pnt0(3) = z1(1)

 anglesum = 0

 do i = 1, npts
   if(abs(x1(1)-x2(i)) < EPSLN10 .and.  &
      abs(y1(1)-y2(i)) < EPSLN10 .and.  &
      abs(z1(1)-z2(i)) < EPSLN10 ) then ! same as the corner point
     inside_a_polygon = .true.
     return
   endif
   ip1 = i+1
   if(ip1>npts) ip1 = 1
   pnt1(1) = x2(i)
   pnt1(2) = y2(i)
   pnt1(3) = z2(i)
   pnt2(1) = x2(ip1)
   pnt2(2) = y2(ip1)
   pnt2(3) = z2(ip1)
   angle = spherical_angle(pnt0, pnt2, pnt1);
   anglesum = anglesum + angle
 enddo

 if(abs(anglesum-2*PI) < EPSLN8) then
   inside_a_polygon = .true.
 else
   inside_a_polygon = .false.
 endif

 end function inside_a_polygon

!> Find the point on the model grid tile closest to the
!! north and south pole.
!!
!! @param[in] geolat Latitude on the supergrid.
!! @param[in] nx i-dimension of the supergrid.
!! @param[in] ny j-dimension of the supergrid.
!! @param[out] i_north_pole 'i' index of north pole. '0' if
!! pole is outside of grid.
!! @param[out] j_north_pole 'j' index of north pole. '0' if
!! pole is outside of grid.
!! @param[out] i_south_pole 'i' index of south pole. '0' if
!! pole is outside of grid.
!! @param[out] j_south_pole 'j' index of south pole. '0' if
!! pole is outside of grid.
!! @author GFDL Programmer
 subroutine find_poles(geolat, nx, ny, i_north_pole, j_north_pole, &
                       i_south_pole, j_south_pole)

 implicit none

 integer, intent(in)  :: nx, ny

 real, intent(in)     :: geolat(nx+1,ny+1)

 integer, intent(out) :: i_north_pole, j_north_pole
 integer, intent(out) :: i_south_pole, j_south_pole

 integer              :: i, j

 real                 :: maxlat, minlat

 print*,'- CHECK IF THE TILE CONTAINS A POLE.'

!--- figure out pole location.

 maxlat = -90
 minlat = 90
 i_north_pole = 0
 j_north_pole = 0
 i_south_pole = 0
 j_south_pole = 0
 do j = 1, ny+1; do i = 1, nx+1
   if( geolat(i,j) > maxlat ) then
      i_north_pole=i
      j_north_pole=j
      maxlat = geolat(i,j)
   endif
   if( geolat(i,j) < minlat ) then
      i_south_pole=i
      j_south_pole=j
      minlat = geolat(i,j)
   endif
 enddo ; enddo

!--- only when maxlat is close to 90. the point is north pole

 if(maxlat < 89.9 ) then
   i_north_pole = 0
   j_north_pole = 0
 endif
 if(minlat > -89.9 ) then
   i_south_pole = 0
   j_south_pole = 0
 endif

 print*, "- MINLAT=", minlat, "MAXLAT=", maxlat
 print*, "- NORTH POLE SUPERGRID INDEX IS ",  &
                i_north_pole, j_north_pole
 print*, "- SOUTH POLE SUPERGRID INDEX IS ",  &
                i_south_pole, j_south_pole

 end subroutine find_poles

!> Find the point on the model grid tile closest to the
!! north and south pole.
!!
!! @param[in] i_north_pole 'i' index of north pole. '0' if
!! pole is outside of grid.
!! @param[in] j_north_pole 'j' index of north pole. '0' if
!! pole is outside of grid.
!! @param[in] i_south_pole 'i' index of south pole. '0' if
!! pole is outside of grid.
!! @param[in] j_south_pole 'j' index of south pole. '0' if
!! pole is outside of grid.
!! @param[in] im i-dimension of model tile
!! @param[in] jm j-dimension of model tile
!! @param[out] is_north_pole 'true' for points surrounding the north pole.
!! @param[out] is_south_pole 'true' for points surrounding the south pole.
!! @author GFDL Programmer

 subroutine find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 implicit none

 integer, intent(in)     :: im, jm
 integer, intent(in)     :: i_north_pole, j_north_pole
 integer, intent(in)     :: i_south_pole, j_south_pole

 logical, intent(out)    :: is_north_pole(im,jm)
 logical, intent(out)    :: is_south_pole(im,jm)

 integer                 :: i, j

 print*,'- FIND NEAREST POLE POINTS.'

 is_north_pole=.false.
 is_south_pole=.false.

 if(i_south_pole >0 .and. j_south_pole > 0) then
   if(mod(i_south_pole,2)==0) then ! stretched grid
      do j = 1, JM; do i = 1, IM
        if(i==i_south_pole/2 .and. (j==j_south_pole/2  &
                    .or. j==j_south_pole/2+1) ) then
          is_south_pole(i,j) = .true.
          print*, "- SOUTH POLE AT I,J= ", i, j
        endif
      enddo; enddo
   else
      do j = 1, JM; do i = 1, IM
        if((i==i_south_pole/2 .or. i==i_south_pole/2+1) &
             .and. (j==j_south_pole/2 .or. &
             j==j_south_pole/2+1) ) then
          is_south_pole(i,j) = .true.
          print*, "- SOUTH POLE AT I,J= ", i, j
        endif
      enddo; enddo
    endif
 endif

 if(i_north_pole >0 .and. j_north_pole > 0) then
   if(mod(i_north_pole,2)==0) then ! stretched grid
      do j = 1, JM; do i = 1, IM
        if(i==i_north_pole/2 .and. (j==j_north_pole/2 .or.  &
           j==j_north_pole/2+1) ) then
          is_north_pole(i,j) = .true.
          print*, "- NORTH POLE AT I,J= ", i, j
        endif
      enddo; enddo
    else
      do j = 1, JM; do i = 1, IM
        if((i==i_north_pole/2 .or. i==i_north_pole/2+1) &
            .and. (j==j_north_pole/2 .or.  &
                   j==j_north_pole/2+1) ) then
          is_north_pole(i,j) = .true.
          print*, "- NORTH POLE AT I,J= ", i, j
        endif
      enddo; enddo
   endif
 endif

 end subroutine find_nearest_pole_points

!> Remove isolated model points.
!!
!! @param[in] im 'i' dimension of a model grid tile.
!! @param[in] jm 'j' dimension of a model grid tile.
!! @param[inout] slm Land-mask on the model tile.
!! @param[inout] oro Orography on the model tile.
!! @param[inout] var Standard deviation of orography on the model tile.
!! @param[inout] var4 Convexity on the model tile.
!! @param[inout] oa Orographic asymmetry on the model tile.
!! @param[inout] ol Orographic length scale on the model tile.
!! @author Jordan Alpert NOAA/EMC

 subroutine remove_isolated_pts(im,jm,slm,oro,var,var4,oa,ol)

 implicit none

 integer, intent(in)       :: im, jm

 real, intent(inout)       :: slm(im,jm)
 real, intent(inout)       :: oro(im,jm)
 real, intent(inout)       :: var(im,jm)
 real, intent(inout)       :: var4(im,jm)
 real, intent(inout)       :: oa(im,jm,4)
 real, intent(inout)       :: ol(im,jm,4)

 integer                   :: i, j, jn, js, k
 integer                   :: iw, ie, wgta, is, ise
 integer                   :: in, ine, inw, isw

 real                      :: slma, oroa, vara, var4a, xn, xs
 real, allocatable         :: oaa(:), ola(:)

!  REMOVE ISOLATED POINTS

 print*,"- REMOVE ISOLATED POINTS."

 allocate (oaa(4),ola(4))

 iso_loop : DO J=2,JM-1
   JN=J-1
   JS=J+1
   i_loop : DO I=1,IM
     IW=MOD(I+IM-2,IM)+1
     IE=MOD(I,IM)+1
     SLMA=SLM(IW,J)+SLM(IE,J)
     OROA=ORO(IW,J)+ORO(IE,J)
     VARA=VAR(IW,J)+VAR(IE,J)
     VAR4A=VAR4(IW,J)+VAR4(IE,J)
     DO K=1,4
       OAA(K)=OA(IW,J,K)+OA(IE,J,K)
! --- (*j*) fix typo:
       OLA(K)=OL(IW,J,K)+OL(IE,J,K)
     ENDDO
     WGTA=2
     XN=(I-1)+1
     IF(ABS(XN-NINT(XN)).LT.1.E-2) THEN
       IN=MOD(NINT(XN)-1,IM)+1
       INW=MOD(IN+IM-2,IM)+1
       INE=MOD(IN,IM)+1
       SLMA=SLMA+SLM(INW,JN)+SLM(IN,JN)+SLM(INE,JN)
       OROA=OROA+ORO(INW,JN)+ORO(IN,JN)+ORO(INE,JN)
       VARA=VARA+VAR(INW,JN)+VAR(IN,JN)+VAR(INE,JN)
       VAR4A=VAR4A+VAR4(INW,JN)+VAR4(IN,JN)+VAR4(INE,JN)
       DO K=1,4
         OAA(K)=OAA(K)+OA(INW,JN,K)+OA(IN,JN,K)+OA(INE,JN,K)
         OLA(K)=OLA(K)+OL(INW,JN,K)+OL(IN,JN,K)+OL(INE,JN,K)
       ENDDO
       WGTA=WGTA+3
     ELSE
       INW=INT(XN)
       INE=MOD(INW,IM)+1
       SLMA=SLMA+SLM(INW,JN)+SLM(INE,JN)
       OROA=OROA+ORO(INW,JN)+ORO(INE,JN)
       VARA=VARA+VAR(INW,JN)+VAR(INE,JN)
       VAR4A=VAR4A+VAR4(INW,JN)+VAR4(INE,JN)
       DO K=1,4
         OAA(K)=OAA(K)+OA(INW,JN,K)+OA(INE,JN,K)
         OLA(K)=OLA(K)+OL(INW,JN,K)+OL(INE,JN,K)
       ENDDO
       WGTA=WGTA+2
     ENDIF
     XS=(I-1)+1
     IF(ABS(XS-NINT(XS)).LT.1.E-2) THEN
       IS=MOD(NINT(XS)-1,IM)+1
       ISW=MOD(IS+IM-2,IM)+1
       ISE=MOD(IS,IM)+1
       SLMA=SLMA+SLM(ISW,JS)+SLM(IS,JS)+SLM(ISE,JS)
       OROA=OROA+ORO(ISW,JS)+ORO(IS,JS)+ORO(ISE,JS)
       VARA=VARA+VAR(ISW,JS)+VAR(IS,JS)+VAR(ISE,JS)
       VAR4A=VAR4A+VAR4(ISW,JS)+VAR4(IS,JS)+VAR4(ISE,JS)
       DO K=1,4
         OAA(K)=OAA(K)+OA(ISW,JS,K)+OA(IS,JS,K)+OA(ISE,JS,K)
         OLA(K)=OLA(K)+OL(ISW,JS,K)+OL(IS,JS,K)+OL(ISE,JS,K)
       ENDDO
       WGTA=WGTA+3
     ELSE
       ISW=INT(XS)
       ISE=MOD(ISW,IM)+1
       SLMA=SLMA+SLM(ISW,JS)+SLM(ISE,JS)
       OROA=OROA+ORO(ISW,JS)+ORO(ISE,JS)
       VARA=VARA+VAR(ISW,JS)+VAR(ISE,JS)
       VAR4A=VAR4A+VAR4(ISW,JS)+VAR4(ISE,JS)
       DO K=1,4
         OAA(K)=OAA(K)+OA(ISW,JS,K)+OA(ISE,JS,K)
         OLA(K)=OLA(K)+OL(ISW,JS,K)+OL(ISE,JS,K)
       ENDDO
       WGTA=WGTA+2
     ENDIF
     OROA=OROA/WGTA
     VARA=VARA/WGTA
     VAR4A=VAR4A/WGTA
     DO K=1,4
       OAA(K)=OAA(K)/WGTA
       OLA(K)=OLA(K)/WGTA
     ENDDO
     IF(SLM(I,J).EQ.0..AND.SLMA.EQ.WGTA) THEN
       PRINT '(" - SEA ",2F8.0," MODIFIED TO LAND",2F8.0,  &
               " AT ",2I8)',ORO(I,J),VAR(I,J),OROA,VARA,I,J
       SLM(I,J)=1.
       ORO(I,J)=OROA
       VAR(I,J)=VARA
       VAR4(I,J)=VAR4A
       DO K=1,4
         OA(I,J,K)=OAA(K)
         OL(I,J,K)=OLA(K)
       ENDDO
     ELSEIF(SLM(I,J).EQ.1..AND.SLMA.EQ.0.) THEN
       PRINT '(" - LAND",2F8.0," MODIFIED TO SEA ",2F8.0,  &
               " AT ",2I8)',ORO(I,J),VAR(I,J),OROA,VARA,I,J
       SLM(I,J)=0.
       ORO(I,J)=OROA
       VAR(I,J)=VARA
       VAR4(I,J)=VAR4A
       DO K=1,4
         OA(I,J,K)=OAA(K)
         OL(I,J,K)=OLA(K)
       ENDDO
     ENDIF
   ENDDO i_loop
 ENDDO iso_loop

 deallocate (oaa,ola)

 end subroutine remove_isolated_pts

!> Determine the location of a cubed-sphere point within
!! the high-resolution orography data.  The location is
!! described by the range of i/j indices on the high-res grid.
!!
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data set.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data set.
!! @param[in] npts Number of vertices to describe the cubed-sphere point.
!! @param[in] lonO The longitudes of the cubed-sphere vertices.
!! @param[in] latO The latitudes of the cubed-sphere vertices.
!! @param[in] delxn Resolution of the high-resolution orography
!! data set.
!! @param[out] jst Starting 'j' index on the high-resolution grid.
!! @param[out] jen Ending 'j' index on the high-resolution grid.
!! @param[out] ilist List of 'i' indices on the high-resolution grid.
!! @param[out] numx The number of 'i' indices on the high-resolution
!! grid.
!! @author GFDL programmer
 subroutine get_index(imn,jmn,npts,lonO,latO,delxn,  &
                jst,jen,ilist,numx)

 implicit none
 integer, intent(in)  :: IMN,JMN
 integer, intent(in)  :: npts
 real,    intent(in)  :: LONO(npts), LATO(npts)
 real,    intent(in)  :: DELXN
 integer, intent(out) :: jst,jen
 integer, intent(out) :: ilist(IMN)
 integer, intent(out) :: numx

 integer              :: i2, ii, ist, ien
 real                 :: minlat,maxlat,minlon,maxlon
        
 minlat = minval(LATO)
 maxlat = maxval(LATO)
 minlon = minval(LONO)
 maxlon = maxval(LONO)
 ist = minlon/DELXN+1
 ien = maxlon/DELXN+1
 jst = (minlat+90)/DELXN+1 
 jen = (maxlat+90)/DELXN 
!--- add a few points to both ends of j-direction
 jst = jst - 5
 if(jst<1) jst = 1
 jen = jen + 5
 if(jen>JMN) jen = JMN

!--- when around the pole, just search through all the points.
 if((jst == 1 .OR. jen == JMN) .and. &
    (ien-ist+1 > IMN/2) )then
   numx = IMN
   do i2 = 1, IMN
     ilist(i2) = i2
   enddo
 else if( ien-ist+1 > IMN/2 ) then  ! cross longitude = 0
!--- find the minimum that greater than IMN/2 
!--- and maximum that less than IMN/2
   ist = 0
   ien = IMN+1
   do i2 = 1, npts
     ii = LONO(i2)/DELXN+1
     if(ii <0 .or. ii>IMN) then
        print*,"- II=",ii,IMN,LONO(i2),DELXN
     endif
     if( ii < IMN/2 ) then
       ist = max(ist,ii)
     else if( ii > IMN/2 ) then
       ien = min(ien,ii)
     endif
   enddo 
   if(ist<1 .OR. ist>IMN) then
     print*, "FATAL ERROR: ist<1 .or. ist>IMN"
     call ABORT()
   endif
   if(ien<1 .OR. ien>IMN) then
     print*, "FATAL ERROR: iend<1 .or. iend>IMN"
     call ABORT()
   endif

   numx = IMN - ien + 1
   do i2 = 1, numx
     ilist(i2) = ien + (i2-1)           
   enddo
   do i2 = 1, ist
     ilist(numx+i2) = i2
   enddo
   numx = numx+ist
 else
   numx = ien-ist+1
   do i2 = 1, numx
     ilist(i2) = ist + (i2-1)
   enddo
 endif

 end subroutine get_index

!> Count the number of high-resolution orography points that
!! are higher than the model grid box average orography height.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] zslm  The high-resolution land mask.
!! @param[in] delxn Resolution of the high-res orography data.
!! @return get_xnsum The number of high-res points above the
!! mean orography.
!! @author GFDL Programmer

 function get_xnsum(lon1,lat1,lon2,lat2,imn,jmn, &
                    glat,zavg,zslm,delxn)

 implicit none

 real                    :: get_xnsum

 integer, intent(in)     :: imn,jmn
 integer, intent(in)     :: zavg(imn,jmn),zslm(imn,jmn)
 real, intent(in)        :: lon1,lat1,lon2,lat2,delxn
 real, intent(in)        :: glat(jmn)

 integer                 :: i, j, ist, ien, jst, jen, i1
 real                    :: oro, height
 real                    :: xland,xwatr,xl1,xs1,slm,xnsum

!-- Figure out ist,ien,jst,jen

 do j = 1, jmn
   if( glat(j) .gt. lat1 ) then
     jst = j
     exit
   endif
 enddo

 do j = 1, jmn
   if( glat(j) .gt. lat2 ) then
     jen = j
     exit
   endif
 enddo

 ist = lon1/delxn + 1
 ien = lon2/delxn
 if(ist .le.0) ist = ist + imn
 if(ien < ist) ien = ien + imn

!--- Compute average oro

 oro = 0.0
 xnsum = 0
 xland = 0
 xwatr = 0
 xl1 = 0
 xs1 = 0
 do j = jst,jen
   do i1 = 1, ien - ist + 1
     i = ist + i1 -1
     if( i .le. 0) i = i + imn
     if( i .gt. imn) i = i - imn
     xland = xland + float(zslm(i,j))
     xwatr = xwatr + float(1-zslm(i,j))
     xnsum = xnsum + 1.
     height = float(zavg(i,j))
     if(height.lt.-990.) height = 0.0
     xl1 = xl1 + height * float(zslm(i,j))
     xs1 = xs1 + height * float(1-zslm(i,j))
   enddo
 enddo

 if( xnsum > 1.) then
   slm = float(nint(xland/xnsum))
   if(slm.ne.0.) then
     oro= xl1 / xland
   else
     oro = xs1 / xwatr
   endif
 endif

 get_xnsum = 0
 do j = jst, jen
 do i1= 1, ien-ist+1
   i = ist + i1 -1
   if( i .le. 0) i = i + imn
   if( i .gt. IMN) i = i - imn
   height = float(zavg(i,j))
   if(height.lt.-990.) height = 0.0
   if ( height .gt. oro ) get_xnsum = get_xnsum + 1
 enddo
 enddo

 end function get_xnsum

!> Count the number of high-resolution orography points that
!! are higher than a critical value inside a model grid box
!! (or a portion of a model grid box). The critical value is a
!! function of the standard deviation of orography.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] delxn Resolution of the high-res orography data.
!! @param[out] xnsum1 The number of high-resolution orography
!! above the critical value inside a model grid box.
!! @param[out] xnsum2 The number of high-resolution orography
!! points inside a model grid box.
!! @param[out] hc Critical height.
!! @author GFDL Programmer

 subroutine get_xnsum2(lon1,lat1,lon2,lat2,imn,jmn, &
                   glat,zavg,delxn,xnsum1,xnsum2,hc)

 implicit none

 integer, intent(in) :: imn,jmn
 integer, intent(in) :: zavg(imn,jmn)

 real, intent(in)    :: lon1,lat1,lon2,lat2,delxn
 real, intent(in)    :: glat(jmn)
 real, intent(out)   :: xnsum1,xnsum2,hc

 integer             :: i, j, ist, ien, jst, jen, i1

 real                :: height, var
 real                :: xw1,xw2,xnsum

!-- Figure out ist,ien,jst,jen

 do j = 1, jmn
   if( glat(j) .gt. lat1 ) then
     jst = j
     exit
   endif
 enddo

 do j = 1, jmn
   if( glat(j) .gt. lat2 ) then
     jen = j
     exit
   endif
 enddo

 ist = lon1/delxn + 1
 ien = lon2/delxn
 if(ist .le.0) ist = ist + imn
 if(ien < ist) ien = ien + imn

!--- Compute average oro

 xnsum = 0
 xw1 = 0
 xw2 = 0
 do j = jst,jen
   do i1 = 1, ien - ist + 1
     i = ist + i1 -1
     if( i .le. 0) i = i + imn
     if( i .gt. imn) i = i - imn
     xnsum = xnsum + 1.
     height = float(zavg(i,j))
     if(height.lt.-990.) height = 0.0
     xw1 = xw1 + height
     xw2 = xw2 + height ** 2
   enddo
 enddo

 var = sqrt(max(xw2/xnsum-(xw1/xnsum)**2,0.))
 hc = 1116.2 - 0.878 * var
 xnsum1 = 0
 xnsum2 = 0
 do j = jst, jen
   do i1= 1, ien-ist+1
     i = ist + i1 -1
     if( i .le. 0) i = i + imn
     if( i .gt. imn) i = i - imn
     height = float(zavg(i,j))
     if ( height .gt. hc ) xnsum1 = xnsum1 + 1
     xnsum2 = xnsum2 + 1
   enddo
 enddo

 end subroutine get_xnsum2

!> Count the number of high-resolution orography points that
!! are higher than a critical value inside a model grid box
!! (or a portion of a model grid box). Unlike routine
!! get_xnsum2(), this routine does not compute the critical
!! value. Rather, it is passed in.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] delxn Resolution of the high-res orography data.
!! @param[out] xnsum1 The number of high-resolution orography
!! above the critical value inside a model grid box.
!! @param[out] xnsum2 The number of high-resolution orography
!! points inside a model grid box.
!! @param[in] hc Critical height.
!! @author GFDL Programmer

 subroutine get_xnsum3(lon1,lat1,lon2,lat2,imn,jmn, &
                      glat,zavg,delxn,xnsum1,xnsum2,HC)
 implicit none

 integer, intent(in) :: imn,jmn
 integer, intent(in) :: zavg(imn,jmn)

 real, intent(in)    :: glat(jmn)
 real, intent(in)    :: lon1,lat1,lon2,lat2,delxn
 real, intent(out)   :: xnsum1,xnsum2,hc

 integer             :: i, j, ist, ien, jst, jen, i1

 real                :: height

!-- Figure out ist,ien,jst,jen

! if lat1 or lat 2 is 90 degree. set jst = JMN

 jst = jmn
 jen = jmn
 do j = 1, jmn
   if( glat(j) .gt. lat1 ) then
     jst = j
     exit
   endif
 enddo

 do j = 1, jmn
   if( glat(j) .gt. lat2 ) then
     jen = j
     exit
   endif
 enddo

 ist = lon1/delxn + 1
 ien = lon2/delxn
 if(ist .le.0) ist = ist + imn
 if(ien < ist) ien = ien + imn

 xnsum1 = 0
 xnsum2 = 0
 do j = jst, jen
   do i1= 1, ien-ist+1
     i = ist + i1 -1
     if( i .le. 0) i = i + imn
     if( i .gt. imn) i = i - imn
     height = float(zavg(i,j))
     if ( height .gt. hc ) xnsum1 = xnsum1 + 1
     xnsum2 = xnsum2 + 1
   enddo
 enddo

 end subroutine get_xnsum3
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
