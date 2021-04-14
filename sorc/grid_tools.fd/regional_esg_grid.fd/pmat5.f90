!> @file
!! @brief Handy geographical transformations.
!! @author R. J. Purser @date 1996
 
!> Constants for orientation and stretching of map.
!! @author R. J. Purser
module cstgeo ! Constants for orientation and stretching of map
use pkind, only: sp
implicit none
real(sp),dimension(3,3):: rotm !< Orthogonal rotation matrix
real(sp)               :: sc !< Schmidt stretching factor
real(sp)               :: sci !< Schmidt inverse stretching factor
end module cstgeo

!> Constants for orientation and stretching of map.
!! @author R. J. Purser
module dcstgeo ! Constants for orientation and stretching of map
use pkind, only: dp
implicit none
real(dp),dimension(3,3):: rotm !< Orthogonal rotation matrix
real(dp)               :: sc !< Schmidt stretching factor
real(dp)               :: sci !< Schmidt inverse stretching factor
end module dcstgeo

!> Utility routines for orienting the globe and basic geographical mappings.
!! @author R. J. Purser
module pmat5
use pkind, only: spi,sp,dp
implicit none
private
public :: ininmap,inivmap,ctogr,grtoc,ctog,gtoc,&
          gtoframe,paraframe,frametwist,&
          ctoc_schm,plrot,plroti,plctoc
interface ininmap;   module procedure sininmap,dininmap;          end interface
interface inivmap;   module procedure sinivmap,dinivmap;          end interface
interface ctogr;     module procedure sctogr,   dctogr;           end interface
interface grtoc
   module procedure sgrtoc,dgrtoc, sgrtocd,dgrtocd, sgrtocdd,dgrtocdd
                                                                  end interface
interface ctog;      module procedure sctog,   dctog;             end interface
interface gtoc
   module procedure sgtoc,dgtoc, sgtocd,dgtocd, sgtocdd,dgtocdd;  end interface
interface gtoframe
   module procedure sgtoframev,gtoframev,sgtoframem,gtoframem;    end interface
interface paraframe; module procedure sparaframe,paraframe;       end interface
interface frametwist;module procedure sframetwist,frametwist;     end interface
interface ctoc_schm 
   module procedure sctoc,dctoc, sctocd,dctocd, sctocdd,dctocdd;  end interface
interface plrot;     module procedure plrot, dplrot;              end interface
interface plroti;    module procedure plroti,dplroti;             end interface
interface plctoc;    module procedure plctoc;                     end interface
contains

!> Initialize the rotation matrix ROT3 needed to transform standard
!! earth-centered cartesian components to the alternative cartesian frame
!! oriented so as to put geographical point (ALAT0,ALON0) on the projection
!! axis. Single precision version.
!!
!! @param[in] alon0 geographical longitude (degrees) of projection center.
!! @param[in] alat0 geographical latitude (degrees) of projection center.
!! @param[out] rot3 rotation matrix.
!! @author R. J. Purser @date 1995
subroutine sininmap(alon0,alat0,rot3)!                               [ininmap]
use pietc_s, only: u0,dtor
implicit none
real(sp),               intent(IN ):: alon0,alat0
real(sp),dimension(3,3),intent(OUT):: rot3
real(sp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
blon0=dtor*alon0; clon0=cos(blon0); slon0=sin(blon0)
blat0=dtor*alat0; clat0=cos(blat0); slat0=sin(blat0)
rot3(1,1)=slat0*clon0; rot3(1,2)=slat0*slon0; rot3(1,3)=-clat0
rot3(2,1)=-slon0;      rot3(2,2)=clon0;       rot3(2,3)=u0
rot3(3,1)=clat0*clon0; rot3(3,2)=clat0*slon0; rot3(3,3)=slat0
end subroutine sininmap

!> Initialize the rotation matrix ROT3 needed to transform standard
!! earth-centered cartesian components to the alternative cartesian frame
!! oriented so as to put geographical point (ALAT0,ALON0) on the projection
!! axis. Double precision version.
!!
!! @param[in] alon0 geographical longitude (degrees) of projection center.
!! @param[in] alat0 geographical latitude (degrees) of projection center.
!! @param[out] rot3 rotation matrix.
!! @author R. J. Purser
subroutine dininmap(alon0,alat0,rot3)!                               [ininmap]
use pietc, only: u0,dtor
implicit none
real(dp),               intent(IN ):: alon0,alat0
real(dp),dimension(3,3),intent(OUT):: rot3
real(dp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
blon0=dtor*alon0; clon0=cos(blon0); slon0=sin(blon0)
blat0=dtor*alat0; clat0=cos(blat0); slat0=sin(blat0)
rot3(1,1)=slat0*clon0; rot3(1,2)=slat0*slon0; rot3(1,3)=-clat0
rot3(2,1)=-slon0;      rot3(2,2)=clon0;       rot3(2,3)=u0
rot3(3,1)=clat0*clon0; rot3(3,2)=clat0*slon0; rot3(3,3)=slat0
end subroutine dininmap

!> Initialize the rotation matrix ROT3 needed to transform standard
!! earth-centered cartesian components to the alternative cartesian frame
!! oriented so as to put geographical point (ALAT0,ALON0) at the viewing
!! nadir. Single precision version.
!!
!! @param[in] alon0 geographical longitude (degrees) of viewing nadir.
!! @param[in] alat0 geographical latitude (degrees) of viewing nadir.
!! @param[out] rot3 rotation matrix.
!! @author R. J. Purser @date 1995
subroutine sinivmap(alon0,alat0,rot3)!                               [inivmap]
use pietc_s, only: u0,dtor
implicit none
real(sp),               intent(IN ):: alon0,alat0
real(sp),dimension(3,3),intent(OUT):: rot3
real(sp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
blon0=dtor*alon0
blat0=dtor*alat0
clon0=cos(blon0)
slon0=sin(blon0)
clat0=cos(blat0)
slat0=sin(blat0)
rot3(1,1)=-slon0
rot3(1,2)=clon0
rot3(1,3)=u0
rot3(2,1)=-slat0*clon0
rot3(2,2)=-slat0*slon0
rot3(2,3)=clat0
rot3(3,1)=clat0*clon0
rot3(3,2)=clat0*slon0
rot3(3,3)=slat0
end subroutine sinivmap

!> Initialize the rotation matrix ROT3 needed to transform standard
!! earth-centered cartesian components to the alternative cartesian frame
!! oriented so as to put geographical point (ALAT0,ALON0) at the viewing
!! nadir. Double precision version.
!!
!! @param[in] alon0 geographical longitude (degrees) of viewing nadir.
!! @param[in] alat0 geographical latitude (degrees) of viewing nadir.
!! @param[out] rot3 rotation matrix.
!! @author R. J. Purser
subroutine dinivmap(alon0,alat0,rot3)!                               [inivmap]
use pietc, only: u0,dtor
implicit none
real(dp),               intent(IN ):: alon0,alat0
real(dp),dimension(3,3),intent(OUT):: rot3
real(dp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
blon0=dtor*alon0
blat0=dtor*alat0
clon0=cos(blon0)
slon0=sin(blon0)
clat0=cos(blat0)
slat0=sin(blat0)
rot3(1,1)=-slon0
rot3(1,2)=clon0
rot3(1,3)=u0
rot3(2,1)=-slat0*clon0
rot3(2,2)=-slat0*slon0
rot3(2,3)=clat0
rot3(3,1)=clat0*clon0
rot3(3,2)=clat0*slon0
rot3(3,3)=slat0
end subroutine dinivmap

!> Transform "Cartesian" to "Geographical" coordinates, where the
!! geographical coordinates refer to latitude and longitude (radians)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Single precision version.
!!
!! @param[in] xe Earth-centered cartesian unit 3-vector.
!! @param[out] rlat radians latitude.
!! @param[out] rlon radians longitude.
!! @author R. J. Purser
subroutine sctogr(xe,rlat,rlon)!                                       [ctogr]
use pietc_s, only: u0
implicit none
real(sp),dimension(3),intent(IN ):: xe
real(sp),             intent(OUT):: rlat,rlon
real(sp)                         :: r
r=sqrt(xe(1)**2+xe(2)**2)
rlat=atan2(xe(3),r)
if(r==u0)then
   rlon=u0
else
   rlon=atan2(xe(2),xe(1))
endif
end subroutine sctogr

!> Transform "Cartesian" to "Geographical" coordinates, where the
!! geographical coordinates refer to latitude and longitude (radians)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Double precision version.
!!
!! @param[in] xe Earth-centered cartesian unit 3-vector.
!! @param[out] rlat radians latitude.
!! @param[out] rlon radians longitude.
!! @author R. J. Purser
subroutine dctogr(xe,rlat,rlon)!                                       [ctogr]
use pietc, only: u0
implicit none
real(dp),dimension(3),intent(IN ):: xe
real(dp),             intent(OUT):: rlat,rlon
real(dp)                         :: r
r=sqrt(xe(1)**2+xe(2)**2)
rlat=atan2(xe(3),r)
if(r==u0)then
   rlon=u0
else
   rlon=atan2(xe(2),xe(1))
endif
end subroutine dctogr

!> Transform "Geographical" to "Cartesian" coordinates.
!! Single proecision version.
!!
!! @param[in] rlat Latitude (radians) of point.
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe Earth-centered cartesian unit 3-vector of point.
!! @author R. J. Purser
subroutine sgrtoc(rlat,rlon,xe)!                                       [grtoc]
implicit none
real(sp),             intent(IN ):: rlat,rlon
real(sp),dimension(3),intent(OUT):: xe
real(sp)                         :: sla,cla,slo,clo
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine sgrtoc

!> Transform "Geographical" to "Cartesian" coordinates.
!! Double proecision version.
!!
!! @param[in] rlat Latitude (radians) of point. 
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe Earth-centered cartesian unit 3-vector of point.
!! @author R. J. Purser
subroutine dgrtoc(rlat,rlon,xe)!                                       [grtoc]
implicit none
real(dp),             intent(IN ):: rlat,rlon
real(dp),dimension(3),intent(OUT):: xe
real(dp)                         :: sla,cla,slo,clo
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine dgrtoc

!> Transform "Geographical" to "Cartesian" coordinates, together with the
!! partial derivatives of cartesians wrt latitude and longitude.
!! Single precision version.
!!
!! @param[in] rlat Latitude (radians) of point.
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe Earth-centered cartesian unit 3-vector of point.
!! @param[out] dxedlat Derivative, d(xe)/d(rlat).
!! @param[out] dxedlon Derivative, d(xe)/d(rlon).
!! @author R. J. Purser
subroutine sgrtocd(rlat,rlon,xe,dxedlat,dxedlon)!                      [grtoc]
implicit none
real(sp),             intent(IN ):: rlat,rlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
real(dp)             :: rlat_d,rlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d
rlat_d=rlat; rlon_d=rlon
call dgrtocd(rlat_d,rlon_d,xe_d,dxedlat_d,dxedlon_d)
xe     =xe_d
dxedlat=dxedlat_d
dxedlon=dxedlon_d
end subroutine sgrtocd

!> Transform "Geographical" to "Cartesian" coordinates, together with the
!! partial derivatives of cartesians wrt latitude and longitude.
!! Double precision version.
!!
!! @param[in] rlat Latitude (radians) of point.
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe Earth-centered cartesian unit 3-vector of point.
!! @param[out] dxedlat Derivative, d(xe)/d(rlat).
!! @param[out] dxedlon Derivative, d(xe)/d(rlon).
!! @author R. J. Purser
subroutine dgrtocd(rlat,rlon,xe,dxedlat,dxedlon)!                      [grtoc]
use pietc, only: u0
implicit none
real(dp),             intent(IN ):: rlat,rlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
real(dp)                         :: sla,cla,slo,clo
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=u0 
end subroutine dgrtocd

!> Transform "Geographical" to "Cartesian" coordinates, together with the
!! 1st and 2nd partial derivatives of cartesians wrt latitude and longitude.
!! Single precision version.
!!
!! @param[in] rlat Latitude (radians) of point.
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe Earth-centered cartesian unit 3-vector of point.
!! @param[out] dxedlat Derivative, d(xe)/d(rlat).
!! @param[out] dxedlon Derivative, d(xe)/d(rlon).
!! @param[out] ddxedlatdlat Derivative, d^2(xe)/(d(rlat)d(rlat)).
!! @param[out] ddxedlatdlon Derivative, d^2(xe)/(d(rlat)d(rlon)).
!! @param[out] ddxedlondlon Derivative, d^2(xe)/(d(rlon)d(rlon)).
!! @author R. J. Purser
subroutine sgrtocdd(rlat,rlon,xe,dxedlat,dxedlon, &!                   [grtoc]
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)
implicit none
real(sp),             intent(IN ):: rlat,rlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
real(dp)             :: rlat_d,rlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d, &
                        ddxedlatdlat_d,ddxedlatdlon_d,ddxedlondlon_d
rlat_d=rlat; rlon_d=rlon
call dgtocdd(rlat_d,rlon_d,xe_d,dxedlat_d,dxedlon_d, &
     ddxedlatdlat_d,ddxedlatdlon_d,ddxedlondlon_d)
xe          =xe_d
dxedlat     =dxedlat_d
dxedlon     =dxedlon_d
ddxedlatdlat=ddxedlatdlat_d
ddxedlatdlon=ddxedlatdlon_d
ddxedlondlon=ddxedlondlon_d
end subroutine sgrtocdd

!> Transform "Geographical" to "Cartesian" coordinates, together with the
!! 1st and 2nd partial derivatives of cartesians wrt latitude and longitude.
!! Double precision version.
!!
!! @param[in] rlat Latitude (radians) of point.
!! @param[in] rlon Longitude (radians) of point.
!! @param[out] xe  Earth-centered cartesian unit 3-vector of point.
!! @param[out] dxedlat Derivative, d(xe)/d(rlat).
!! @param[out] dxedlon Derivative, d(xe)/d(rlon).
!! @param[out] ddxedlatdlat Derivative, d^2(xe)/(d(rlat)d(rlat)).
!! @param[out] ddxedlatdlon Derivative, d^2(xe)/(d(rlat)d(rlon)).
!! @param[out] ddxedlondlon Derivative, d^2(xe)/(d(rlon)d(rlon)).
!! @author R. J. Purser
subroutine dgrtocdd(rlat,rlon,xe,dxedlat,dxedlon, &!                   [grtoc]
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)
use pietc, only: u0
implicit none
real(dp),             intent(IN ):: rlat,rlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
real(dp)                         :: sla,cla,slo,clo
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=u0 
ddxedlatdlat(1)=-cla*clo
ddxedlatdlat(2)=-cla*slo
ddxedlatdlat(3)=-sla
ddxedlatdlon(1)= sla*slo
ddxedlatdlon(2)=-sla*clo
ddxedlatdlon(3)= u0
ddxedlondlon(1)=-cla*clo
ddxedlondlon(2)=-cla*slo
ddxedlondlon(3)= u0
end subroutine dgrtocdd

!> Transform "Cartesian" to "Geographical" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Single precision version.
!!
!! @param[in] xe Earth-centered cartesian unit 3-vector.
!! @param[out] dlat degrees latitude.
!! @param[out] dlon degrees longitude.
!! @author R. J. Purser @date 1994
subroutine sctog(xe,dlat,dlon)!                                         [ctog]
use pietc_s, only: u0,rtod
implicit none
real(sp),dimension(3),intent(IN ):: xe
real(sp),             intent(OUT):: dlat,dlon
real(sp)                         :: r
r=sqrt(xe(1)**2+xe(2)**2)
dlat=atan2(xe(3),r)*rtod
if(r==u0)then
   dlon=u0
else
   dlon=atan2(xe(2),xe(1))*rtod
endif
end subroutine sctog

!> Transform "Cartesian" to "Geographical" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Double precision version.
!!
!! @param[in] xe Earth-centered cartesian unit 3-vector.
!! @param[out] dlat degrees latitude.
!! @param[out] dlon degrees longitude.
!! @author R. J. Purser
subroutine dctog(xe,dlat,dlon)!                                         [ctog]
use pietc, only: u0,rtod
implicit none
real(dp),dimension(3),intent(IN ):: xe
real(dp),             intent(OUT):: dlat,dlon
real(dp)                         :: r
r=sqrt(xe(1)**2+xe(2)**2)
dlat=atan2(xe(3),r)*rtod
if(r==u0)then
   dlon=u0
else
   dlon=atan2(xe(2),xe(1))*rtod
endif
end subroutine dctog

!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Single precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe Earth-centered cartesian unit 3-vector.
!! @author R. J. Purser @date 1994
subroutine sgtoc(dlat,dlon,xe)!                                         [gtoc]
use pietc_s, only: dtor
implicit none
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe
real(sp)                         :: rlat,rlon,sla,cla,slo,clo
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine sgtoc

!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Double precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe cartesian unit 3-vector.
!! @author R. J. Purser
subroutine dgtoc(dlat,dlon,xe)!                                         [gtoc]
use pietc, only: dtor
implicit none
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine dgtoc
 
!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Also, return the partial derivatives of xe wrt latitude and longitude.
!! Single precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe cartesian unit 3-vector.
!! @param[out] dxedlat derivative, d(xe)/d(dlat).
!! @param[out] dxedlon derivative, d(xe)/d(dlon).
!! @author R. J. Purser
subroutine sgtocd(dlat,dlon,xe,dxedlat,dxedlon)!                        [gtoc]
implicit none
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
real(dp)             :: dlat_d,dlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d
dlat_d=dlat; dlon_d=dlon
call dgtocd(dlat_d,dlon_d,xe_d,dxedlat_d,dxedlon_d)
xe     =xe_d
dxedlat=dxedlat_d
dxedlon=dxedlon_d
end subroutine sgtocd

!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Also, return the partial derivatives of xe wrt latitude and longitude.
!! Double precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe cartesian unit 3-vector.
!! @param[out] dxedlat derivative, d(xe)/d(dlat).
!! @param[out] dxedlon derivative, d(xe)/d(dlon).
!! @author R. J. Purser
subroutine dgtocd(dlat,dlon,xe,dxedlat,dxedlon)!                        [gtoc]
use pietc, only: u0,dtor
implicit none
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla; dxedlat=dxedlat*dtor
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=u0 ; dxedlon=dxedlon*dtor
end subroutine dgtocd

!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Also, return the 1st and 2nd partial derivatives of xe wrt latitude and 
!! longitude.
!! Single precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe cartesian unit 3-vector.
!! @param[out] dxedlat derivative, d(xe)/d(dlat).
!! @param[out] dxedlon derivative, d(xe)/d(dlon).
!! @param[out] ddxedlatdlat derivative, d^2(xe)/(d(dlat)d(dlat)).
!! @param[out] ddxedlatdlon derivative, d^2(xe)/(d(dlat)d(dlon)).
!! @param[out] ddxedlondlon derivative, d^2(xe)/(d(dlon)d(dlon)).
!! @author R. J. Purser
subroutine sgtocdd(dlat,dlon,xe,dxedlat,dxedlon, &
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)!                           [gtoc]
implicit none
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
real(dp)             :: dlat_d,dlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d, &
                        ddxedlatdlat_d,ddxedlatdlon_d,ddxedlondlon_d
dlat_d=dlat; dlon_d=dlon
call dgtocdd(dlat_d,dlon_d,xe_d,dxedlat_d,dxedlon_d, &
     ddxedlatdlat_d,ddxedlatdlon_d,ddxedlondlon_d)
xe          =xe_d
dxedlat     =dxedlat_d
dxedlon     =dxedlon_d
ddxedlatdlat=ddxedlatdlat_d
ddxedlatdlon=ddxedlatdlon_d
ddxedlondlon=ddxedlondlon_d
end subroutine sgtocdd

!> Transform "Geographical" to "Cartesian" coordinates, where the
!! geographical coordinates refer to latitude and longitude (degrees)
!! and cartesian coordinates are standard earth-centered cartesian
!! coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!! Also, return the 1st and 2nd partial derivatives of xe wrt latitude and 
!! longitude.
!! Double precision version.
!!
!! @param[in] dlat degrees latitude.
!! @param[in] dlon degrees longitude.
!! @param[out] xe cartesian unit 3-vector.
!! @param[out] dxedlat d(xe)/d(dlat).
!! @param[out] dxedlon d(xe)/d(dlon).
!! @param[out] ddxedlatdlat derivative, d^2(xe)/(d(dlat)d(dlat)).
!! @param[out] ddxedlatdlon derivative, d^2(xe)/(d(dlat)d(dlon)).
!! @param[out] ddxedlondlon derivative, d^2(xe)/(d(dlon)d(dlon)).
!! @author R. J. Purser
subroutine dgtocdd(dlat,dlon,xe,dxedlat,dxedlon, &
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)!                           [gtoc]
use pietc, only: u0,dtor
implicit none
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla; dxedlat=dxedlat*dtor
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=u0 ; dxedlon=dxedlon*dtor
ddxedlatdlat(1)=-cla*clo
ddxedlatdlat(2)=-cla*slo
ddxedlatdlat(3)=-sla
ddxedlatdlon(1)= sla*slo
ddxedlatdlon(2)=-sla*clo
ddxedlatdlon(3)= u0
ddxedlondlon(1)=-cla*clo
ddxedlondlon(2)=-cla*slo
ddxedlondlon(3)= u0
ddxedlatdlat=ddxedlatdlat*dtor**2
ddxedlatdlon=ddxedlatdlon*dtor**2
ddxedlondlon=ddxedlondlon*dtor**2
end subroutine dgtocdd

!> From the degree lat and lon (plat and plon) return the standard orthogonal
!! 3D frame at this location as an orthonormal matrix, orth.
!! Single precision version.
!!
!! @param[in] splat latitude (degrees) of point.
!! @param[in] splon longitude (degrees) of point.
!! @param[out] sorth orthonormal matrix.
!! @author R. J. Purser
subroutine sgtoframem(splat,splon,sorth)!                            [gtoframe]
implicit none
real(sp),               intent(in ):: splat,splon
real(sp),dimension(3,3),intent(out):: sorth
real(dp):: plat,plon
real(dp),dimension(3,3):: orth
plat=splat; plon=splon; call gtoframem(plat,plon,orth); sorth=orth
end subroutine sgtoframem

!> From the degree lat and lon (plat and plon) return the standard orthogonal
!! 3D frame at this location as an orthonormal matrix, orth.
!! Double precision version.
!!
!! @param[in] plat latitude (degrees) of point.
!! @param[in] plon longitude (degrees) of point.
!! @param[out] orth orthonormal matrix.
!! @author R. J. Purser
subroutine gtoframem(plat,plon,orth)!                                [gtoframe]
implicit none
real(dp),               intent(in ):: plat,plon
real(dp),dimension(3,3),intent(out):: orth
real(dp),dimension(3):: xp,yp,zp
call gtoframev(plat,plon, xp,yp,zp)
orth(:,1)=xp; orth(:,2)=yp; orth(:,3)=zp
end subroutine gtoframem

!> Given a geographical point by its degrees lat and lon, plat and
!! plon, return its standard orthogonal cartesian frame, {xp,yp,zp} in
!! earth-centered coordinates.
!! Single precision version.
!!
!! @param[in] splat latitude (degrees) of point.
!! @param[in] splon longitude (degrees) of point.
!! @param[out] sxp xp unit X-basis vector of cartesian frame.
!! @param[out] syp yp unit Y-basis vector of cartesian frame.
!! @param[out] szp zp unit Z-basis vector of cartesian frame.
!! @author R. J. Purser
subroutine sgtoframev(splat,splon,sxp,syp,szp)!                       [gtoframe]
!==============================================================================
implicit none
real(sp),             intent(in ):: splat,splon
real(sp),dimension(3),intent(out):: sxp,syp,szp
!------------------------------------------------------------------------------
real(dp)             :: plat,plon
real(dp),dimension(3):: xp,yp,zp
!==============================================================================
plat=splat; plon=splon
call gtoframev(plat,plon, xp,yp,zp)
sxp=xp; syp=yp; szp=zp
end subroutine sgtoframev

!> Given a geographical point by its degrees lat and lon, plat and
!! plon, return its standard orthogonal cartesian frame, {xp,yp,zp} in
!! earth-centered coordinates.
!! Double precision version.
!!
!! @param[in] plat latitude (degrees) of point.
!! @param[in] plon longitude (degrees) of point.
!! @param[out] xp unit X-basis vector of cartesian frame.
!! @param[out] yp unit Y-basis vector of cartesian frame.
!! @param[out] zp unit Z-basis vector of cartesian frame.
!! @author R. J. Purser
subroutine gtoframev(plat,plon, xp,yp,zp)!                           [gtoframe]
use pietc, only: u0,u1
implicit none
real(dp),             intent(in ):: plat,plon
real(dp),dimension(3),intent(out):: xp,yp,zp
real(dp),dimension(3):: dzpdlat,dzpdlon
if(plat==90)then ! is this the north pole?
   xp=(/ u0,u1, u0/) ! Assume the limiting case lat-->90 along the 0-meridian
   yp=(/-u1,u0, u0/) !
   zp=(/ u0,u0, u1/)
elseif(plat==-90)then
   xp=(/ u0,u1, u0/) ! Assume the limiting case lat-->90 along the 0-meridian
   yp=(/ u1,u0, u0/) !
   zp=(/ u0,u0,-u1/)
else
   call gtoc(plat,plon,zp,dzpdlat,dzpdlon)
   xp=dzpdlon/sqrt(dot_product(dzpdlon,dzpdlon))
   yp=dzpdlat/sqrt(dot_product(dzpdlat,dzpdlat))
endif
end subroutine gtoframev

!> Take a principal reference orthonormal frame, {xp,yp,zp} and a dependent
!! point defined by unit vector, zv, and complete the V-frame cartesian
!! components, {xv,yv}, that are the result of parallel-transport of {xp,yp}
!! along the geodesic between P and V.
!! Single precision version.
!!
!! @param[in] sxp reference orthonormal P-frame cartesian X-vector.
!! @param[in] syp reference orthonormal P-frame cartesian Y-vector.
!! @param[in] szp reference orthonormal P-frame cartesian Z-vector.
!! @param[out] sxv V-frame cartesian X-vector.
!! @param[out] syv V-frame cartesian Y-vector.
!! @param[in] szv  dependent point zenith, V-frame cartesian Z-vector.
!! @author R. J. Purser
subroutine sparaframe(sxp,syp,szp, sxv,syv,szv)!                    [paraframe]
implicit none
real(sp),dimension(3),intent(in ):: sxp,syp,szp, szv
real(sp),dimension(3),intent(out):: sxv,syv
real(dp),dimension(3):: xp,yp,zp, xv,yv,zv
xp=sxp; yp=syp; zp=szp
call paraframe(xp,yp,zp, xv,yv,zv)
sxv=xv; syv=yv
end subroutine sparaframe

!> Take a principal reference orthonormal frame, {xp,yp,zp} and a dependent
!! point defined by unit vector, zv, and complete the V-frame cartesian
!! components, {xv,yv}, that are the result of parallel-transport of {xp,yp}
!! along the geodesic between P and V.
!! Double precision version.
!!
!! @param[in] xp reference orthonormal P-frame cartesian X-vector.
!! @param[in] yp reference orthonormal P-frame cartesian Y-vector.
!! @param[in] zp reference orthonormal P-frame cartesian Z-vector.
!! @param[out] xv V-frame cartesian X-vector.
!! @param[out] yv V-frame cartesian Y-vector.
!! @param[in] zv dependent point zenith, V-frame cartesian Z-vector.
!! @author R. J. Purser
subroutine paraframe(xp,yp,zp, xv,yv,zv)!                           [paraframe]
use pmat4,  only: cross_product,normalized
implicit none
real(dp),dimension(3),intent(in ):: xp,yp,zp, zv
real(dp),dimension(3),intent(out):: xv,yv
real(dp)             :: xpofv,ypofv,theta,ctheta,stheta
real(dp),dimension(3):: xq,yq
xpofv=dot_product(xp,zv)
ypofv=dot_product(yp,zv)
theta=atan2(ypofv,xpofv); ctheta=cos(theta); stheta=sin(theta)
xq=zv-zp; xq=xq-zv*dot_product(xq,zv); xq=normalized(xq)
yq=cross_product(zv,xq)
xv=xq*ctheta-yq*stheta
yv=xq*stheta+yq*ctheta
end subroutine paraframe

!> Given a principal cartesian orthonormal frame, {xp,yp,zp} (i.e., at P with
!! Earth-centered cartesians, zp), and another similar frame {xv,yv,zv} at V
!! with Earth-centered cartesians zv, find the relative rotation angle, "twist"
!! by which the frame at V is rotated in the counterclockwise sense relative
!! to the parallel-transportation of P's frame to V.
!! Note that, by symmetry, transposing P and V leads to the opposite twist.
!! Single precision version.
!!
!! @param[in] sxp P-frame cartesian X-vector.
!! @param[in] syp P-frame cartesian Y-vector.
!! @param[in] szp P-frame cartesian Z-vector.
!! @param[in] sxv V-frame cartesian X-vector.
!! @param[in] syv V-frame cartesian Y-vector.
!! @param[in] szv V-frame cartesian Z-vector.
!! @param[out] stwist relative rotation angle (radians) of frames.
!! @author R. J. Purser
subroutine sframetwist(sxp,syp,szp, sxv,syv,szv, stwist)!          [frametwist]
implicit none
real(sp),dimension(3),intent(in ):: sxp,syp,szp, sxv,syv,szv
real(sp),             intent(out):: stwist
real(dp),dimension(3):: xp,yp,zp, xv,yv,zv
real(dp)             :: twist
xp=sxp;yp=syp; zp=szp; xv=sxv; yv=syv; zv=szv
call frametwist(xp,yp,zp, xv,yv,zv, twist)
stwist=twist
end subroutine sframetwist

!> Given a principal cartesian orthonormal frame, {xp,yp,zp} (i.e., at P with
!! Earth-centered cartesians, zp), and another similar frame {xv,yv,zv} at V
!! with Earth-centered cartesians zv, find the relative rotation angle, "twist"
!! by which the frame at V is rotated in the counterclockwise sense relative
!! to the parallel-transportation of P's frame to V.
!! Note that, by symmetry, transposing P and V leads to the opposite twist.
!! Double precision version.
!!
!! @param[in] xp P-frame cartesian X-vector.
!! @param[in] yp P-frame cartesian Y-vector.
!! @param[in] zp P-frame cartesian Z-vector.
!! @param[in] xv V-frame cartesian X-vector.
!! @param[in] yv V-frame cartesian Y-vector.
!! @param[in] zv V-frame cartesian Z-vector.
!! @param[out] twist relative rotation angle (radians) of frames.
!! @author R. J. Purser
subroutine frametwist(xp,yp,zp, xv,yv,zv, twist)!                  [frametwist]
implicit none
real(dp),dimension(3),intent(in ):: xp,yp,zp, xv,yv,zv
real(dp),             intent(out):: twist
real(dp),dimension(3):: xxv,yyv
real(dp)             :: c,s
call paraframe(xp,yp,zp, xxv,yyv,zv)
c=dot_product(xv,xxv); s=dot_product(xv,yyv)
twist=atan2(s,c)
end subroutine frametwist

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s.
!!
!! @param[in] s Schmidt scaling parameter 
!! @param[inout] xc1 input cartesian 3-vector
!! @param[inout] xc2 output cartesian unit 3-vector
!! @author R. J. Purser
subroutine sctoc(s,xc1,xc2)!                                       [ctoc_schm]
use pietc_s, only: u1,u2
implicit none
real(sp),             intent(IN   ):: s
real(sp),dimension(3),intent(INOUT):: xc1,xc2
real(sp)                           :: x,y,z,a,b,d,e,ab2,aa,bb,di,aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
end subroutine sctoc

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s,
!! and its jacobian, dxc2.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[inout] xc1 input cartesian 3-vector.
!! @param[inout] xc2 output cartesian unit 3-vector.
!! @param[out] dxc2 jacobian 1st derivative, d(xc2)/d(xc1).
!! @author R. J. Purser
subroutine sctocd(s,xc1,xc2,dxc2)!                                 [ctoc_schm]
use pietc_s, only: u0,u1,u2
implicit none
real(sp),intent(IN)                  :: s
real(sp),dimension(3),  intent(INOUT):: xc1,xc2
real(sp),dimension(3,3),intent(  OUT):: dxc2
real(sp)                             :: x,y,z,a,b,d,e, &
                                        ab2,aa,bb,di,ddi,aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di

dxc2=u0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi
end subroutine sctocd

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s,
!! its jacobian, dxc2, and its 2nd derivative, ddxc2.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[in] xc1 input cartesian 3-vector.
!! @param[in] xc2 output cartesian unit 3-vector.
!! @param[out] dxc2 jacobian 1st derivative, d(xc2)/d(xc1).
!! @param[out] ddxc2 2nd derivative, d^2(xc2)/(d(xc1)d(xc1)).
!! @author R. J. Purser
subroutine sctocdd(s,xc1,xc2,dxc2,ddxc2)!                          [ctoc_schm]
use pietc_s, only: u0,u1,u2
implicit none
real(sp),                 intent(IN   ):: s 
real(sp),dimension(3),    intent(INOUT):: xc1,xc2
real(sp),dimension(3,3),  intent(  OUT):: dxc2
real(sp),dimension(3,3,3),intent(  OUT):: ddxc2
real(sp)                               :: x,y,z,a,b,d,e, &
                                          ab2,aa,bb,di,ddi,dddi, &
                                          aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di
dddi=ddi*di

dxc2=u0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi

ddxc2=u0
ddxc2(1,1,3)=ab2*aambb*ddi
ddxc2(1,3,1)=ddxc2(1,1,3)
ddxc2(1,3,3)=u2*ab2**2*aambb*x*dddi
ddxc2(2,2,3)=ab2*aambb*ddi
ddxc2(2,3,2)=ddxc2(2,2,3)
ddxc2(2,3,3)=u2*ab2**2*aambb*y*dddi
ddxc2(3,3,3)=u2*ab2*(aapbb*ddi+ab2*e*dddi)
end subroutine sctocdd

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[inout] xc1 input cartesian 3-vector.
!! @param[inout] xc2 output cartesian unit 3-vector.
!! @author R. J. Purser
subroutine dctoc(s,xc1,xc2)!                                       [ctoc_schm]
use pietc, only: u1,u2
implicit none
real(dp),             intent(IN   ):: s
real(dp),dimension(3),intent(INOUT):: xc1,xc2
real(dp)                           :: x,y,z,a,b,d,e, &
                                      ab2,aa,bb,di,aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
end subroutine dctoc

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s,
!! and its jacobian, dxc2.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[inout] xc1 input cartesian 3-vector.
!! @param[inout] xc2 output cartesian unit 3-vector.
!! @param[out] dxc2 jacobian 1st derivative, d(xc2)/d(xc1).
!! @author R. J. Purser
subroutine dctocd(s,xc1,xc2,dxc2)!                                 [ctoc_schm]
use pietc, only: u0,u1,u2
implicit none
real(dp),               intent(IN   ):: s
real(dp),dimension(3),  intent(INOUT):: xc1,xc2
real(dp),dimension(3,3),intent(  OUT):: dxc2 
real(dp)                             :: x,y,z,a,b,d,e, &
                                        ab2,aa,bb,di,ddi,aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di

dxc2=u0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi
end subroutine dctocd

!> Evaluate Schmidt transformation, xc1 --> xc2, with scaling parameter s,
!! its jacobian, dxc2, and its 2nd derivative, ddxc2.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[inout] xc1 input cartesian 3-vector.
!! @param[inout] xc2 output cartesian unit 3-vector.
!! @param[out] dxc2 jacobian 1st derivative, d(xc2)/d(xc1).
!! @param[out] ddxc2 2nd derivative, d^2(xc2)/(d(xc1)d(xc1)).
!! @author R. J. Purser
subroutine dctocdd(s,xc1,xc2,dxc2,ddxc2)!                          [ctoc_schm]
use pietc, only: u0,u1,u2
implicit none
real(dp),intent(IN)                    :: s
real(dp),dimension(3),    intent(INOUT):: xc1,xc2
real(dp),dimension(3,3),  intent(OUT  ):: dxc2
real(dp),dimension(3,3,3),intent(OUT  ):: ddxc2
real(dp)                               :: x,y,z,a,b,d,e, &
                                          ab2,aa,bb,di,ddi,dddi, &
                                          aapbb,aambb
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+u1
b=s-u1
ab2=a*b*u2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=u1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di
dddi=ddi*di

dxc2=u0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi

ddxc2=u0
ddxc2(1,1,3)=ab2*aambb*ddi
ddxc2(1,3,1)=ddxc2(1,1,3)
ddxc2(1,3,3)=u2*ab2**2*aambb*x*dddi
ddxc2(2,2,3)=ab2*aambb*ddi
ddxc2(2,3,2)=ddxc2(2,2,3)
ddxc2(2,3,3)=u2*ab2**2*aambb*y*dddi
ddxc2(3,3,3)=u2*ab2*(aapbb*ddi+ab2*e*dddi)
end subroutine dctocdd

!> Apply a constant rotation to a three dimensional polyline.
!!
!! @param[in] rot3 rotation matrix.
!! @param[in] n number of points in the polyline.
!! @param[inout] x cartesian components of the three dimensional polyline.
!! @param[inout] y cartesian components of the three dimensional polyline.
!! @param[inout] z cartesian components of the three dimensional polyline.
!! @author R. J. Purser
subroutine plrot(rot3,n,x,y,z)!                                        [plrot]
implicit none
integer,                intent(IN   ):: n
real(sp),dimension(3,3),intent(IN   ):: rot3
real(sp),dimension(n),  intent(INOUT):: x,y,z
real(sp),dimension(3)                :: t
integer                              :: k
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(rot3,t) 
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine plrot

!> Invert the rotation of a three-dimensional polyline.
!!
!! @param[in] rot3 rotation to be inverted.
!! @param[in] n number of points in the polyline.
!! @param[inout] x cartesian components of the three dimensional polyline.
!! @param[inout] y cartesian components of the three dimensional polyline.
!! @param[inout] z cartesian components of the three dimensional polyline.
!! @author R. J. Purser
subroutine plroti(rot3,n,x,y,z)!                                      [plroti]
implicit none
integer,                intent(IN   ):: n
real(sp),dimension(3,3),intent(IN   ):: rot3
real(sp),dimension(n),  intent(INOUT):: x,y,z
real(sp),dimension(3)                :: t
integer                              :: k
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(t,rot3)
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine plroti

!> Apply a constant rotation to a three dimensional polyline.
!!
!! @param[in] rot3 rotation matrix.
!! @param[in] n number of points in the polyline.
!! @param[inout] x cartesian components of the three dimensional polyline.
!! @param[inout] y cartesian components of the three dimensional polyline.
!! @param[inout] z cartesian components of the three dimensional polyline.
!! @author R. J. Purser
subroutine dplrot(rot3,n,x,y,z)!                                        [plrot]
implicit none
integer,                intent(IN   ):: n
real(dP),dimension(3,3),intent(IN   ):: rot3
real(dP),dimension(n),  intent(INOUT):: x,y,z
real(dP),dimension(3)                :: t
integer                              :: k
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(rot3,t) 
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine dplrot

!> Invert the rotation of a three-dimensional polyline.
!!
!! @param[in] rot3 rotation to be inverted.
!! @param[in] n number of points in the polyline.
!! @param[inout] x cartesian components of the three dimensional polyline.
!! @param[inout] y cartesian components of the three dimensional polyline.
!! @param[inout] z cartesian components of the three dimensional polyline.
!! @author R. J. Purser
subroutine dplroti(rot3,n,x,y,z)!                                      [plroti]
implicit none
integer,                intent(IN   ):: n
real(dP),dimension(3,3),intent(IN   ):: rot3
real(dP),dimension(n),  intent(INOUT):: x,y,z
real(dP),dimension(3)                :: t
integer                              :: k
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(t,rot3)
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine dplroti

!> Perform Schmidt transformation with scaling parameter s to a polyline.
!!
!! @param[in] s Schmidt scaling parameter.
!! @param[in] n number of points in the polyline.
!! @param[inout] x cartesian components of the three dimensional polyline.
!! @param[inout] y cartesian components of the three dimensional polyline.
!! @param[inout] z cartesian components of the three dimensional polyline.
!! @author R. J. Purser
subroutine plctoc(s,n,x,y,z)!                                         [plctoc]
use pietc_s, only: u1
implicit none
integer,              intent(IN   ):: n
real(sp),             intent(IN   ):: s
real(sp),dimension(n),intent(INOUT):: x,y,z
real(sp)                           :: a,b,d,e,ab2,aa,bb,di,aapbb,aambb
integer                            :: i
a=s+u1
b=s-u1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
do i=1,n
   d=aapbb-ab2*z(i)
   e=aapbb*z(i)-ab2
   di=1/d
   x(i)=(aambb*x(i))*di
   y(i)=(aambb*y(i))*di
   z(i)=e*di
enddo
end subroutine plctoc

end module pmat5



