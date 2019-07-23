!
!                                **********************************************
!                                *  MODULES cstgeo, dcstgeo, pmat5            *
!                                *  R. J. Purser, NOAA/NCEP/EMC          1996 * 
!                                *                              26 Sep   2012 *
!                                *  jim.purser@noaa.gov                       *
!                                *                                            *
!                                **********************************************
! Handy geographical transformations
!
! DEPENDENCIES
! Modules:  pkind, pietc, pmat4
!=============================================================================
module cstgeo ! Constants for orientation and stretching of map
!=============================================================================
use pkind, only: sp
implicit none
real(sp),dimension(3,3):: rotm
real(sp)               :: sc,sci
!=============================================================================
end module cstgeo

!=============================================================================
module dcstgeo ! Constants for orientation and stretching of map
!=============================================================================
use pkind, only: dp
implicit none
real(dp),dimension(3,3):: rotm
real(dp)               :: sc,sci
!=============================================================================
end module dcstgeo

! Utility routines for orienting the globe and basic geographical mappings
!=============================================================================
module pmat5
!=============================================================================
use pkind, only: sp,dp
implicit none
private
public :: ininmap,inivmap,ctog,gtoc,gtoframe,paraframe,frametwist,&
          ctoc_schm,plrot,plroti,plctoc
interface ininmap;   module procedure sininmap,dininmap;          end interface
interface inivmap;   module procedure sinivmap,dinivmap;          end interface
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
!=============================================================================
subroutine sininmap(alon0,alat0,rot3)!                               [ininmap]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1995	       
!  Initialize the rotation matrix ROT3 needed to transform standard	       
!  earth-centered cartesian components to the alternative cartesian frame      
!  oriented so as to put geographical point (ALAT0,ALON0) on the projection    
!  axis.								       
!=============================================================================
use pietc, only: dtor
real(sp),               intent(IN ):: alon0,alat0
real(sp),dimension(3,3),intent(OUT):: rot3
!-----------------------------------------------------------------------------
real(sp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
!=============================================================================
blon0=dtor*alon0; clon0=cos(blon0); slon0=sin(blon0)
blat0=dtor*alat0; clat0=cos(blat0); slat0=sin(blat0)

rot3(1,1)=slat0*clon0; rot3(1,2)=slat0*slon0; rot3(1,3)=-clat0
rot3(2,1)=-slon0;      rot3(2,2)=clon0;       rot3(2,3)=0
rot3(3,1)=clat0*clon0; rot3(3,2)=clat0*slon0; rot3(3,3)=slat0
end subroutine sininmap
!=============================================================================
subroutine dininmap(alon0,alat0,rot3)!                               [ininmap]
!=============================================================================
use pietc, only: dtor
real(dp),               intent(IN ):: alon0,alat0
real(dp),dimension(3,3),intent(OUT):: rot3
!-----------------------------------------------------------------------------
real(dp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
!=============================================================================
blon0=dtor*alon0; clon0=cos(blon0); slon0=sin(blon0)
blat0=dtor*alat0; clat0=cos(blat0); slat0=sin(blat0)

rot3(1,1)=slat0*clon0; rot3(1,2)=slat0*slon0; rot3(1,3)=-clat0
rot3(2,1)=-slon0;      rot3(2,2)=clon0;       rot3(2,3)=0
rot3(3,1)=clat0*clon0; rot3(3,2)=clat0*slon0; rot3(3,3)=slat0
end subroutine dininmap

!=============================================================================
subroutine sinivmap(alon0,alat0,rot3)!                               [inivmap]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1995	       
!  Initialize the rotation matrix ROT3 needed to transform standard	       
!  earth-centered cartesian components to the alternative cartesian frame      
!  oriented so as to put geographical point (ALAT0,ALON0) at the viewing
!  nadir.								       
!=============================================================================
use pietc, only: dtor
real(sp),               intent(IN ):: alon0,alat0
real(sp),dimension(3,3),intent(OUT):: rot3
!-----------------------------------------------------------------------------
real(sp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
!=============================================================================
blon0=dtor*alon0
blat0=dtor*alat0
clon0=cos(blon0)
slon0=sin(blon0)
clat0=cos(blat0)
slat0=sin(blat0)
rot3(1,1)=-slon0
rot3(1,2)=clon0
rot3(1,3)=0
rot3(2,1)=-slat0*clon0
rot3(2,2)=-slat0*slon0
rot3(2,3)=clat0
rot3(3,1)=clat0*clon0
rot3(3,2)=clat0*slon0
rot3(3,3)=slat0
end subroutine sinivmap
!=============================================================================
subroutine dinivmap(alon0,alat0,rot3)!                               [inivmap]
!=============================================================================
use pietc, only: dtor
real(dp),               intent(IN ):: alon0,alat0
real(dp),dimension(3,3),intent(OUT):: rot3
!-----------------------------------------------------------------------------
real(dp)                           :: blon0,blat0,clon0,clat0,slon0,slat0
!=============================================================================
blon0=dtor*alon0
blat0=dtor*alat0
clon0=cos(blon0)
slon0=sin(blon0)
clat0=cos(blat0)
slat0=sin(blat0)
rot3(1,1)=-slon0
rot3(1,2)=clon0
rot3(1,3)=0
rot3(2,1)=-slat0*clon0
rot3(2,2)=-slat0*slon0
rot3(2,3)=clat0
rot3(3,1)=clat0*clon0
rot3(3,2)=clat0*slon0
rot3(3,3)=slat0
end subroutine dinivmap

!=============================================================================
subroutine sctog(xe,dlat,dlon)!                                         [ctog]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994	      
!                   SUBROUTINE CTOG
!   Transform "Cartesian" to "Geographical" coordinates, where the
!   geographical coordinates refer to latitude and longitude (degrees) 
!   and cartesian coordinates are standard earth-centered cartesian    
!   coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!  --> XE	three cartesian components.				      
!  <-- DLAT	degrees latitude					      
!  <-- DLON	degrees longitude					      
!=============================================================================
use pietc, only: u0,rtod
real(sp),dimension(3),intent(IN ):: xe
real(sp),             intent(OUT):: dlat,dlon
!-----------------------------------------------------------------------------
real(sp)                         :: r
!=============================================================================
r=sqrt(xe(1)**2+xe(2)**2)
dlat=atan2(xe(3),r)*rtod
if(r==u0)then
   dlon=u0
else
   dlon=atan2(xe(2),xe(1))*rtod
endif
end subroutine sctog

!=============================================================================
subroutine dctog(xe,dlat,dlon)!                                         [ctog]
!=============================================================================
use pietc, only: u0,rtod
real(dp),dimension(3),intent(IN ):: xe
real(dp),             intent(OUT):: dlat,dlon
!-----------------------------------------------------------------------------
real(dp)                         :: r
!=============================================================================
r=sqrt(xe(1)**2+xe(2)**2)
dlat=atan2(xe(3),r)*rtod
if(r==u0)then
   dlon=u0
else
   dlon=atan2(xe(2),xe(1))*rtod
endif
end subroutine dctog

!=============================================================================
subroutine sgtoc(dlat,dlon,xe)!                                         [gtoc]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994	      
!		    SUBROUTINE GTOC					      
!   Transform "Geographical" to "Cartesian" coordinates, where the
!   geographical coordinates refer to latitude and longitude (degrees) 
!   and cartesian coordinates are standard earth-centered cartesian    
!   coordinates: xe(3) pointing north, xe(1) pointing to the 0-meridian.
!  --> DLAT	degrees latitude					      
!  --> DLON	degrees longitude					      
!  <-- XE	three cartesian components.				      
!=============================================================================
use pietc, only: dtor
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe
!-----------------------------------------------------------------------------
real(sp)                         :: rlat,rlon,sla,cla,slo,clo
!=============================================================================
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine sgtoc
!=============================================================================
subroutine dgtoc(dlat,dlon,xe)!                                         [gtoc]
!=============================================================================
use pietc, only: dtor
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe
!-----------------------------------------------------------------------------
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
!=============================================================================
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine dgtoc
!=============================================================================
subroutine sgtocd(dlat,dlon,xe,dxedlat,dxedlon)!                        [gtoc]
!=============================================================================
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
!-----------------------------------------------------------------------------
real(dp)             :: dlat_d,dlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d
!=============================================================================
dlat_d=dlat; dlon_d=dlon
call dgtocd(dlat_d,dlon_d,xe_d,dxedlat_d,dxedlon_d)
xe     =xe_d
dxedlat=dxedlat_d
dxedlon=dxedlon_d
end subroutine sgtocd
!=============================================================================
subroutine dgtocd(dlat,dlon,xe,dxedlat,dxedlon)!                        [gtoc]
!=============================================================================
use pietc, only: dtor
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon
!-----------------------------------------------------------------------------
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
!=============================================================================
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla; dxedlat=dxedlat*dtor
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=0  ; dxedlon=dxedlon*dtor
end subroutine dgtocd
!=============================================================================
subroutine sgtocdd(dlat,dlon,xe,dxedlat,dxedlon, &
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)!                           [gtoc]
!=============================================================================
use pietc, only: dtor
real(sp),             intent(IN ):: dlat,dlon
real(sp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
!-----------------------------------------------------------------------------
real(dp)             :: dlat_d,dlon_d
real(dp),dimension(3):: xe_d,dxedlat_d,dxedlon_d, &
                        ddxedlatdlat_d,ddxedlatdlon_d,ddxedlondlon_d
!=============================================================================
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
!=============================================================================
subroutine dgtocdd(dlat,dlon,xe,dxedlat,dxedlon, &
     ddxedlatdlat,ddxedlatdlon,ddxedlondlon)!                           [gtoc]
!=============================================================================
use pietc, only: dtor
real(dp),             intent(IN ):: dlat,dlon
real(dp),dimension(3),intent(OUT):: xe,dxedlat,dxedlon, &
                                    ddxedlatdlat,ddxedlatdlon,ddxedlondlon
!-----------------------------------------------------------------------------
real(dp)                         :: rlat,rlon,sla,cla,slo,clo
!=============================================================================
rlat=dtor*dlat; rlon=dtor*dlon
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
dxedlat(1)=-sla*clo; dxedlat(2)=-sla*slo; dxedlat(3)=cla; dxedlat=dxedlat*dtor
dxedlon(1)=-cla*slo; dxedlon(2)= cla*clo; dxedlon(3)=0  ; dxedlon=dxedlon*dtor
ddxedlatdlat(1)=-cla*clo
ddxedlatdlat(2)=-cla*slo
ddxedlatdlat(3)=-sla
ddxedlatdlon(1)= sla*slo
ddxedlatdlon(2)=-sla*clo
ddxedlatdlon(3)= 0
ddxedlondlon(1)=-cla*clo
ddxedlondlon(2)=-cla*slo
ddxedlondlon(3)= 0
ddxedlatdlat=ddxedlatdlat*dtor**2
ddxedlatdlon=ddxedlatdlon*dtor**2
ddxedlondlon=ddxedlondlon*dtor**2
end subroutine dgtocdd

!==============================================================================
subroutine sgtoframem(splat,splon,sorth)!                            [gtoframe]
!==============================================================================
real(sp),               intent(in ):: splat,splon
real(sp),dimension(3,3),intent(out):: sorth
!------------------------------------------------------------------------------
real(dp):: plat,plon
real(dp),dimension(3,3):: orth
!==============================================================================
plat=splat; plon=splon; call gtoframem(plat,plon,orth); sorth=orth
end subroutine sgtoframem
!==============================================================================
subroutine gtoframem(plat,plon,orth)!                                [gtoframe]
!==============================================================================
! From the degree lat and lo (plat and plon) return the standard orthogonal
! 3D frame at this location as an orthonormal matrix, orth.
!==============================================================================
real(dp),               intent(in ):: plat,plon
real(dp),dimension(3,3),intent(out):: orth
!------------------------------------------------------------------------------
real(dp),dimension(3):: xp,yp,zp
!==============================================================================
call gtoframev(plat,plon, xp,yp,zp)
orth(:,1)=xp; orth(:,2)=yp; orth(:,3)=zp
end subroutine gtoframem
!==============================================================================
subroutine sgtoframev(splat,splon,sxp,syp,szp)!                       [gtoframe]
!==============================================================================
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
!==============================================================================
subroutine gtoframev(plat,plon, xp,yp,zp)!                           [gtoframe]
!==============================================================================
! Given a geographical point by its degrees lat and lon, plat and plon,
! return its standard orthogonal cartesian frame, {xp,yp,zp} in earth-centered
! coordinates.
!=============================================================================
use pietc, only: u0,u1
real(dp),             intent(in ):: plat,plon
real(dp),dimension(3),intent(out):: xp,yp,zp
!------------------------------------------------------------------------------
real(dp),dimension(3):: dzpdlat,dzpdlon
!==============================================================================
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

!==============================================================================
subroutine sparaframe(sxp,syp,szp, sxv,syv,szv)!                    [paraframe]
!==============================================================================
real(sp),dimension(3),intent(in ):: sxp,syp,szp, szv
real(sp),dimension(3),intent(out):: sxv,syv
!-----------------------------------------------------------------------------
real(dp),dimension(3):: xp,yp,zp, xv,yv,zv
!=============================================================================
xp=sxp; yp=syp; zp=szp
call paraframe(xp,yp,zp, xv,yv,zv)
sxv=xv; syv=yv
end subroutine sparaframe
!==============================================================================
subroutine paraframe(xp,yp,zp, xv,yv,zv)!                           [paraframe]
!==============================================================================
! Take a principal reference orthonormal frame, {xp,yp,zp} and a dependent
! point defined by unit vector, zv, and complete the V-frame cartesian
! components, {xv,yv}, that are the result of parallel-transport of {xp,yp}
! along the geodesic between P and V
!==============================================================================
use pmat4,  only: cross_product,normalized
real(dp),dimension(3),intent(in ):: xp,yp,zp, zv
real(dp),dimension(3),intent(out):: xv,yv
!------------------------------------------------------------------------------
real(dp)             :: xpofv,ypofv,theta,ctheta,stheta
real(dp),dimension(3):: xq,yq
!==============================================================================
xpofv=dot_product(xp,zv)
ypofv=dot_product(yp,zv)
theta=atan2(ypofv,xpofv); ctheta=cos(theta); stheta=sin(theta)
xq=zv-zp; xq=xq-zv*dot_product(xq,zv); xq=normalized(xq)
yq=cross_product(zv,xq)
xv=xq*ctheta-yq*stheta
yv=xq*stheta+yq*ctheta
end subroutine paraframe

!==============================================================================
subroutine sframetwist(sxp,syp,szp, sxv,syv,szv, stwist)!          [frametwist]
!==============================================================================
real(sp),dimension(3),intent(in ):: sxp,syp,szp, sxv,syv,szv
real(sp),             intent(out):: stwist
!------------------------------------------------------------------------------
real(dp),dimension(3):: xp,yp,zp, xv,yv,zv
real(dp)             :: twist
!==============================================================================
xp=sxp;yp=syp; zp=szp; xv=sxv; yv=syv; zv=szv
call frametwist(xp,yp,zp, xv,yv,zv, twist)
stwist=twist
end subroutine sframetwist
!==============================================================================
subroutine frametwist(xp,yp,zp, xv,yv,zv, twist)!                  [frametwist]
!==============================================================================
! Given a principal cartesian orthonormal frame, {xp,yp,zp} (i.e., at P with
! Earth-centered cartesians, zp), and another similar frame {xv,yv,zv} at V
! with Earth-centered cartesians zv, find the relative rotation angle, "twist"
! by which the frame at V is rotated in the counterclockwise sense relative
! to the parallel-transportation of P's frame to V.
! Note that, by symmetry, transposing P and V leads to the opposite twist.
!==============================================================================
real(dp),dimension(3),intent(in ):: xp,yp,zp, xv,yv,zv
real(dp),             intent(out):: twist
!------------------------------------------------------------------------------
real(dp),dimension(3):: xxv,yyv
real(dp)             :: c,s
!==============================================================================
call paraframe(xp,yp,zp, xxv,yyv,zv)
c=dot_product(xv,xxv); s=dot_product(xv,yyv)
twist=atan2(s,c)
end subroutine frametwist

!=============================================================================
subroutine sctoc(s,xc1,xc2)!                                       [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s
!=============================================================================
real(sp),             intent(IN   ):: s
real(sp),dimension(3),intent(INOUT):: xc1,xc2
!-----------------------------------------------------------------------------
real(sp)                           :: x,y,z,a,b,d,e,ab2,aa,bb,di,aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
end subroutine sctoc

!=============================================================================
subroutine sctocd(s,xc1,xc2,dxc2)!                                 [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s,
!  and its jacobian, dxc2.
!=============================================================================
real(sp),intent(IN)                  :: s
real(sp),dimension(3),  intent(INOUT):: xc1,xc2
real(sp),dimension(3,3),intent(  OUT):: dxc2
!-----------------------------------------------------------------------------
real(sp)                             :: x,y,z,a,b,d,e, &
                                        ab2,aa,bb,di,ddi,aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di

dxc2=0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi
end subroutine sctocd

!=============================================================================
subroutine sctocdd(s,xc1,xc2,dxc2,ddxc2)!                          [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s,
!  its jacobian, dxc2, and its 2nd derivative, ddxc2.
!=============================================================================
real(sp),                 intent(IN   ):: s
real(sp),dimension(3),    intent(INOUT):: xc1,xc2
real(sp),dimension(3,3),  intent(  OUT):: dxc2
real(sp),dimension(3,3,3),intent(  OUT):: ddxc2
!-----------------------------------------------------------------------------
real(sp)                               :: x,y,z,a,b,d,e, &
                                          ab2,aa,bb,di,ddi,dddi, &
                                          aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di
dddi=ddi*di

dxc2=0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi

ddxc2=0
ddxc2(1,1,3)=ab2*aambb*ddi
ddxc2(1,3,1)=ddxc2(1,1,3)
ddxc2(1,3,3)=2*ab2**2*aambb*x*dddi
ddxc2(2,2,3)=ab2*aambb*ddi
ddxc2(2,3,2)=ddxc2(2,2,3)
ddxc2(2,3,3)=2*ab2**2*aambb*y*dddi
ddxc2(3,3,3)=2*ab2*(aapbb*ddi+ab2*e*dddi)
end subroutine sctocdd

!=============================================================================
subroutine dctoc(s,xc1,xc2)!                                       [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s
!=============================================================================
real(dp),             intent(IN   ):: s
real(dp),dimension(3),intent(INOUT):: xc1,xc2
!-----------------------------------------------------------------------------
real(dp)                           :: x,y,z,a,b,d,e, &
                                      ab2,aa,bb,di,aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
end subroutine dctoc

!=============================================================================
subroutine dctocd(s,xc1,xc2,dxc2)!                                 [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s,
!  and its jacobian, dxc2.
!=============================================================================
real(dp),               intent(IN   ):: s
real(dp),dimension(3),  intent(INOUT):: xc1,xc2
real(dp),dimension(3,3),intent(  OUT):: dxc2
!-----------------------------------------------------------------------------
real(dp)                             :: x,y,z,a,b,d,e, &
                                        ab2,aa,bb,di,ddi,aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di

dxc2=0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi
end subroutine dctocd

!=============================================================================
subroutine dctocdd(s,xc1,xc2,dxc2,ddxc2)!                          [ctoc_schm]
!=============================================================================
!  Evaluate schmidt transformation, xc1 --> xc2, with scaling parameter s,
!  its jacobian, dxc2, and its 2nd derivative, ddxc2.
!=============================================================================
real(dp),intent(IN)                    :: s
real(dp),dimension(3),    intent(INOUT):: xc1,xc2
real(dp),dimension(3,3),  intent(OUT  ):: dxc2
real(dp),dimension(3,3,3),intent(OUT  ):: ddxc2
!-----------------------------------------------------------------------------
real(dp)                               :: x,y,z,a,b,d,e, &
                                          ab2,aa,bb,di,ddi,dddi, &
                                          aapbb,aambb
!=============================================================================
x=xc1(1); y=xc1(2); z=xc1(3)
a=s+1
b=s-1
ab2=a*b*2
aa=a*a
bb=b*b
aapbb=aa+bb
aambb=aa-bb
d=aapbb-ab2*z
e=aapbb*z-ab2
di=1/d
xc2(1)=(aambb*x)*di
xc2(2)=(aambb*y)*di
xc2(3)=e*di
ddi=di*di
dddi=ddi*di

dxc2=0
dxc2(1,1)=aambb*di
dxc2(1,3)=ab2*aambb*x*ddi
dxc2(2,2)=aambb*di
dxc2(2,3)=ab2*aambb*y*ddi
dxc2(3,3)=aapbb*di +ab2*e*ddi

ddxc2=0
ddxc2(1,1,3)=ab2*aambb*ddi
ddxc2(1,3,1)=ddxc2(1,1,3)
ddxc2(1,3,3)=2*ab2**2*aambb*x*dddi
ddxc2(2,2,3)=ab2*aambb*ddi
ddxc2(2,3,2)=ddxc2(2,2,3)
ddxc2(2,3,3)=2*ab2**2*aambb*y*dddi
ddxc2(3,3,3)=2*ab2*(aapbb*ddi+ab2*e*dddi)
end subroutine dctocdd

!=============================================================================
subroutine plrot(rot3,n,x,y,z)!                                        [plrot]
!=============================================================================
! Apply a constant rotation to a three dimensional polyline
!=============================================================================
integer,                intent(IN   ):: n
real(sp),dimension(3,3),intent(IN   ):: rot3
real(sp),dimension(n),  intent(INOUT):: x,y,z
!-----------------------------------------------------------------------------
real(sp),dimension(3)                :: t
integer                              :: k
!=============================================================================
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(rot3,t) 
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine plrot

!=============================================================================
subroutine plroti(rot3,n,x,y,z)!                                      [plroti]
!=============================================================================
! Invert the rotation of a three-dimensional polyline
!=============================================================================
integer,                intent(IN   ):: n
real(sp),dimension(3,3),intent(IN   ):: rot3
real(sp),dimension(n),  intent(INOUT):: x,y,z
!-----------------------------------------------------------------------------
real(sp),dimension(3)                :: t
integer                              :: k
!=============================================================================
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(t,rot3)
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine plroti

!=============================================================================
subroutine dplrot(rot3,n,x,y,z)!                                        [plrot]
!=============================================================================
! Apply a constant rotation to a three dimensional polyline
!=============================================================================
integer,                intent(IN   ):: n
real(dP),dimension(3,3),intent(IN   ):: rot3
real(dP),dimension(n),  intent(INOUT):: x,y,z
!-----------------------------------------------------------------------------
real(dP),dimension(3)                :: t
integer                              :: k
!=============================================================================
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(rot3,t) 
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine dplrot

!=============================================================================
subroutine dplroti(rot3,n,x,y,z)!                                      [plroti]
!=============================================================================
! Invert the rotation of a three-dimensional polyline
!=============================================================================
integer,                intent(IN   ):: n
real(dP),dimension(3,3),intent(IN   ):: rot3
real(dP),dimension(n),  intent(INOUT):: x,y,z
!-----------------------------------------------------------------------------
real(dP),dimension(3)                :: t
integer                              :: k
!=============================================================================
do k=1,n
   t(1)=x(k); t(2)=y(k); t(3)=z(k)
   t=matmul(t,rot3)
   x(k)=t(1); y(k)=t(2); z(k)=t(3)
enddo
end subroutine dplroti

!=============================================================================
subroutine plctoc(s,n,x,y,z)!                                         [plctoc]
!=============================================================================
!  Perform schmidt transformation with scaling parameter s to a polyline
!=============================================================================
integer,              intent(IN   ):: n
real(sp),             intent(IN   ):: s
real(sp),dimension(n),intent(INOUT):: x,y,z
!-----------------------------------------------------------------------------
real(sp)                           :: a,b,d,e,ab2,aa,bb,di,aapbb,aambb
integer                            :: i
!=============================================================================
a=s+1
b=s-1
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



