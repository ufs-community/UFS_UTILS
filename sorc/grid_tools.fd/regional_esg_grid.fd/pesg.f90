!> @file
!! @brief Routines to perform ESG regional grid mappings.
!! @author R. J. Purser @date May 2020

!> Suite of routines to perform the 2-parameter family of Extended
!! Schmidt Gnomonic (ESG) regional grid mappings, and to optimize the
!! the two parameters, A and K, of those mappings for a given
!! rectangular domain's principal (median) semi-arcs with respect to a
!! domain-averaged measure of distortion. This criterion is itself
!! endowed with a parameter, lam (for "lambda" in [0,1) ) which gives
!! weight to additional weight areal inhomogeneities instead of
!! treating all distortion components equally.
!!
!! @author R. J. Purser
module pesg
!=============================================================================
use pkind, only: spi,dp
use pietc, only: F,T,u0,u1,u2,o2,rtod,dtor,pih,pi2
implicit none
private
public :: xctoxm_ak,xmtoxc_ak,get_edges,bestesg_geo,bestesg_map, &
          hgrid_ak_rr,hgrid_ak_rc,hgrid_ak_dd,hgrid_ak_dc,hgrid_ak, &
          gtoxm_ak_rr,gtoxm_ak_dd,xmtog_ak_rr,xmtog_ak_dd

interface xctoxs;         module procedure xctoxs;                end interface
interface xstoxc;         module procedure xstoxc,xstoxc1;        end interface
interface xstoxt;         module procedure xstoxt;                end interface
interface xttoxs;         module procedure xttoxs,xttoxs1;        end interface
interface xttoxm;         module procedure xttoxm;                end interface
interface zttozm;         module procedure zttozm;                end interface
interface xmtoxt;         module procedure xmtoxt,xmtoxt1;        end interface
interface zmtozt;         module procedure zmtozt,zmtozt1;        end interface
interface xctoxm_ak;      module procedure xctoxm_ak;             end interface
interface xmtoxc_ak
   module procedure xmtoxc_ak,xmtoxc_vak,xmtoxc_vak1;             end interface
interface get_edges;      module procedure get_edges;             end interface
interface get_qx;         module procedure get_qx,get_qxd;        end interface
interface get_qofv;module procedure get_qofv,get_qofvd,get_qsofvs;end interface
interface get_meanq;      module procedure get_meanqd,get_meanqs; end interface
interface guessak_map;    module procedure guessak_map;           end interface
interface guessak_geo;    module procedure guessak_geo;           end interface
interface bestesg_geo;    module procedure bestesg_geo;           end interface
interface bestesg_map;    module procedure bestesg_map;           end interface
interface hgrid_ak_rr;module procedure hgrid_ak_rr,hgrid_ak_rr_c; end interface
interface hgrid_ak_rc;    module procedure hgrid_ak_rc;           end interface
interface hgrid_ak_dd;module procedure hgrid_ak_dd,hgrid_ak_dd_c; end interface
interface hgrid_ak_dc;    module procedure hgrid_ak_dc;           end interface
interface hgrid_ak;       module procedure hgrid_ak,hgrid_ak_c;   end interface
interface gtoxm_ak_rr
   module procedure gtoxm_ak_rr_m,gtoxm_ak_rr_g;                  end interface
interface gtoxm_ak_dd
   module procedure gtoxm_ak_dd_m,gtoxm_ak_dd_g;                  end interface
interface xmtog_ak_rr
   module procedure xmtog_ak_rr_m,xmtog_ak_rr_g;                  end interface
interface xmtog_ak_dd
   module procedure xmtog_ak_dd_m,xmtog_ak_dd_g;                  end interface

interface gaulegh;        module procedure gaulegh;               end interface

contains

!> Inverse of xstoxc. I.e., cartesians to stereographic.
!!
!! @param[in] xc Earth-centered cartesian unit 3-vector
!! @param[out] xs Stereographic map coordinates
!! @author R. J. Purser
subroutine xctoxs(xc,xs)!                                             [xctoxs]
implicit none
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xs
real(dp):: zp
zp=u1+xc(3); xs=xc(1:2)/zp
end subroutine xctoxs

!> Standard transformation from polar stereographic map coordinates,
!! xs, to cartesian, xc, assuming the projection axis is polar.
!! xcd=d(xc)/d(xs) is the jacobian matrix, encoding distortion and
!! metric.
!!
!! @param[in] xs Stereographic map coordinates
!! @param[out] xc Cartesian earth-centered 3-vector
!! @param[out] xcd Value of jacobian matrix, encoding distortion and metric
!! @author R. J. Purser
subroutine xstoxc(xs,xc,xcd)!                                         [xstoxc]
use pmat4, only: outer_product
implicit none
real(dp),dimension(2),  intent(in ):: xs
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
real(dp):: zp
zp=u2/(u1+dot_product(xs,xs)); xc(1:2)=xs*zp; xc(3)=zp
xcd=-outer_product(xc,xs)*zp; xcd(1,1)=xcd(1,1)+zp; xcd(2,2)=xcd(2,2)+zp
xc(3)=xc(3)-u1
end subroutine xstoxc

!> Standard transformation from polar stereographic map coordinates,
!! xs, to cartesian, xc, assuming the projection axis is polar.
!! xcd=d(xc)/d(xs) is the jacobian matrix, encoding distortion and
!! metric. xcdd is the further derivative, wrt xs, of xcd.
!!
!! @param[in] xs Stereographic map coordinates
!! @param[out] xc Cartesian earth-centered 3-vector
!! @param[out] xcd Jacobian matrix, encoding distortion and metric
!! @param[out] xcdd Further derivative, wrt xs, of xcd
!! @author R. J. Purser
subroutine xstoxc1(xs,xc,xcd,xcdd)!                                   [xstoxc]
use pmat4, only: outer_product
implicit none
real(dp),dimension(2),    intent(in ):: xs
real(dp),dimension(3),    intent(out):: xc
real(dp),dimension(3,2),  intent(out):: xcd
real(dp),dimension(3,2,2),intent(out):: xcdd
real(dp),dimension(3,2):: zpxcdxs
real(dp),dimension(3)  :: zpxc
real(dp)               :: zp
integer(spi)           :: i
zp=u2/(u1+dot_product(xs,xs)); xc(1:2)=xs*zp; xc(3)=zp
xcd=-outer_product(xc,xs)*zp
zpxc=zp*xc; xc(3)=xc(3)-u1; xcdd=u0
do i=1,2
   zpxcdxs=xcd*xc(i)
   xcdd(:,i,i)=xcdd(:,i,i)-zpxc
   xcdd(:,i,:)=xcdd(:,i,:)-zpxcdxs
   xcdd(:,:,i)=xcdd(:,:,i)-zpxcdxs
   xcdd(i,:,i)=xcdd(i,:,i)-zpxc(1:2)
   xcdd(i,i,:)=xcdd(i,i,:)-zpxc(1:2)
enddo
do i=1,2; xcd(i,i)=xcd(i,i)+zp; enddo
end subroutine xstoxc1

!> Inverse of xttoxs.
!!
!! @param[in] k Gaussian curvature parameter of Schmidt mapping 
!! @param[in] xs Stereographic plane coordinates
!! @param[out] xt Scaled gnomonic plane coordinates 
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xstoxt(k,xs,xt,ff)!                                        [xstoxt]
implicit none
real(dp),             intent(in ):: k
real(dp),dimension(2),intent(in ):: xs
real(dp),dimension(2),intent(out):: xt
logical,              intent(out):: ff
real(dp):: s,sc
s=k*(xs(1)*xs(1)+xs(2)*xs(2)); sc=u1-s
ff=abs(s)>=u1; if(ff)return
xt=u2*xs/sc
end subroutine xstoxt

!> Scaled gnomonic plane xt to standard stereographic plane xs.
!!
!! @param[in] k Gaussian curvature parameter of Schmidt mapping 
!! @param[in] xt Scaled gnomonic plane
!! @param[out] xs Standard stereographic plane
!! @param[out] xsd Jacobian matrix, d(xs)/d(xt). 
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xttoxs(k,xt,xs,xsd,ff)!                                     [xttoxs
use pmat4, only: outer_product
implicit none
real(dp),               intent(in ):: k
real(dp),dimension(2),  intent(in ):: xt
real(dp),dimension(2),  intent(out):: xs
real(dp),dimension(2,2),intent(out):: xsd
logical,                intent(out):: ff
real(dp),dimension(2):: rspd
real(dp)             :: s,sp,rsp,rsppi,rsppis
integer(spi)         :: i
s=k*dot_product(xt,xt); sp=u1+s
ff=(sp<=u0); if(ff)return
rsp=sqrt(sp)
rsppi=u1/(u1+rsp)
rsppis=rsppi**2
xs=xt*rsppi
rspd=k*xt/rsp
xsd=-outer_product(xt,rspd)*rsppis
do i=1,2; xsd(i,i)=xsd(i,i)+rsppi; enddo
end subroutine xttoxs

!> Like xttoxs, but also, return the derivatives, wrt K, of xs and
!! xsd.
!!
!! @param[in] k Gaussian curvature parameter of the Schmidt mapping
!! @param[in] xt Scaled gnomonic plane
!! @param[out] xs Standard stereographic plane
!! @param[out] xsd Jacobian matrix, d(xs)/d(xt)
!! @param[out] xsdd Second partial derivatives, d^2(xs)/(d(xt)d(xt))
!! @param[out] xs1 Derivative of xs wrt mapping parameter, d(xs)/dk
!! @param[out] xsd1 Derivative of Jacobian wrt k: d^2(xs)/(d(xt)dk)
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xttoxs1(k,xt,xs,xsd,xsdd,xs1,xsd1,ff)!                     [xttoxs]
use pmat4, only: outer_product
implicit none
real(dp),                 intent(in ):: k
real(dp),dimension(2),    intent(in ):: xt
real(dp),dimension(2),    intent(out):: xs ,xs1
real(dp),dimension(2,2),  intent(out):: xsd,xsd1
real(dp),dimension(2,2,2),intent(out):: xsdd
logical,                  intent(out):: ff
real(dp),dimension(2,2):: rspdd
real(dp),dimension(2)  :: rspd,rspd1,rsppid
real(dp)               :: s,sp,rsp,rsppi,rsppis,s1,rsp1
integer(spi)           :: i
s1=dot_product(xt,xt); s=k*s1; sp=u1+s
ff=(sp<=u0); if(ff)return
rsp=sqrt(sp);      rsp1=o2*s1/rsp
rsppi=u1/(u1+rsp); rsppis=rsppi**2
xs=xt*rsppi;       xs1=-xt*rsp1*rsppis
rspd=k*xt/rsp;     rspd1=(xt*rsp-k*xt*rsp1)/sp
rsppid=-rspd*rsppis
xsd1=-outer_product(xt,rspd1-u2*rspd*rsp1*rsppi)
do i=1,2; xsd1(i,i)=xsd1(i,i)-rsp1; enddo; xsd1=xsd1*rsppis

xsd=-outer_product(xt,rspd)*rsppis
do i=1,2; xsd(i,i)=xsd(i,i)+rsppi; enddo

rspdd=-outer_product(xt,rspd)*rsppi
xsdd=u0
do i=1,2; xsdd(i,:,i)=            rsppid;                 enddo
do i=1,2; xsdd(i,i,:)=xsdd(i,i,:)+rsppid;                 enddo
do i=1,2; xsdd(:,:,i)=xsdd(:,:,i)+u2*rspdd*rsppid(i);     enddo
do i=1,2; rspdd(i,i)=rspdd(i,i)+rsp*rsppi;                enddo
do i=1,2; xsdd(i,:,:)=xsdd(i,:,:)-xt(i)*rspdd*rsppi*k/sp; enddo
end subroutine xttoxs1

!> Inverse of xmtoxt.
!!
!! @param[in] a Mapping parameter controlling grid line spacing profile  
!! @param[in] xt Gnomonic plane coordinates  
!! @param[out] xm Map coordinates
!! @param[out] ff Failure flag 
!! @author R. J. Purser
subroutine xttoxm(a,xt,xm,ff)!                                       [xttoxm]
implicit none
real(dp),             intent(in ):: a
real(dp),dimension(2),intent(in ):: xt
real(dp),dimension(2),intent(out):: xm
logical              ,intent(out):: ff
integer(spi):: i
do i=1,2; call zttozm(a,xt(i),xm(i),ff); if(ff)return; enddo
end subroutine xttoxm

!> Like zmtozt, but for 2-vector xm and xt, and 2*2 diagonal Jacobian
!! xtd.
!!
!! @param[in] a Mapping parameter controlling grid line spacing profile 
!! @param[in] xm Vector value of map coordinates
!! @param[out] xt Vector value of gnomonic plane coordinates
!! @param[out] xtd 2*2 diagonal Jacobian, d(xt)/d(xm)
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xmtoxt(a,xm,xt,xtd,ff)!                                    [xmtoxt]
implicit none
real(dp),               intent(in ):: a
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(2),  intent(out):: xt
real(dp),dimension(2,2),intent(out):: xtd
logical,                intent(out):: ff
integer(spi):: i
xtd=u0; do i=1,2; call zmtozt(a,xm(i),xt(i),xtd(i,i),ff); if(ff)return; enddo
end subroutine xmtoxt

!> Like zmtozt1, but for 2-vector xm and xt, and 2*2 diagonal Jacobian
!! xtd Also, the derivatives, wrt a, of these quantities.
!!
!! @param[in] a Mapping parameter controlling grid line spacing profile
!! @param[in] xm Vector value of map plane coordinates
!! @param[out] xt Vector value of gnomonic plane coordinates
!! @param[out] xtd 2*2 diagonal Jacobian, d(xt)/d(xm)
!! @param[out] xt1 Derivative wrt a of xt, d(xt)/da 
!! @param[out] xtd1 Derivative wrt a of Jacobian xtd, d^2(xt)/(d(xm)da)
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xmtoxt1(a,xm,xt,xtd,xt1,xtd1,ff)!                          [xmtoxt]
implicit none
real(dp),                 intent(in ):: a
real(dp),dimension(2),    intent(in ):: xm
real(dp),dimension(2),    intent(out):: xt,xt1
real(dp),dimension(2,2),  intent(out):: xtd,xtd1
logical,                  intent(out):: ff
integer(spi):: i
xtd=u0
xtd1=u0
do i=1,2
   call zmtozt1(a,xm(i),xt(i),xtd(i,i),xt1(i),xtd1(i,i),ff)
   if(ff)return
enddo
end subroutine xmtoxt1

!> Inverse of zmtozt
!!
!! @param[in] a Mapping parameter controlling grid line spacing profile 
!! @param[in] zt Scalar value of single gnomonic plane coordinate
!! @param[out] zm Scalar value of single map plane coordinate
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine zttozm(a,zt,zm,ff)!                                        [zttozm]
implicit none
real(dp),intent(in ):: a,zt
real(dp),intent(out):: zm
logical, intent(out):: ff
real(dp):: ra,razt
ff=F
if    (a>u0)then; ra=sqrt( a); razt=ra*zt; zm=atan (razt)/ra
elseif(a<u0)then; ra=sqrt(-a); razt=ra*zt; ff=abs(razt)>=u1; if(ff)return
                                           zm=atanh(razt)/ra
else                                     ; zm=zt
endif
end subroutine zttozm

!> Evaluate the function, zt = tan(sqrt(A)*z)/sqrt(A), and its
!! derivative, ztd, for positive and negative A and for the limiting
!! case, A --> 0.
!!
!! @param[in] a Mapping parameter controlling grid line spacing profile 
!! @param[in] zm Scalar value of single map plane coordinate 
!! @param[out] zt Scalar value of single gnomonic plane coordinate
!! @param[out] ztd Derivative of gnomonic coordinate, d(zt)/d(zm)
!! @param[out] ff Failure flag 
!! @author R. J. Purser
subroutine zmtozt(a,zm,zt,ztd,ff)!                                    [zmtozt]
implicit none
real(dp),intent(in ):: a,zm
real(dp),intent(out):: zt,ztd
logical, intent(out):: ff
real(dp):: ra
ff=f
if    (a>u0)then; ra=sqrt( a); zt=tan (ra*zm)/ra; ff=abs(ra*zm)>=pih
elseif(a<u0)then; ra=sqrt(-a); zt=tanh(ra*zm)/ra
else                         ; zt=zm
endif
ztd=u1+a*zt*zt
end subroutine zmtozt

!> Like zmtozt, but also, get the derivative with respect to a, zt1 of
!! zt, and ztd1 of ztd.
!! 
!! @param[in] a Mapping parameter controlling grid line spacing profile 
!! @param[in] zm Single map plane coordinate
!! @param[in] zt Single gnomonic plane coordinate
!! @param[in] ztd Derivative wrt zm of zt, d(zt)/d(zm)
!! @param[in] zt1 Derivative wrt a of zt, d(zt)/da 
!! @param[in] ztd1 Derivative wrt a of ztd, d^2(zt)/(d(zm)da) 
!! @param[in] ff Failure flag 
!! @author R. J. Purser
subroutine zmtozt1(a,zm,zt,ztd,zt1,ztd1,ff)!                          [zmtozt]
use pietc, only: o3
use pfun,  only: sinoxm,sinox,sinhoxm,sinhox
implicit none
real(dp),intent(in ):: a,zm
real(dp),intent(out):: zt,ztd,zt1,ztd1
logical, intent(out):: ff
real(dp):: ra,rad,razm
ff=f
if    (a>u0)then;ra=sqrt( a);razm=ra*zm; zt=tan(razm)/ra; ff=abs(razm)>=pih
rad=o2/ra
zt1=(rad*zm/ra)*((-u2*sin(razm*o2)**2-sinoxm(razm))/cos(razm)+(tan(razm))**2)
elseif(a<u0)then;ra=sqrt(-a);razm=ra*zm; zt=tanh(razm)/ra
rad=-o2/ra
zt1=(rad*zm/ra)*((u2*sinh(razm*o2)**2-sinhoxm(razm))/cosh(razm)-(tanh(razm))**2)
else                        ;zt=zm; zt1=zm**3*o3
endif
ztd=u1+a*zt*zt
ztd1=zt*zt +u2*a*zt*zt1
end subroutine zmtozt1

!> Assuming the vector AK parameterization of the Extended Schmidt-transformed
!! Gnomonic (ESG) mapping with parameter vector, and given a map-space
!! 2-vector, xm, find the corresponding cartesian unit 3-vector and its
!! derivative wrt xm, the Jacobian matrix, xcd.
!! If for any reason the mapping cannot be done, return a raised failure flag,z
!! FF.
!! @param [in] ak 2-vector parameterization of the ESG mapping
!! @param [in] xm 2-vector of map plane coordinates
!! @param [out] xc Earth-centered cartesian unit 3-vector
!! @param [out] xcd Jacobian, d(xc)/d(xm)
!! @param [out] ff Failure flag
!! @author R. J. Purser
subroutine xmtoxc_vak(ak,xm,xc,xcd,ff)!                            [xmtoxc_ak]
implicit none
real(dp),dimension(2),  intent(in ):: ak,xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
logical,                intent(out):: ff
call xmtoxc_ak(ak(1),ak(2),xm,xc,xcd,ff)
end subroutine xmtoxc_vak

!> Like xmtoxc_vak, _ak, but also return derivatives wrt ak.
!!
!! @param[in] ak 2-vector parameterization of the ESG mapping 
!! @param[in] xm 2-vector of map plane coordinates
!! @param[out] xc Earth-centered cartesian unit 3-vector
!! @param[out] xcd Jacobian of xc wrt xm, d(xc)/d(xm)
!! @param[out] xc1 Partial derivatives wrt ak of xc, d(xc)/d(ak)
!! @param[out] xcd1 Second derivative wrt xm and ak of xc, d^2(xc)/(d(xm)d(ak))
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xmtoxc_vak1(ak,xm,xc,xcd,xc1,xcd1,ff)!                  [xmtoxc_ak]
implicit none
real(dp),dimension(2),    intent(in ):: ak,xm
real(dp),dimension(3),    intent(out):: xc
real(dp),dimension(3,2),  intent(out):: xcd
real(dp),dimension(3,2),  intent(out):: xc1
real(dp),dimension(3,2,2),intent(out):: xcd1
logical,                  intent(out):: ff
real(dp),dimension(3,2,2):: xcdd
real(dp),dimension(2,2,2):: xsd1,xsdd
real(dp),dimension(2,2)  :: xtd,xsd,xs1,xtd1,xsdk
real(dp),dimension(2)    :: xt,xt1,xs,xsk
integer(spi)             :: i
call xmtoxt1(ak(1),xm,xt,xtd,xt1,xtd1,ff);      if(ff)return
call xttoxs1(ak(2),xt,xs,xsd,xsdd,xsk,xsdk,ff); if(ff)return
xs1(:,2)=xsk; xs1(:,1)=matmul(xsd,xt1)
xsd1(:,:,1)=matmul(xsd,xtd1)+matmul(xsdd(:,:,1)*xt1(1)+xsdd(:,:,2)*xt1(2),xtd)
xsd1(:,:,2)=matmul(xsdk,xtd)
xsd=matmul(xsd,xtd)
call xstoxc(xs,xc,xcd,xcdd)
xc1=matmul(xcd,xs1)
do i=1,3; xcd1(i,:,:)=matmul(transpose(xsd),matmul(xcdd(i,:,:),xs1)); enddo
do i=1,2; xcd1(:,:,i)=xcd1(:,:,i)+matmul(xcd,xsd1(:,:,i));            enddo
xcd=matmul(xcd,xsd)
end subroutine xmtoxc_vak1

!> Assuming the A-K parameterization of the Extended
!! Schmidt-transformed Gnomonic (ESG) mapping, and given a map-space
!! 2-vector, xm, find the corresponding cartesian unit 3-vector and
!! its derivative wrt xm, jacobian matrix, xcd. If for any reason the
!! mapping cannot be done, return a raised failure flag, FF.
!!
!! @param[in] a ESG mapping parameter for line spacing profile
!! @param[in] k ESG mapping parameter for Gauss curvature in Schmidt mapping
!! @param[in] xm map-space 2-vector
!! @param[out] xc Earth-centered cartesian unit 3-vector
!! @param[out] xcd Jacobian matrix, d(xc)/d(xm)
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xmtoxc_ak(a,k,xm,xc,xcd,ff)!                            [xmtoxc_ak]
implicit none
real(dp),               intent(in ):: a,k
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
logical,                intent(out):: ff
real(dp),dimension(2,2):: xtd,xsd
real(dp),dimension(2)  :: xt,xs
call xmtoxt(a,xm,xt,xtd,ff); if(ff)return
call xttoxs(k,xt,xs,xsd,ff); if(ff)return
xsd=matmul(xsd,xtd)
call xstoxc(xs,xc,xcd)
xcd=matmul(xcd,xsd)
end subroutine xmtoxc_ak

!> Inverse mapping of xmtoxc_ak. That is, go from given cartesian unit
!! 3-vector, xc, to map coordinate 2-vector xm (or return a raised
!! failure flag, FF, if the attempt fails).
!!
!! @param[in] a ESG mapping parameter for line spacing profile
!! @param[in] k ESG mapping parameter for Gauss curvature in Schmidt mapping
!! @param[in] xc Earth-centered cartesian unit 3-vector
!! @param[out] xm 2-vector map coordinate
!! @param[out] ff Failure flag
!! @author R. J. Purser
subroutine xctoxm_ak(a,k,xc,xm,ff)!                                [xctoxm_ak]
implicit none
real(dp),             intent(in ):: a,k
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
real(dp),dimension(2):: xs,xt
ff=F
call xctoxs(xc,xs)
call xstoxt(k,xs,xt,ff); if(ff)return
call xttoxm(a,xt,xm,ff)
end subroutine xctoxm_ak

!> For angles (degrees) of the arcs spanning the halfwidths between
!! the region's center and its x and y edges, get the two cartesian
!! vectors that represent the locations of these edge midpoints in the
!! positive x and y directions.
!!
!! @param[in] arcx Center-relative angle (degrees) of edge midpoint in +x
!! @param[in] arcy Center-relative angle (degrees) of edge midpoint in +y
!! @param[out] edgex region's +x edge midpoint as cartesian unit 3-vector
!! @param[out] edgey region's +y edge midpoint as cartesian unit 3-vector
!! @author R. J. Purser
subroutine get_edges(arcx,arcy,edgex,edgey)!                       [get_edges]
implicit none
real(dp),             intent(in ):: arcx,arcy
real(dp),dimension(3),intent(out):: edgex,edgey
real(dp):: cx,sx,cy,sy,rarcx,rarcy
rarcx=arcx*dtor;   rarcy=arcy*dtor
cx=cos(rarcx); sx=sin(rarcx)
cy=cos(rarcy); sy=sin(rarcy)
edgex=(/sx,u0,cx/); edgey=(/u0,sy,cy/)
end subroutine get_edges

!> From a jacobian matrix, j0, get a sufficient set of v.. diagnostics
!! such that, from averages of these v, we can later compute the
!! collective variance of Q(lam) that they imply for any choice of the
!! "lambda" parameter, lam.  Note that v1 and v4 are quadratic
!! diagnostics of EL, while v2 and v3 are linear.
!!
!! @param[in] j0 jacobian matrix
!! @param[out] v1 quadratic diagnostics of EL
!! @param[out] v2 linear diagnostics of EL
!! @param[out] v3 linear diagnostics of EL
!! @param[out] v4 quadratic diagnostics of EL
!! @author R. J. Purser
subroutine get_qx(j0, v1,v2,v3,v4)!                                   [get_qx]
use psym2, only: logsym2
implicit none
real(dp),dimension(3,2),intent(in ):: j0
real(dp),               intent(out):: v1,v2,v3,v4
real(dp),dimension(2,2):: el,g
g=matmul(transpose(j0),j0)
call logsym2(g,el); el=el*o2
v1=el(1,1)**2+u2*el(1,2)**2+el(2,2)**2
v2=el(1,1)
v3=el(2,2)
v4=(el(1,1)+el(2,2))**2
end subroutine get_qx

!> From a jacobian matrix, j0, and its derivative, j0d, get a
!! sufficient set of v.. diagnostics such that, from average of these
!! diagnostics, we can later compute the collective variance of Q and
!! its derivative.
!!
!! @param[in] j0 jacobian matrix
!! @param[in] j0d derivative of j0
!! @param[in] v1
!! @param[in] v2
!! @param[in] v3
!! @param[in] v4
!! @param[in] v1d
!! @param[in] v2d
!! @param[in] v3d
!! @param[in] v4d
!! @author R. J. Purser
subroutine get_qxd(j0,j0d, v1,v2,v3,v4,v1d,v2d,v3d,v4d)!              [get_qx]
use psym2, only: logsym2
implicit none
real(dp),dimension(3,2),  intent(in ):: j0
real(dp),dimension(3,2,2),intent(in ):: j0d
real(dp),                 intent(out):: v1,v2,v3,v4
real(dp),dimension(2),    intent(out):: v1d,v2d,v3d,v4d
real(dp),dimension(2,2,2,2):: deldg
real(dp),dimension(2,2,2)  :: eld,gd
real(dp),dimension(2,2)    :: el,g
integer(spi)               :: i,j,k
g=matmul(transpose(j0),j0)
do i=1,2
   gd(:,:,i)=matmul(transpose(j0d(:,:,i)),j0)+matmul(transpose(j0),j0d(:,:,i))
enddo
call logsym2(g,el,deldg); el=el*o2; deldg=deldg*o2
eld=u0
do i=1,2; do j=1,2; do k=1,2
   eld(:,:,k)=eld(:,:,k)+deldg(:,:,i,j)*gd(i,j,k)
enddo   ; enddo   ; enddo
v1=el(1,1)**2+u2*el(1,2)**2+el(2,2)**2
v2=el(1,1)
v3=el(2,2)
v4=(el(1,1)+el(2,2))**2
v1d=u2*(el(1,1)*eld(1,1,:)+u2*el(1,2)*eld(1,2,:)+el(2,2)*eld(2,2,:))
v2d=eld(1,1,:)
v3d=eld(2,2,:)
v4d=u2*(el(1,1)+el(2,2))*(eld(1,1,:)+eld(2,2,:))
end subroutine get_qxd

!> For a parameter vector, ak and a map-space domain-parameter vector,
!! ma, return the lambda-parameterized quality diagnostic, Q, and the
!! geographic domain-parameter vector ga. Lambda is given by lam
!! <1. Also, return the derivatives, qdak and qdma, of Q wrt ak and
!! ma, and the derivatives gadak and gadma, of ga wrt ak and ma.
!!
!! The domain averages of Q are accurately computed by
!! bi-Gauss-Legendre quadrature over the positive quadrant of the
!! domain (exploiting the symmetry) of the four constituent terms, v1,
!! v2, v3, v4, from which the mean Q is computed using a quadratic
!! formula of these constituents.  The number of Gauss points in eaxh
!! half-interval is ngh, and the nodes themselves are, in proportion
!! to the half-interval, at xg.  the normalized gauss weights are wg.
!!
!! If a failure occurs, colmputations cease immediately and a failure
!! flag, FF, is raised on return.
!!
!! @param[in] ngh
!! @param[in] lam Lambda
!! @param[in] xg
!! @param[in] wg
!! @param[in] ak parameter vector
!! @param[in] ma map-space domain-parameter vector
!! @param[out] q lambda-parameterized quality diagnostic
!! @param[out] qdak derivatives value
!! @param[out] qdma derivatives value
!! @param[out] ga geographic domain-parameter vector
!! @param[out] gadak
!! @param[out] gadma
!! @param[out] ff error flag
!! @author R. J. Purser
subroutine get_meanqd(ngh,lam,xg,wg,ak,ma, q,qdak,qdma, & !        [get_meanq]
     ga,gadak,gadma, ff)
implicit none
integer(spi),           intent(in ):: ngh
real(dp),               intent(in ):: lam
real(dp),dimension(ngh),intent(in ):: xg,wg
real(dp),dimension(2)  ,intent(in ):: ak,ma
real(dp),               intent(out):: q
real(dp),dimension(2),  intent(out):: qdak,qdma
real(dp),dimension(2),  intent(out):: ga
real(dp),dimension(2,2),intent(out):: gadak,gadma
logical,                intent(out):: ff
real(dp),dimension(3,2,2):: xcd1
real(dp),dimension(3,2)  :: xcd,xc1
real(dp),dimension(3)    :: xc
real(dp),dimension(2)    :: xm, v1dxy,v2dxy,v3dxy,v4dxy,                &
                            v1dL,v2dL,v3dL,v4dL, v1d,v2d,v3d,v4d
real(dp)                 :: wx,wy,                                  &
                            v1xy,v2xy,v3xy,v4xy, v1L,v2L,v3L,v4L, v1,v2,v3,v4
integer(spi)             :: i,ic,ix,iy
v1 =u0; v2 =u0; v3 =u0; v4 =u0
v1d=u0; v2d=u0; v3d=u0; v4d=u0
do iy=1,ngh
   wy=wg(iy)
   xm(2)=ma(2)*xg(iy)
   v1L =u0; v2L =u0; v3L =u0; v4L =u0
   v1dL=u0; v2dL=u0; v3dL=u0; v4dL=u0
   do ix=1,ngh
      wx=wg(ix)
      xm(1)=ma(1)*xg(ix)
      call xmtoxc_ak(ak,xm,xc,xcd,xc1,xcd1,ff); if(ff)return
      call get_qx(xcd,xcd1, v1xy,v2xy,v3xy,v4xy, v1dxy,v2dxy,v3dxy,v4dxy)
      v1L =v1L +wx*v1xy; v2L =v2L +wx*v2xy 
      v3L =v3L +wx*v3xy; v4L =v4L +wx*v4xy
      v1dL=v1dL+wx*v1dxy;v2dL=v2dL+wx*v2dxy 
      v3dL=v3dL+wx*v3dxy;v4dL=v4dL+wx*v4dxy
   enddo
   v1 =v1 +wy*v1L;  v2 =v2 +wy*v2L;  v3 =v3 +wy*v3L;  v4 =v4 +wy*v4L
   v1d=v1d+wy*v1dL; v2d=v2d+wy*v2dL; v3d=v3d+wy*v3dL; v4d=v4d+wy*v4dL
enddo
call get_qofv(lam,v1,v2,v3,v4, q)! <- Q(lam) based on the v1,v2,v3,v4 
call get_qofv(lam,v2,v3, v1d,v2d,v3d,v4d, qdak)! <- Derivative of Q wrt ak
! Derivatives of ga wrt ak, and of q and ga wrt ma:
gadma=u0! <- needed because only diagonal elements are filled
do i=1,2
   ic=3-i
   xm=0; xm(i)=ma(i)
   call xmtoxc_ak(ak,xm,xc,xcd,xc1,xcd1,ff); if(ff)return
   ga(i)=atan2(xc(i),xc(3))*rtod
   gadak(i,:)=(xc(3)*xc1(i,:)-xc(i)*xc1(3,:))*rtod
   gadma(i,i)=(xc(3)*xcd(i,i)-xc(i)*xcd(3,i))*rtod

   v1L=u0; v2L=u0; v3L=u0; v4L=u0
   do iy=1,ngh
      wy=wg(iy)
      xm(ic)=ma(ic)*xg(iy)
      call xmtoxc_ak(ak,xm,xc,xcd,ff); if(ff)return
      call get_qx(xcd, v1xy,v2xy,v3xy,v4xy)
      v1L=v1L+wy*v1xy; v2L=v2L+wy*v2xy; v3L=v3L+wy*v3xy; v4L=v4L+wy*v4xy
   enddo
   v1d(i)=(v1L-v1)/ma(i); v2d(i)=(v2L-v2)/ma(i)
   v3d(i)=(v3L-v3)/ma(i); v4d(i)=(v4L-v4)/ma(i)
enddo
call get_qofv(lam,v2,v3, v1d,v2d,v3d,v4d, qdma)
end subroutine get_meanqd

!> Like getmeanqd, except for n different values, aks, of ak and n
!! different values, mas of ma, and without any of the derivatives.
!!
!! @param[in] n
!! @param[in] ngh
!! @param[in] lam
!! @param[in] xg
!! @param[in] wg
!! @param[in] aks
!! @param[in] mas
!! @param[out] qs
!! @param[out] ff
!! @author R. J. Purser
subroutine get_meanqs(n,ngh,lam,xg,wg,aks,mas, qs,ff)!            [get_meanq]
implicit none
integer(spi),           intent(in ):: n,ngh
real(dp),dimension(ngh),intent(in ):: xg,wg
real(dp),               intent(in ):: lam
real(dp),dimension(2,n),intent(in ):: aks,mas
real(dp),dimension(n),  intent(out):: qs
logical,                intent(out):: ff
real(dp),dimension(n)  :: v1s,v2s,v3s,v4s
real(dp),dimension(n)  :: v1sL,v2sL,v3sL,v4sL
real(dp),dimension(3,2):: xcd
real(dp),dimension(3)  :: xc
real(dp),dimension(2)  :: xm,xgs
real(dp)               :: wx,wy, v1xy,v2xy,v3xy,v4xy
integer(spi)           :: i,ix,iy
v1s=u0; v2s=u0; v3s=u0; v4s=u0
do iy=1,ngh
   wy=wg(iy)
   v1sL=u0; v2sL=u0; v3sL=u0; v4sL=u0
   do ix=1,ngh
      wx=wg(ix)
      xgs=(/xg(ix),xg(iy)/)
      do i=1,n
         xm=mas(:,i)*xgs
         call xmtoxc_ak(aks(:,i),xm,xc,xcd,ff); if(ff)return
         call get_qx(xcd,v1xy,v2xy,v3xy,v4xy)
         v1sL(i)=v1sL(i)+wx*v1xy; v2sL(i)=v2sL(i)+wx*v2xy
         v3sL(i)=v3sL(i)+wx*v3xy; v4sL(i)=v4sL(i)+wx*v4xy
      enddo
   enddo
   v1s=v1s+wy*v1sL; v2s=v2s+wy*v2sL; v3s=v3s+wy*v3sL; v4s=v4s+wy*v4sL
enddo
call get_qofv(n,lam,v1s,v2s,v3s,v4s, qs)
end subroutine get_meanqs

!> The quadratic quantity Q depends linearly on v1 and v4 (which are
!! already quadratic diagnostics of EL) and quadratically on v2 and v3
!! (which are linear diagnostics of EL). EL = (1/2)log(G), where
!! G=J^T.J, J the jacobian.
!!
!! @param[in] lam
!! @param[in] v1 quadratic diagnostics of EL
!! @param[in] v2 linear diagnostics of EL
!! @param[in] v3 linear diagnostics of EL
!! @param[in] v4 quadratic diagnostics of EL
!! @param[out] q quadratic quantity
!! @author R. J. Purser
subroutine get_qofv(lam,v1,v2,v3,v4, q)!                            [get_qofv]
implicit none
real(dp),intent(in ):: lam,v1,v2,v3,v4
real(dp),intent(out):: q
real(dp):: lamc
lamc=u1-lam
q=lamc*(v1-(v2**2+v3**2)) +lam*(v4 -(v2+v3)**2)
end subroutine get_qofv

!> Like get_qofv, but for (only) the 2-vector derivatives of Q. Note
!! that the quadratic diagnostics v1 and v4 do not participate in this
!! formula.
!!
!! @param[in] lam
!! @param[in] v2
!! @param[in] v3
!! @param[in] v1d
!! @param[in] v2d
!! @param[in] v3d
!! @param[in] v4d
!! @param[out] qd
!! @author R. J. Purser
subroutine get_qofvd(lam, v2,v3, v1d,v2d,v3d,v4d, qd)!              [get_qofv]
implicit none
real(dp),             intent(in ):: lam,v2,v3
real(dp),dimension(2),intent(in ):: v1d,v2d,v3d,v4d
real(dp),dimension(2),intent(out):: qd
real(dp):: lamc
lamc=u1-lam
qd=lamc*(v1d-u2*(v2d*v2+v3d*v3))+lam*(v4d-u2*(v2d+v3d)*(v2+v3))
end subroutine get_qofvd

!> General util to convert value.
!!
!! @param[in] n
!! @param[in] lam
!! @param[in] v1s
!! @param[in] v2s
!! @param[in] v3s
!! @param[in] v4s
!! @param[out] qs
!! @author R. J. Purser
subroutine get_qsofvs(n,lam,v1s,v2s,v3s,v4s, qs)!                   [get_qofv]
implicit none
integer(spi),         intent(in ):: n
real(dp),             intent(in ):: lam
real(dp),dimension(n),intent(in ):: v1s,v2s,v3s,v4s
real(dp),dimension(n),intent(out):: qs
real(dp):: lamc
lamc=u1-lam
qs=lamc*(v1s-(v2s**2+v3s**2)) +lam*(v4s -(v2s+v3s)**2)
end subroutine get_qsofvs

!> Given an aspect ratio, asp<=1, and major semi-axis, arc, in
!! map-space nondimensional units, return a first guess for the
!! parameter vector, ak, approximately optimal for the domain of the
!! given dimensions.
!!
!! @param[in] asp aspect ratio
!! @param[in] tmarcx
!! @param[out] ak first guess for the parameter vector
!! @author R. J. Purser
subroutine guessak_map(asp,tmarcx,ak)!                           [guessak_map]
implicit none
real(dp),             intent(in ):: asp,tmarcx
real(dp),dimension(2),intent(out):: ak
real(dp):: gmarcx
gmarcx=tmarcx*rtod
call guessak_geo(asp,gmarcx,ak)
end subroutine guessak_map

!> Given an aspect ratio, asp<=1, and major semi-axis, arc, in
!! geographical (degree) units measured along the rectangle's median,
!! return a first guess for the parameter vector, ak, approximately
!! optimal for the domain of the given dimensions.
!!
!! @param asp aspect ratio of intended domain
!! @param arc major semi-axis angle in degrees for intended domain
!! @param ak first guess of the parameter vector
!! @author R. J. Purser
subroutine guessak_geo(asp,arc,ak)!                              [guessak_geo]
implicit none
real(dp),             intent(in ):: asp,arc
real(dp),dimension(2),intent(out):: ak
integer(spi),parameter         :: narc=11,nasp=10! <- Table index bounds
real(dp),    parameter         :: eps=1.e-7_dp,darc=10._dp+eps,dasp=.1_dp+eps
real(dp),dimension(nasp,0:narc):: adarc,kdarc
real(dp)                       :: sasp,sarc,wx0,wx1,wa0,wa1
integer(spi)                   :: iasp0,iasp1,iarc0,iarc1
!------------------------
! Tables of approximate A (adarc) and K (kdarc), valid for lam=0.8, for aspect ratio,
! asp, at .1, .2, .3, .4, .5, .6, .7, .8, 1. and major semi-arc, arcx, at values,
! 0 (nominally), 10., ..., 100. degrees (where the nominal "0" angle is actually
! from a computation at 3 degrees, since zero would not make sense).
! The 100 degree rows are repeated as the 110 degree entries deliberately to pad the
! table into the partly forbidden parameter space to allow skinny domains of up to
! 110 degree semi-length to be validly endowed with optimal A and K:
data adarc/ &
-.450_dp,-.328_dp,-.185_dp,-.059_dp,0.038_dp,0.107_dp,0.153_dp,0.180_dp,0.195_dp,0.199_dp,&
-.452_dp,-.327_dp,-.184_dp,-.058_dp,0.039_dp,0.108_dp,0.154_dp,0.182_dp,0.196_dp,0.200_dp,&
-.457_dp,-.327_dp,-.180_dp,-.054_dp,0.043_dp,0.112_dp,0.158_dp,0.186_dp,0.200_dp,0.205_dp,&
-.464_dp,-.323_dp,-.173_dp,-.047_dp,0.050_dp,0.118_dp,0.164_dp,0.192_dp,0.208_dp,0.213_dp,&
-.465_dp,-.313_dp,-.160_dp,-.035_dp,0.060_dp,0.127_dp,0.173_dp,0.202_dp,0.217_dp,0.224_dp,&
-.448_dp,-.288_dp,-.138_dp,-.017_dp,0.074_dp,0.140_dp,0.184_dp,0.213_dp,0.230_dp,0.237_dp,&
-.395_dp,-.244_dp,-.104_dp,0.008_dp,0.093_dp,0.156_dp,0.199_dp,0.227_dp,0.244_dp,0.252_dp,&
-.301_dp,-.177_dp,-.057_dp,0.042_dp,0.119_dp,0.175_dp,0.215_dp,0.242_dp,0.259_dp,0.269_dp,&
-.185_dp,-.094_dp,0.001_dp,0.084_dp,0.150_dp,0.199_dp,0.235_dp,0.260_dp,0.277_dp,0.287_dp,&
-.069_dp,-.006_dp,0.066_dp,0.132_dp,0.186_dp,0.227_dp,0.257_dp,0.280_dp,0.296_dp,0.308_dp,&
0.038_dp,0.081_dp,0.134_dp,0.185_dp,0.226_dp,0.258_dp,0.283_dp,0.303_dp,0.319_dp,0.333_dp,&
0.038_dp,0.081_dp,0.134_dp,0.185_dp,0.226_dp,0.258_dp,0.283_dp,0.303_dp,0.319_dp,0.333_dp/

data kdarc/ &
-.947_dp,-.818_dp,-.668_dp,-.535_dp,-.433_dp,-.361_dp,-.313_dp,-.284_dp,-.269_dp,-.264_dp,&
-.946_dp,-.816_dp,-.665_dp,-.533_dp,-.431_dp,-.359_dp,-.311_dp,-.282_dp,-.267_dp,-.262_dp,&
-.942_dp,-.806_dp,-.655_dp,-.524_dp,-.424_dp,-.353_dp,-.305_dp,-.276_dp,-.261_dp,-.255_dp,&
-.932_dp,-.789_dp,-.637_dp,-.509_dp,-.412_dp,-.343_dp,-.296_dp,-.266_dp,-.250_dp,-.244_dp,&
-.909_dp,-.759_dp,-.609_dp,-.486_dp,-.394_dp,-.328_dp,-.283_dp,-.254_dp,-.237_dp,-.230_dp,&
-.863_dp,-.711_dp,-.569_dp,-.456_dp,-.372_dp,-.310_dp,-.266_dp,-.238_dp,-.220_dp,-.212_dp,&
-.779_dp,-.642_dp,-.518_dp,-.419_dp,-.343_dp,-.287_dp,-.247_dp,-.220_dp,-.202_dp,-.192_dp,&
-.661_dp,-.556_dp,-.456_dp,-.374_dp,-.310_dp,-.262_dp,-.226_dp,-.200_dp,-.182_dp,-.171_dp,&
-.533_dp,-.462_dp,-.388_dp,-.325_dp,-.274_dp,-.234_dp,-.203_dp,-.179_dp,-.161_dp,-.150_dp,&
-.418_dp,-.373_dp,-.322_dp,-.275_dp,-.236_dp,-.204_dp,-.178_dp,-.156_dp,-.139_dp,-.127_dp,&
-.324_dp,-.296_dp,-.261_dp,-.229_dp,-.200_dp,-.174_dp,-.152_dp,-.133_dp,-.117_dp,-.104_dp,&
-.324_dp,-.296_dp,-.261_dp,-.229_dp,-.200_dp,-.174_dp,-.152_dp,-.133_dp,-.117_dp,-.104_dp/
!=============================================================================
sasp=asp/dasp
iasp0=floor(sasp); wa1=sasp-iasp0
iasp1=iasp0+1;     wa0=iasp1-sasp
sarc=arc/darc
iarc0=floor(sarc); wx1=sarc-iarc0
iarc1=iarc0+1;      wx0=iarc1-sarc
if(iasp0<1 .or. iasp1>nasp)stop 'Guessak_geo; Aspect ratio out of range'
if(iarc0<0 .or. iarc1>narc)stop 'Guessak_geo; Major semi-arc is out of range'

! Bilinearly interpolate A and K from crude table into a 2-vector:
ak=(/wx0*(wa0*adarc(iasp0,iarc0)+wa1*adarc(iasp1,iarc0))+ &
     wx1*(wa0*adarc(iasp0,iarc1)+wa1*adarc(iasp1,iarc1)), &
     wx0*(wa0*kdarc(iasp0,iarc0)+wa1*kdarc(iasp1,iarc0))+ &
     wx1*(wa0*kdarc(iasp0,iarc1)+wa1*kdarc(iasp1,iarc1))/)
end subroutine guessak_geo

!> Get the best Extended Schmidt Gnomonic parameter, (a,k), for the given
!! geographical half-spans, garcx and garcy, as well as the corresponding
!! map-space half-spans, garcx and garcy (in degrees) and the quality
!! diagnostic, Q(lam) for this optimal parameter choice. If this process
!! fails for any reason, the failure is alerted by a raised flag, FF, and
!! the other output arguments must then be taken to be meaningless.
!!
!! The diagnostic Q measures the variance over the domain of a local measure
!! of grid distortion. A logarithmic measure of local grid deformation is
!! give by L=log(J^t.J)/2, where J is the mapping Jacobian matrix, dX/dx,
!! where X is the cartesian unit 3-vector representation of the image of the
!! mapping of the map-coordinate 2-vector, x.
!! The Frobenius squared-norm, Trace(L*L), of L is the basis for the simplest
!! (lam=0) definition of the variance of L, but (Trace(L))**2 is another.
!! Here, we weight both contributions, by lam and (1-lam) respectively, with
!! 0 <= lam <1, to compute the variance Q(lam,a,k), and search for the (a,k)
!! that minimizes this Q.
!!
!! The domain averages are computed by double Gauss-Legendre quadrature (i.e.,
!! in both the x and y directions), but restricted to a mere quadrant of the
!! domain (since bilateral symmetry pertains across both domain medians,
!! yielding a domain mean L that is strictly diagonal.
!!
!! @param[in] lam
!! @param[in] garcx map-space half-spans
!! @param[in] garcy map-space half-spans
!! @param[out] a Extended Schmidt Gnomonic parameter
!! @param[out] k Extended Schmidt Gnomonic parameter
!! @param[out] marcx
!! @param[out] marcy
!! @param[out] q
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine bestesg_geo(lam,garcx,garcy, a,k,marcx,marcy,q,ff)!   [bestesg_geo]
use pietc, only: u5,o5,s18,s36,s54,s72,ms18,ms36,ms54,ms72
use pmat,  only: inv
use psym2, only: chol2
implicit none
real(dp),intent(in ):: lam,garcx,garcy
real(dp),intent(out):: a,k,marcx,marcy,q
logical ,intent(out):: FF
integer(spi),parameter     :: nit=200,mit=20,ngh=25
real(dp)    ,parameter     :: u2o5=u2*o5,&
                              f18=u2o5*s18,f36=u2o5*s36,f54=u2o5*s54,&
                              f72=u2o5*s72,mf18=-f18,mf36=-f36,mf54=-f54,&
                              mf72=-f72,& !<- (Fourier transform coefficients)
                              r=0.001_dp,rr=r*r,dang=pi2*o5,crit=1.e-14_dp
real(dp),dimension(ngh)    :: wg,xg
real(dp),dimension(0:4,0:4):: em5 ! <- Fourier matrix for 5 points
real(dp),dimension(0:4)    :: qs
real(dp),dimension(2,0:4)  :: aks,mas
real(dp),dimension(2,2)    :: basis0,basis,hess,el,gadak,gadma,madga,madak
real(dp),dimension(2)      :: ak,dak,dma,vec2,grad,qdak,qdma,ga,ma,gat
real(dp)                   :: s,tgarcx,tgarcy,asp,ang
integer(spi)               :: i,it
logical                    :: flip
data em5/o5,u2o5,  u0,u2o5,  u0,& ! <-The Fourier matrix for 5 points. Applied
         o5, f18, f72,mf54, f36,& ! to the five 72-degree spaced values in a
         o5,mf54, f36, f18,mf72,& ! column-vector, the product vector has the
         o5,mf54,mf36, f18, f72,& ! components, wave-0, cos and sin wave-1,
         o5, f18,mf72,mf54,mf36/  ! cos and sin wave-2.
! First guess upper-triangular basis regularizing the samples used to
! estimate the Hessian of q by finite differencing:
data basis0/0.1_dp,u0,  0.3_dp,0.3_dp/
ff=lam<u0 .or. lam>=u1
if(ff)then; print'("In bestesg_geo; lam out of range")';return; endif
ff= garcx<=u0 .or. garcy<=u0
if(ff)then
   print'("In bestesg_geo; a nonpositive domain parameter, garcx or garcy")'
   return
endif
call gaulegh(ngh,xg,wg)! <- Prepare Gauss-Legendre nodes xg and weights wg
flip=garcy>garcx
if(flip)then; tgarcx=garcy; tgarcy=garcx! <- Switch
else        ; tgarcx=garcx; tgarcy=garcy! <- Don't switch
endif
ga=(/tgarcx,tgarcy/)
asp=tgarcy/tgarcx
basis=basis0

call guessak_geo(asp,tgarcx,ak)
ma=ga*dtor*0.9_dp ! Shrink first estimate, to start always within bounds

! Perform a Newton iteration (except with imperfect Hessian) to find the
! parameter vector, ak, at which the derivative of Q at constant geographical
! semi-axes, ga, as given, goes to zero. The direct evaluation of the
! Q-derivative at constant ma (which is what is actually computed in
! get_meanq) therefore needs modification to obtain Q-derivarive at constant
! ga:
! dQ/d(ak)|_ga = dQ/d(ak)|_ma - dQ/d(ma)|_ak*d(ma)/d(ga)|_ak*d(ga)/d(ak)|_ma
!
! Since the Hessian evaluation is only carried out at constant map-space
! semi-axes, ma, it is not ideal for this problem; consequently, the allowance
! of newton iterations, nit, is much more liberal than we allow for the
! companion routine, bestesg_map, where the constant ma condition WAS
! appropriate.
do it=1,nit
   call get_meanq(ngh,lam,xg,wg,ak,ma,q,qdak,qdma,gat,gadak,gadma,ff)
   if(ff)return
   madga=gadma; call inv(madga)! <- d(ma)/d(ga)|_ak ("at constant ak")
   madak=-matmul(madga,gadak)
   qdak=qdak+matmul(qdma,madak)! dQ/d(ak)|_ga
   if(it<=mit)then ! <- Only recompute aks if the basis is new
! Place five additional sample points around the stencil-ellipse:
      do i=0,4
         ang=i*dang                     ! steps of 72 degrees
         vec2=(/cos(ang),sin(ang)/)*r   ! points on a circle of radius r ...
         dak=matmul(basis,vec2)
         dma=matmul(madak,dak)
         aks(:,i)=ak+dak ! ... become points on an ellipse.
         mas(:,i)=ma+dma
      enddo
      call get_meanq(5,ngh,lam,xg,wg,aks,mas, qs,ff)
   endif
   grad=matmul(qdak,basis)/q ! <- New grad estimate, accurate to near roundoff
   if(it<mit)then ! Assume the hessian will have significantly changed
! Recover an approximate Hessian estimate, normalized by q, from all
! 6 samples, q, qs. These are wrt the basis, not (a,k) directly. Note that
! this finite-difference Hessian is wrt q at constant ak, which is roughly
! approximating the desired Hessian wrt q at constant ga; the disparity
! is the reason that we allow more iterations in this routine than in
! companion routine bestesg_map, since we cannot achieve superlinear
! convergence in the present case.
! Make qs the 5-pt discrete Fourier coefficients of the ellipse pts:
      qs=matmul(em5,qs)/q
!      grad=qs(1:2)/r ! Old finite difference estimate is inferior to new grad
      qs(0)=qs(0)-u1! <- cos(2*ang) coefficient relative to the central value.
      hess(1,1)=qs(0)+qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
      hess(1,2)=qs(4)      ! <- sin(2*ang) coefficient
      hess(2,1)=qs(4)      !
      hess(2,2)=qs(0)-qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
      hess=hess*u2/rr      ! <- rr is r**2

! Perform a Cholesky decomposition of the hessian:
      call chol2(hess,el,ff)
      if(ff)then
         print'(" In bestESG_geo, Hessian is not positive; cholesky fails")'
         return
      endif
! Invert the cholesky factor in place:
      el(1,1)=u1/el(1,1); el(2,2)=u1/el(2,2); el(2,1)=-el(2,1)*el(1,1)*el(2,2)
   endif
! Estimate a Newton step towards the minimum of function Q(a,k):
   vec2=-matmul(transpose(el),matmul(el,grad))
   dak=matmul(basis,vec2)
   gat=gat+matmul(gadak,dak)! <- increment ga
   ak=ak+dak! <- increment the parameter vector estimate
   dma=-matmul(madga,gat-ga)
   ma=ma+dma

! Use the inverse cholesky factor to re-condition the basis. This is to make
! the next stencil-ellipse more closely share the shape of the elliptical
! contours of Q near its minumum -- essentially a preconditioning of the
! numerical optimization:
   if(it<mit)basis=matmul(basis,transpose(el))

   s=sqrt(dot_product(dak,dak))
   if(s<crit)exit ! <-Sufficient convergence of the Newton iteration
enddo
if(it>nit)print'("WARNING; Relatively inferior convergence in bestesg_geo")'
a=ak(1)
k=ak(2)
if(flip)then; marcx=ma(2); marcy=ma(1)! Remember to switch back
else        ; marcx=ma(1); marcy=ma(2)! don't switch
endif
end subroutine bestesg_geo

!> Get the best Extended Schmidt Gnomonic parameter, (a,k), for the
!! given map-coordinate half-spans, marcx and marcy, as well as the
!! corresponding geographical half-spans, garcx and garcy (in degrees)
!! and the quality diagnostic, Q(lam) for this optimal parameter
!! choice. If this process fails for any reason, the failure is
!! alerted by a raised flag, FF, and the other output arguments must
!! then be taken to be meaningless.
!!
!! The diagnostic Q measures the variance over the domain of a local
!! measure of grid distortion. A logarithmic measure of local grid
!! deformation is give by L=log(J^t.J)/2, where J is the mapping
!! Jacobian matrix, dX/dx, where X is the cartesian unit 3-vector
!! representation of the image of the mapping of the map-coordinate
!! 2-vector, x.  The Frobenius squared-norm, Trace(L*L), of L is the
!! basis for the simplest (lam=0) definition of the variance of L, but
!! (Trace(L))**2 is another.  Here, we weight both contributions, by
!! lam and (1-lam) respectively, with 0 <= lam <1, to compute the
!! variance Q(lam,a,k), and search for the (a,k) that minimizes this
!! Q.
!!
!! The domain averages are computed by double Gauss-Legendre
!! quadrature (i.e., in both the x and y directions), but restricted
!! to a mere quadrant of the domain (since bilateral symmetry pertains
!! across both domain medians, yielding a domain mean L that is
!! strictly diagonal.
!!
!! @param[in] lam
!! @param[in] marcx map-coordinate half-spans
!! @param[in] marcy map-coordinate half-spans
!! @param[out] a Extended Schmidt Gnomonic parameter
!! @param[out] k Extended Schmidt Gnomonic parameter
!! @param[out] garcx geographical half-spans
!! @param[out] garcy geographical half-spans
!! @param[out] q
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine bestesg_map(lam,marcx,marcy, a,k,garcx,garcy,q,ff)   ![bestesg_map]
use pietc, only: u5,o5,s18,s36,s54,s72,ms18,ms36,ms54,ms72
use psym2, only: chol2
implicit none
real(dp),intent(in ):: lam,marcx,marcy
real(dp),intent(out):: a,k,garcx,garcy,q
logical ,intent(out):: FF
integer(spi),parameter     :: nit=25,mit=7,ngh=25
real(dp),parameter         :: u2o5=u2*o5,                                &
                              f18=u2o5*s18,f36=u2o5*s36,f54=u2o5*s54,    &
                              f72=u2o5*s72,mf18=-f18,mf36=-f36,mf54=-f54,&
                              mf72=-f72,& !<- (Fourier)
                              r=0.001_dp,rr=r*r,dang=pi2*o5,crit=1.e-12_dp
real(dp),dimension(ngh)    :: wg,xg
real(dp),dimension(0:4,0:4):: em5 ! <- Fourier matrix for 5 points
real(dp),dimension(0:4)    :: qs ! <-Sampled q, its Fourier coefficients
real(dp),dimension(2,0:4)  :: aks,mas! <- tiny elliptical array of ak
real(dp),dimension(2,2)    :: basis0,basis,hess,el,gadak,gadma
real(dp),dimension(2)      :: ak,dak,vec2,grad,qdak,qdma,ga,ma
real(dp)                   :: s,tmarcx,tmarcy,asp,ang
integer(spi)               :: i,it
logical                    :: flip
data em5/o5,u2o5,  u0,u2o5,  u0,& ! <-The Fourier matrix for 5 points. Applied
         o5, f18, f72,mf54, f36,& ! to the five 72-degree spaced values in a
         o5,mf54, f36, f18,mf72,& ! column-vector, the product vector has the
         o5,mf54,mf36, f18, f72,& ! components, wave-0, cos and sin wave-1,
         o5, f18,mf72,mf54,mf36/  ! cos and sin wave-2.
! First guess upper-triangular basis regularizing the samples used to
! estimate the Hessian of q by finite differencing:
data basis0/0.1_dp,u0,  0.3_dp,0.3_dp/
ff=lam<u0 .or. lam>=u1
if(ff)then; print'("In bestesg_map; lam out of range")';return; endif
ff= marcx<=u0 .or. marcy<=u0
if(ff)then
   print'("In bestesg_map; a nonpositive domain parameter, marcx or marcy")'
   return
endif
call gaulegh(ngh,xg,wg)
flip=marcy>marcx
if(flip)then; tmarcx=marcy; tmarcy=marcx! <- Switch
else        ; tmarcx=marcx; tmarcy=marcy! <- Don't switch
endif
ma=(/tmarcx,tmarcy/); do i=0,4; mas(:,i)=ma; enddo
asp=tmarcy/tmarcx
basis=basis0

call guessak_map(asp,tmarcx,ak)

do it=1,nit
   call get_meanq(ngh,lam,xg,wg,ak,ma,q,qdak,qdma,ga,gadak,gadma,ff)
   if(ff)return
   if(it<=mit)then
! Place five additional sample points around the stencil-ellipse:
      do i=0,4
         ang=i*dang                     ! steps of 72 degrees
         vec2=(/cos(ang),sin(ang)/)*r   ! points on a circle of radius r ...
         aks(:,i)=ak+matmul(basis,vec2) ! ... become points on an ellipse.
      enddo
      call get_meanq(5,ngh,lam,xg,wg,aks,mas, qs,ff)
   endif
   grad=matmul(qdak,basis)/q ! <- New grad estimate, accurate to near roundoff
   if(it<mit)then
! Recover Hessian estimate, normalized by q, from all
! 6 samples, q, qs. These are wrt the basis, not (a,k) directly.
! The Hessian estimate uses a careful finite method, which is accurate
! enough. The gradient is NOT estimated by finite differences because we
! need the gradient to be accurate to near roundoff levels in order that
! the converged Newton iteration is a precise solution. 
! Make qs the 5-pt discrete Fourier coefficients of the ellipse pts:
      qs=matmul(em5,qs)/q
!      grad=qs(1:2)/r ! Old finite difference estimate is inferior to new grad
      qs(0)=qs(0)-u1 !<- cos(2*ang) coefficient relative to the central value.
      hess(1,1)=qs(0)+qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
      hess(1,2)=qs(4)      ! <- sin(2*ang) coefficient
      hess(2,1)=qs(4)      !
      hess(2,2)=qs(0)-qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
      hess=hess*u2/rr      ! <- rr is r**2

! Perform a Cholesky decomposition of the hessian:
      call chol2(hess,el,ff)
      if(ff)then
         print'(" In bestESG_map, hessian is not positive; cholesky fails")'
         return
      endif
! Invert the cholesky factor in place:
      el(1,1)=u1/el(1,1); el(2,2)=u1/el(2,2); el(2,1)=-el(2,1)*el(1,1)*el(2,2)
   endif
! Estimate a Newton step towards the minimum of function Q(a,k):
   vec2=-matmul(transpose(el),matmul(el,grad))
   dak=matmul(basis,vec2)
   ga=ga+matmul(gadak,dak)! <- increment ga
   ak=ak+dak! <- increment the parameter vector estimate

! Use the inverse cholesky factor to re-condition the basis. This is to make
! the next stencil-ellipse more closely share the shape of the elliptical
! contours of Q near its minumum -- essentially a preconditioning of the
! numerical optimization:
   if(it<mit)basis=matmul(basis,transpose(el))

   s=sqrt(dot_product(dak,dak))
   if(s<crit)exit ! <-Sufficient convergence of the Newton iteration
enddo
if(it>nit)print'("WARNING; Relatively inferior convergence in bestesg_map")'
a=ak(1)
k=ak(2)
if(flip)then; garcx=ga(2); garcy=ga(1)! Remember to switch back
else        ; garcx=ga(1); garcy=ga(2)! don't switch
endif
end subroutine bestesg_map

!> Use a and k as the parameters of an Extended Schmidt-transformed
!! Gnomonic (ESG) mapping centered at (plat,plon) and twisted about
!! this center by an azimuth angle of pazi counterclockwise (these
!! angles in radians).
!!
!! Assume the radius of the earth is unity, and using the central
!! mapping point as the coordinate origin, set up the grid with
!! central x-spacing delx and y-spacing dely. The grid index location
!! of the left-lower corner of the domain is (lx,ly) (typically both
!! NEGATIVE). The numbers of the grid spaces in x and y directions are
!! nx and ny.  (Note that, for a centered rectangular grid lx and ly
!! are negative and, in magnitude, half the values of nx and ny
!! respectively.)  Return the latitude and longitude, in radians
!! again, of the grid points thus defined in the arrays, glat and
!! glon, and return a rectangular array, garea, of dimensions nx-1 by
!! ny-1, that contains the areas of each of the grid cells
!!
!! If all goes well, return a lowered failure flag, ff=.false. .  But
!! if, for some reason, it is not possible to complete this task,
!! return the raised failure flag, ff=.TRUE. .
!!
!! @param[in] lx center-relative grid index in x of left edge of the domain
!! @param[in] ly center-relative grid index in y of lower edge of the domain
!! @param[in] nx number of grid spaces in x
!! @param[in] ny number of grid spaces in y
!! @param[in] A parameter of the ESG mapping centered at (plat,plon)
!! @param[in] K parameter of the ESG mapping centered at (plat,plon)
!! @param[in] plat latitude of projection center of mapping (radians)
!! @param[in] plon longitude of projection center of mapping (radians)
!! @param[in] pazi azimuth of orientation of mapping at its center
!! @param[in] delx central x-spacing of the grid (radians)
!! @param[in] dely central y-spacing of the grid (radians)
!! @param[out] glat grid points' latitudes
!! @param[out] glon grid points' longitudes
!! @param[out] garea array of grid-cell areas (steradians)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_rr(lx,ly,nx,ny,A,K,plat,plon,pazi, & !       [hgrid_ak_rr]
     delx,dely,  glat,glon,garea, ff)
use pmat4, only: sarea
use pmat5, only: ctogr
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(3,2):: xcd
real(dp),dimension(3)  :: xc
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,&
                          rlat,drlata,drlatb,drlatc,    &
                          rlon,drlona,drlonb,drlonc
integer(spi)           :: ix,iy,mx,my
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)
mx=lx+nx ! Index of the 'right' edge of the rectangular grid
my=ly+ny ! Index of the 'top' edge of the rectangular grid
do iy=ly,my
   xm(2)=iy*dely
   do ix=lx,mx
      xm(1)=ix*delx
      call xmtoxc_ak(a,k,xm,xc,xcd,ff)
      if(ff)return
      xcd=matmul(prot,xcd)
      xc =matmul(prot,xc )
      call ctogr(xc,glat(ix,iy),glon(ix,iy))
   enddo
enddo

! Compute the areas of the quadrilateral grid cells:
do iy=ly,my-1
   do ix=lx,mx-1
      rlat  =glat(ix  ,iy  )
      drlata=glat(ix+1,iy  )-rlat
      drlatb=glat(ix+1,iy+1)-rlat
      drlatc=glat(ix  ,iy+1)-rlat
      rlon  =glon(ix  ,iy  )
      drlona=glon(ix+1,iy  )-rlon
      drlonb=glon(ix+1,iy+1)-rlon
      drlonc=glon(ix  ,iy+1)-rlon
! If 'I' is the grid point (ix,iy), 'A' is (ix+1,iy); 'B' is (ix+1,iy+1)
! and 'C' is (ix,iy+1), then the area of the grid cell IABC is the sum of
! the areas of the traingles, IAB and IBC (the latter being the negative
! of the signed area of triangle, ICB):
      garea(ix,iy)=sarea(rlat, drlata,drlona, drlatb,drlonb) &
                  -sarea(rlat, drlatc,drlonc, drlatb,drlonb)
   enddo
enddo
end subroutine hgrid_ak_rr

!> Use a and k as the parameters of an extended Schmidt-transformed
!! gnomonic (ESG) mapping centered at (plat,plon) and twisted about
!! this center by an azimuth angle of pazi counterclockwise (these
!! angles in radians).
!!
!! Using the central mapping point as the coordinate origin, set up
!! the grid with central x-spacing delx and y-spacing dely in
!! nondimensional units, (i.e., as if the earth had unit radius) and
!! with the location of the left- lower corner of the grid at
!! center-relative grid index pair, (lx,ly) and with the number of the
!! grid spaces in x and y directions given by nx and ny.  (Note that,
!! for a centered rectangular grid lx and ly are negative and, in
!! magnitude, half the values of nx and ny respectively.)  Return the
!! latitude and longitude, again, in radians, of the grid pts thus
!! defined in the arrays, glat and glon; return a rectangular array,
!! garea, of dimensions nx-1 by ny-1, that contains the areas of each
!! of the grid cells in nondimensional "steradian" units.
!!
!! In this version, these grid cell areas are computed by 2D
!! integrating the scalar jacobian of the transformation, using a
!! 4th-order centered scheme.  The estimated grid steps, dx and dy,
!! are returned at the grid cell edges, using the same 4th-order
!! scheme to integrate the 1D projected jacobian.  The angles,
!! relative to local east and north, are returned respectively as
!! angle_dx and angle_dy at grid cell corners, in radians
!! counterclockwise.
!!
!! if all goes well, return a .FALSE. failure flag, ff. If, for some
!! reason, it is not possible to complete this task, return the
!! failure flag as .TRUE.
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain
!! @param[in] nx numbers of the grid spaces in x
!! @param[in] ny numbers of the grid spaces in y
!! @param[in] a Extended Schmidt Gnomonic parameter
!! @param[in] k Extended Schmidt Gnomonic parameter
!! @param[in] plat latitude of the projection center of the mapping (radians)
!! @param[in] plon longitude of the projection center of the mapping (radians)
!! @param[in] pazi azimuth of the orientation of the mapping at its center
!! @param[in] delx central x-spacing of the grid (radians) 
!! @param[in] dely central y-spacing of the grid (radians)
!! @param[out] glat grid points' latitudes (radians)
!! @param[out] glon grid points' longitudes (radians)
!! @param[out] garea array of grid-cell areas (steradians)
!! @param[out] dx grid steps in x at grid cell edges (radians)
!! @param[out] dy grid steps in y at grid cell edges (radians)
!! @param[out] angle_dx x angles relative to local east (radians)
!! @param[out] angle_dy y angles relative to local north (radians)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_rr_c(lx,ly,nx,ny,a,k,plat,plon,pazi, & !     [hgrid_ak_rr]
                    delx,dely,  glat,glon,garea,dx,dy,angle_dx,angle_dy, ff)
use pmat4, only: cross_product,triple_product
use pmat5, only: ctogr
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
real(dp),dimension(lx:lx+nx-1,ly:ly+ny  ),intent(out):: dx
real(dp),dimension(lx:lx+nx  ,ly:ly+ny-1),intent(out):: dy
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: angle_dx,angle_dy
logical,                                  intent(out):: ff
real(dp),dimension(lx-1:lx+nx+1,ly-1:ly+ny+1):: gat ! Temporary area array
real(dp),dimension(lx-1:lx+nx+1,ly  :ly+ny  ):: dxt ! Temporary dx array
real(dp),dimension(lx  :lx+nx  ,ly-1:ly+ny+1):: dyt ! Temporary dy array
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(3,2):: xcd,eano
real(dp),dimension(2,2):: xcd2
real(dp),dimension(3)  :: xc,east,north
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,delxy
integer(spi)           :: ix,iy,mx,my,lxm,lym,mxp,myp
delxy=delx*dely
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)

mx=lx+nx ! Index of the 'right' edge of the rectangular grid
my=ly+ny ! Index of the 'top' edge of the rectangular grid
lxm=lx-1; mxp=mx+1 ! Indices of extra left and right edges
lym=ly-1; myp=my+1 ! Indices of extra bottom and top edges

!-- main body of horizontal grid:
do iy=ly,my
   xm(2)=iy*dely
   do ix=lx,mx
      xm(1)=ix*delx
      call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
      xcd=matmul(prot,xcd)
      xc =matmul(prot,xc )
      call ctogr(xc,glat(ix,iy),glon(ix,iy))
      east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
      north=cross_product(xc,east)
      eano(:,1)=east; eano(:,2)=north
      xcd2=matmul(transpose(eano),xcd)
      angle_dx(ix,iy)=atan2( xcd2(2,1),xcd2(1,1))
      angle_dy(ix,iy)=atan2(-xcd2(1,2),xcd2(2,2))
      dxt(ix,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
      dyt(ix,iy)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
      gat(ix,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
   enddo
enddo

!-- extra left edge, gat, dxt only:
xm(1)=lxm*delx
do iy=ly,my
   xm(2)=iy*dely
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dxt(lxm,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
   gat(lxm,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra right edge, gat, dxt only:
xm(1)=mxp*delx
do iy=ly,my
   xm(2)=iy*dely
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dxt(mxp,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
   gat(mxp,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra bottom edge, gat, dyt only:
xm(2)=lym*dely
do ix=lx,mx
   xm(1)=ix*delx
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dyt(ix,lym)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
   gat(ix,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra top edge, gat, dyt only:
xm(2)=myp*dely
do ix=lx,mx
   xm(1)=ix*delx
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dyt(ix,myp)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
   gat(ix,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

! Extra four corners, gat only:
xm(2)=lym*dely
!-- extra bottom left corner:
xm(1)=lxm*delx
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(lxm,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- extra bottom right corner:
xm(1)=mxp*delx
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(mxp,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

xm(2)=myp*dely
!-- extra top left corner:
xm(1)=lxm*delx
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(lxm,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- extra top right corner:
xm(1)=mxp*delx
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(mxp,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- 4th-order averaging over each central interval using 4-pt. stencils:
dx            =(13_dp*(dxt(lx :mx-1,:)+dxt(lx+1:mx ,:)) &
                  -(dxt(lxm:mx-2,:)+dxt(lx+2:mxp,:)))/24_dp
dy            =(13_dp*(dyt(:,ly :my-1)+dyt(:,ly+1:my )) &
                  -(dyt(:,lym:my-2)+dyt(:,ly+2:myp)))/24_dp
gat(lx:mx-1,:)=(13_dp*(gat(lx :mx-1,:)+gat(lx+1:mx ,:)) &
                  -(gat(lxm:mx-2,:)+gat(lx+2:mxp,:)))/24_dp
garea         =(13_dp*(gat(lx:mx-1,ly :my-1)+gat(lx:mx-1,ly+1:my )) &
                  -(gat(lx:mx-1,lym:my-2)+gat(lx:mx-1,ly+2:myp)))/24_dp
end subroutine hgrid_ak_rr_c

!> Use a and k as the parameters of an Extended Schmidt-transformed
!! Gnomonic (ESG) mapping centered at (plat,plon) and twisted about
!! this center by an azimuth angle of pazi counterclockwise (these
!! angles in radians).
!!
!! Assume the radius of the earth is unity, and using the central
!! mapping point as the coordinate origin, set up the grid with
!! central x-spacing delx and y-spacing dely. The grid index location
!! of the left-lower corner of the domain is (lx,ly) (typically both
!! NEGATIVE).  The numbers of the grid spaces in x and y directions
!! are nx and ny.  (Note that, for a centered rectangular grid lx and
!! ly are negative and, in magnitude, half the values of nx and ny
!! respectively.)  Return the unit cartesian vectors xc of the grid
!! points and their jacobian matrices xcd wrt the map coordinates, and
!! return a rectangular array, garea, of dimensions nx-1 by ny-1, that
!! contains the areas of each of the grid cells
!!
!! If all goes well, return a lowered failure flag, ff=.false. .  But
!! if, for some reason, it is not possible to complete this task,
!! return the raised failure flag, ff=.TRUE. .
!!
!! @param lx center-relative x grid index for left edge of the domain
!! @param ly center-relative y grid index for lower edge of the domain
!! @param nx numbers of the grid spaces in x
!! @param ny numbers of the grid spaces in y
!! @param a parameters of the ESG mapping centered at (plat,plon)
!! @param k parameters of the ESG mapping centered at (plat,plon)
!! @param plat latitude of the projection center of the mapping (radians)
!! @param plon longitude of the projection center of the mapping (radians)
!! @param pazi azimuth of orientation of mapping at its center 
!! @param delx central x-spacing of the grid (in radians)
!! @param dely central y-spacing of the grid (in radians)
!! @param xc Earth-centered unit cartesian 3-vectors at each grid point
!! @param xcd Jacobian matrices, d(xc)/d(xm), at each grid point
!! @param garea rectangular array of grid-cell areas (steradians)
!! @param ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_rc(lx,ly,nx,ny,A,K,plat,plon,pazi, & !       [hgrid_ak_rc]
     delx,dely, xc,xcd,garea, ff)
use pmat4, only: sarea
use pmat5, only: ctogr
implicit none
integer(spi),intent(in ):: lx,ly,nx,ny
real(dp),    intent(in ):: a,k,plat,plon,pazi,delx,dely
real(dp),dimension(3,  lx:lx+nx  ,ly:ly+ny  ),intent(out):: xc
real(dp),dimension(3,2,lx:lx+nx  ,ly:ly+ny  ),intent(out):: xcd
real(dp),dimension(    lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                      intent(out):: ff
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,                  &
                          rlat,rlata,rlatb,rlatc,drlata,drlatb,drlatc,    &
                          rlon,rlona,rlonb,rlonc,drlona,drlonb,drlonc
integer(spi)           :: ix,iy,mx,my
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)
mx=lx+nx ! Index of the 'right' edge of the rectangular grid
my=ly+ny ! Index of the 'top' edge of the rectangular grid
do iy=ly,my
   xm(2)=iy*dely
   do ix=lx,mx
      xm(1)=ix*delx
      call xmtoxc_ak(a,k,xm,xc(:,ix,iy),xcd(:,:,ix,iy),ff)
      if(ff)return
      xcd(:,:,ix,iy)=matmul(prot,xcd(:,:,ix,iy))
      xc (:,  ix,iy)=matmul(prot,xc (:,  ix,iy))
   enddo
enddo

! Compute the areas of the quadrilateral grid cells:
do iy=ly,my-1
   do ix=lx,mx-1
      call ctogr(xc(:,ix  ,iy  ),rlat ,rlon )
      call ctogr(xc(:,ix+1,iy  ),rlata,rlona)
      call ctogr(xc(:,ix+1,iy+1),rlatb,rlonb)
      call ctogr(xc(:,ix  ,iy+1),rlatc,rlonc)
      drlata=rlata-rlat; drlona=rlona-rlon
      drlatb=rlatb-rlat; drlonb=rlonb-rlon
      drlatc=rlatc-rlat; drlonc=rlonc-rlon

! If 'I' is the grid point (ix,iy), 'A' is (ix+1,iy); 'B' is (ix+1,iy+1)
! and 'C' is (ix,iy+1), then the area of the grid cell IABC is the sum of
! the areas of the triangles, IAB and IBC (the latter being the negative
! of the signed area of triangle, ICB):
      garea(ix,iy)=sarea(rlat, drlata,drlona, drlatb,drlonb) &
                  -sarea(rlat, drlatc,drlonc, drlatb,drlonb)
   enddo
enddo
end subroutine hgrid_ak_rc

!> Use a and k as the parameters of an Extended Schmidt-transformed
!! Gnomonic (ESG) mapping centered at (pdlat,pdlon) and twisted about
!! this center by an azimuth angle of pdazi counterclockwise (these
!! angles in degrees).
!!
!! Like hgrid_ak_rr, return the grid points' lats and lons, except
!! that here the angles are returned in degrees. Garea, the area of
!! each grid cell, is returned as in hgrid_ak_rr, and a failure flag,
!! ff, raised when a failure occurs anywhere in these calculations.
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain
!! @param[in] nx number of the grid spaces in x
!! @param[in] ny number of the grid spaces in y
!! @param[in] a parameter of an ESG mapping
!! @param[in] k parameter of an ESG mapping
!! @param[in] pdlat degrees latitude of the projection center of mapping
!! @param[in] pdlon degrees longitude of the projection center of mapping
!! @param[in] pdazi degrees azimuth of orientation of mapping at its center 
!! @param[in] delx central x-spacing of the grid (in radians)
!! @param[in] dely central y-spacing of the grid (in radians)
!! @param[out] gdlat array of grid point latitudes (in degrees)
!! @param[out] gdlon array of grid point longitudes (in dgrees)
!! @param[out] garea array of grid cell areas (in steradians)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_dd(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, & !    [hgrid_ak_dd]
     delx,dely,  gdlat,gdlon,garea, ff)
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: A,K,pdlat,pdlon,&
                                                        pdazi,delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: gdlat,gdlon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
real(dp):: plat,plon,pazi
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
call hgrid_ak_rr(lx,ly,nx,ny,A,K,plat,plon,pazi, &
    delx,dely,   gdlat,gdlon,garea, ff)
if(ff)return
gdlat=gdlat*rtod ! Convert these angles from radians to degrees
gdlon=gdlon*rtod !
end subroutine hgrid_ak_dd

!> Like hgrid_ak_rr_c, except all the angle arguments (but not
!! delx,dely) are in degrees instead of radians.
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain 
!! @param[in] nx numbers of the grid spaces in x
!! @param[in] ny numbers of the grid spaces in y
!! @param[in] a parameters of an ESG mapping
!! @param[in] k parameters of an ESG mapping
!! @param[in] pdlat latitude defining projection center of the mapping
!! @param[in] pdlon longitude defining projection center of the mapping
!! @param[in] pdazi azimuth of the orientation of the mapping at its center 
!! @param[in] delx central x-spacing of the grid (in radians)
!! @param[in] dely central y-spacing of the grid (in radians)
!! @param[out] gdlat array of grid point degree-latitudes 
!! @param[out] gdlon array of grid point degree-longitudes 
!! @param[out] garea array of grid-cell areas (steradians)
!! @param[out] dx step sizes of the grid-cell edges in x (earth radius=1 unit)
!! @param[out] dy step sizes of the grid-cell edges in y (earth radius=1 unit)
!! @param[out] dangle_dx azimuth rotation of the x grid steps, dx (degrees)
!! @param[out] dangle_dy azimuth rotation of the y grid steps, dy (degrees)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_dd_c(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, &!   [hgrid_ak_dd]
     delx,dely,  gdlat,gdlon,garea,dx,dy,dangle_dx,dangle_dy, ff)
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,pdlat,pdlon,&
                                                        pdazi,delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: gdlat,gdlon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
real(dp),dimension(lx:lx+nx-1,ly:ly+ny  ),intent(out):: dx
real(dp),dimension(lx:lx+nx  ,ly:ly+ny-1),intent(out):: dy
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: dangle_dx,dangle_dy
logical,                                  intent(out):: ff
real(dp):: plat,plon,pazi
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
call hgrid_ak_rr_c(lx,ly,nx,ny,A,K,plat,plon,pazi, &
     delx,dely,  gdlat,gdlon,garea,dx,dy,dangle_dx,dangle_dy, ff)
if(ff)return
gdlat    =gdlat    *rtod ! Convert these angles from radians to degrees
gdlon    =gdlon    *rtod !
dangle_dx=dangle_dx*rtod !
dangle_dy=dangle_dy*rtod !
end subroutine hgrid_ak_dd_c

!> Use a and k as the parameters of an Extended Schmidt-transformed
!! Gnomonic (ESG) mapping centered at (pdlat,pdlon) and twisted about
!! this center by an azimuth angle of pdazi counterclockwise (these
!! angles in degrees).
!!
!! Like hgrid_ak_rx, return the grid points' cartesians xc and
!! Jacobian matrices, xcd. Garea, the area of each grid cell, is also
!! returned as in hgrid_ak_rx, and a failure flag, ff, raised when a
!! failure occurs anywhere in these calculations.
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain
!! @param[in] nx numbers of the grid spaces in x
!! @param[in] ny numbers of the grid spaces in y
!! @param[in] a parameters of an ESG mapping
!! @param[in] k parameters of an ESG mapping
!! @param[in] pdlat degrees latitude of the projection center of the mapping
!! @param[in] pdlon degrees longitude of the projection center of the mapping
!! @param[in] pdazi azimuth of the orientation of the mapping at its center 
!! @param[in] delx central x-spacing of the grid in radians
!! @param[in] dely central y-spacing of the grid in radians
!! @param[out] xc grid points' earth-centered unit cartesians
!! @param[out] xcd Jacobian matrices, d(xc)/d(xm)
!! @param[out] garea array of grid-cell areas (steradians)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_dc(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, & !    [hgrid_ak_dc]
     delx,dely, xc,xcd,garea, ff)
implicit none
integer(spi),intent(in ):: lx,ly,nx,ny
real(dp),    intent(in ):: A,K,pdlat,pdlon,pdazi,delx,dely
real(dp),dimension(3,  lx:lx+nx  ,ly:ly+ny  ),intent(out):: xc
real(dp),dimension(3,2,lx:lx+nx  ,ly:ly+ny  ),intent(out):: xcd
real(dp),dimension(    lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                      intent(out):: ff
real(dp):: plat,plon,pazi
plat=pdlat*dtor
plon=pdlon*dtor
pazi=pdazi*dtor
call hgrid_ak_rc(lx,ly,nx,ny,A,K,plat,plon,pazi, &
    delx,dely,   xc,xcd,garea, ff)
end subroutine hgrid_ak_dc

!> Like hgrid_ak_rr_c except the argument list includes the earth
!! radius, re, and this is used to express the map-space grid
!! increments in the dimensional units, delxre, delyre on entry, and
!! to express the grid cell areas, garea, in dimensional units upon
!! return.
!!
!! The gridded lats and lons, glat and glon, remain in radians.
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain
!! @param[in] nx numbers of the grid spaces in x
!! @param[in] ny numbers of the grid spaces in y
!! @param[in] a parameters of an ESG mapping
!! @param[in] k parameters of an ESG mapping
!! @param[in] plat radians latitude of the projection center of the mapping
!! @param[in] plon radians longitude of the projection center of the mapping
!! @param[in] pazi Azimuth of map orientation at its center
!! @param[in] re earth radius
!! @param[in] delxre map-space grid increments in the dimensional units
!! @param[in] delyre map-space grid increments in the dimensional units
!! @param[out] glat grid points for latitude
!! @param[out] glon grid points for longitude
!! @param[out] garea array of grid-cell areas in dimensional units
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak(lx,ly,nx,ny,a,k,plat,plon,pazi, & !             [hgrid_ak]
     re,delxre,delyre,  glat,glon,garea, ff)
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        re,delxre,delyre
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
real(dp):: delx,dely,rere
delx=delxre/re ! <- set nondimensional grid step delx
dely=delyre/re ! <- set nondimensional grid step dely
call hgrid_ak_rr(lx,ly,nx,ny,a,k,plat,plon,pazi, &
     delx,dely,  glat,glon,garea, ff)
if(ff)return
rere=re*re
garea=garea*rere ! <- Convert from steradians to physical area units.
end subroutine hgrid_ak

!> Like hgrid_ak_rr_c except the argument list includes the earth
!! radius, re, and this is used to express the map-space grid
!! increments in the dimensional units, delxre, delyre on entry, and
!! to express the grid cell areas, garea, and the x- and y- grid
!! steps, dx and dy, in dimensional units upon return.  The gridded
!! lats and lons, glat and glon, remain in radians.  Also, in order
!! for the argument list to remain compatible with an earlier version
!! of this routine, the relative rotations of the steps, dangle_dx and
!! dangle_dy, are returned as degrees instead of radians (all other
!! angles in the argument list, i.e., plat,plon,pazi,glat,glon, remain
!! radians).
!!
!! @param[in] lx center-relative x grid index for left edge of the domain
!! @param[in] ly center-relative y grid index for lower edge of the domain
!! @param[in] nx number of grid spaces in x 
!! @param[in] ny number of grid spaces in y 
!! @param[in] a Extended Schmidt Gnomonic parameter
!! @param[in] k Extended Schmidt Gnomonic parameter
!! @param[in] plat latitude of projection center of the mapping (radians)
!! @param[in] plon longitude of projection center of the mapping (radians)
!! @param[in] pazi Azimuth of map orientation at its center (radians) 
!! @param[in] re earth radius in dimensional length units 
!! @param[in] delxre map-space grid increments in the dimensional units
!! @param[in] delyre map-space grid increments in the dimensional units
!! @param[out] glat gridded lats (radians)
!! @param[out] glon gridded lons (radians)
!! @param[out] garea grid cell areas in dimensional units
!! @param[out] dx x- grid steps in dimensional units
!! @param[out] dy y- grid steps in dimensional units
!! @param[out] dangle_dx azimuth rotations of the steps dx (in degrees)
!! @param[out] dangle_dy azimuth rotations of the steps dy (in degrees)
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine hgrid_ak_c(lx,ly,nx,ny,a,k,plat,plon,pazi, & !           [hgrid_ak]
     re,delxre,delyre,  glat,glon,garea,dx,dy,dangle_dx,dangle_dy, ff)
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        re,delxre,delyre
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
real(dp),dimension(lx:lx+nx-1,ly:ly+ny  ),intent(out):: dx
real(dp),dimension(lx:lx+nx  ,ly:ly+ny-1),intent(out):: dy
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: dangle_dx,dangle_dy
logical,                                  intent(out):: ff
real(dp):: delx,dely,rere
delx=delxre/re ! <- set nondimensional grid step delx
dely=delyre/re ! <- set nondimensional grid step dely
call hgrid_ak_rr_c(lx,ly,nx,ny,a,k,plat,plon,pazi, &
     delx,dely,  glat,glon,garea,dx,dy,dangle_dx,dangle_dy, ff)
if(ff)return
rere=re*re
garea=garea*rere ! <- Convert from steradians to physical area units.
dx=dx*re         ! <- Convert from nondimensional to physical length units.
dy=dy*re         ! <-
dangle_dx=dangle_dx*rtod ! <-Convert from radians to degrees
dangle_dy=dangle_dy*rtod ! <-Convert from radians to degrees
end subroutine hgrid_ak_c

!> This Gauss-Legendre quadrature integrates exactly any even
!! polynomial up to degree m*4-2 in the half-interval [0,1]. This code
!! is liberally adapted from the algorithm given in Press et al.,
!! Numerical Recipes.
!!
!! @param m number of nodes in half-interval
!! @param x nodes and weights
!! @param w nodes and weights
!! @author R. J. Purser
subroutine gaulegh(m,x,w)!                                           [gaulegh]
implicit none
integer(spi),         intent(IN ):: m   ! <- number of nodes in half-interval 
real(dp),dimension(m),intent(OUT):: x,w ! <- nodes and weights
integer(spi),parameter:: nit=8
real(dp),    parameter:: eps=3.e-14_dp
integer(spi)          :: i,ic,j,jm,it,m2,m4p,m4p3
real(dp)              :: z,zzm,p1,p2,p3,pp,z1
m2=m*2; m4p=m*4+1; m4p3=m4p+2
do i=1,m; ic=m4p3-4*i
   z=cos(pih*ic/m4p); zzm=z*z-u1
   do it=1,nit
      p1=u1; p2=u0
      do j=1,m2; jm=j-1; p3=p2; p2=p1; p1=((j+jm)*z*p2-jm*p3)/j; enddo
      pp=m2*(z*p1-p2)/zzm; z1=z; z=z1-p1/pp; zzm=z*z-u1
      if(abs(z-z1) <= eps)exit
   enddo
   x(i)=z; w(i)=-u2/(zzm*pp*pp)
enddo
end subroutine gaulegh

!> Given the map specification (angles in radians), the grid spacing
!! (in map-space units) and the sample lat-lon (in radian), return the
!! the image in map space in a 2-vector in grid units. If the
!! transformation is invalid, return a .true. failure flag.
!!
!! @param[in] a parameters of an ESG mapping
!! @param[in] k parameters of an ESG mapping
!! @param[in] plat radians latitude defining mapping projection center
!! @param[in] plon radians longitude defining mapping projection center
!! @param[in] pazi Aximuth of mapping orientation at its center 
!! @param[in] lat radians latitude of a point to be mapped
!! @param[in] lon radians longitude of a point to be mapped
!! @param[out] xm 2-vector center-relative map-space image of mapped point 
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine gtoxm_ak_rr_m(A,K,plat,plon,pazi,lat,lon,xm,ff)!      [gtoxm_ak_rr]
use pmat5, only: grtoc
implicit none
real(dp),             intent(in ):: a,k,plat,plon,pazi,lat,lon
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
real(dp),dimension(3,3):: prot,azirot
real(dp)               :: clat,slat,clon,slon,cazi,sazi
real(dp),dimension(3)  :: xc
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)

call grtoc(lat,lon,xc)
xc=matmul(transpose(prot),xc)
call xctoxm_ak(a,k,xc,xm,ff)
end subroutine gtoxm_ak_rr_m

!> Given the map specification (angles in radians), the grid spacing
!! (in map-space units) and the sample lat-lon (in radian), return the
!! the image in map space in a 2-vector in grid units. If the
!! transformation is invalid, return a .true. failure flag.
!!
!! @param[in] a parameter of the ESG mapping
!! @param[in] k parameter of the ESG mapping
!! @param[in] plat radians latitude defining mapping projection center
!! @param[in] plon radians longitude defining mapping projection center
!! @param[in] pazi Azimuth of mapping orientation at its center 
!! @param[in] delx central x-spacing of the grid in radians
!! @param[in] dely central y-spacing of the grid in radians
!! @param[in] lat radians latitude of a point to be mapped
!! @param[in] lon radians longitude of a point to be mapped
!! @param[out] xm 2-vector map space image in center-relative grid units 
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine gtoxm_ak_rr_g(A,K,plat,plon,pazi,delx,dely,lat,lon,&! [gtoxm_ak_rr]
     xm,ff)
implicit none
real(dp),             intent(in ):: a,k,plat,plon,pazi,delx,dely,lat,lon
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
call gtoxm_ak_rr_m(A,K,plat,plon,pazi,lat,lon,xm,ff); if(ff)return
xm(1)=xm(1)/delx; xm(2)=xm(2)/dely
end subroutine gtoxm_ak_rr_g

!> Like gtoxm_ak_rr_m, except lat, lon, azimuth, are expressed in degrees.
!!
!! @param[in] a parameter of the ESG mapping
!! @param[in] k parameter of the ESG mapping
!! @param[in] pdlat degrees latitude defining mapping center
!! @param[in] pdlon degrees longitude defining mapping center
!! @param[in] pdazi Azimuth of mapping orientation at its center 
!! @param[in] dlat degrees latitude of point to be mapped
!! @param[in] dlon degrees longitude of point to be mapped
!! @param[out] xm 2-vector center-relative map space image of the point 
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine gtoxm_ak_dd_m(A,K,pdlat,pdlon,pdazi,dlat,dlon,&!      [gtoxm_ak_dd]
     xm,ff)
implicit none
real(dp),             intent(in ):: a,k,pdlat,pdlon,pdazi,dlat,dlon
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
real(dp):: plat,plon,pazi,lat,lon
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
lat=dlat*dtor
lon=dlon*dtor
call gtoxm_ak_rr_m(A,K,plat,plon,pazi,lat,lon,xm,ff)
end subroutine gtoxm_ak_dd_m

!> Like gtoxm_ak_rr_g, except lat, lon, azimuth, are expressed in degrees.
!!
!! @param[in] a parameter of the ESG mapping
!! @param[in] k parameter of the ESG mapping
!! @param[in] pdlat degrees latitude defining mapping projection center
!! @param[in] pdlon degrees longitude defining mapping projection center
!! @param[in] pdazi Azimuth of mapping orientation at its center 
!! @param[in] delx central x-spacing of the grid in radians
!! @param[in] dely central y-spacing of the grid in radians
!! @param[in] dlat degrees latitude of a point to be mapped
!! @param[in] dlon degrees longitude of a point to be mapped
!! @param[out] xm 2-vector image of the point in center-relative grid units 
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine gtoxm_ak_dd_g(A,K,pdlat,pdlon,pdazi,delx,dely,&!      [gtoxm_ak_dd]
dlat,dlon,     xm,ff)
implicit none
real(dp),             intent(in ):: a,k,pdlat,pdlon,pdazi,delx,dely,dlat,dlon
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
real(dp):: plat,plon,pazi,lat,lon
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
lat=dlat*dtor
lon=dlon*dtor
call gtoxm_ak_rr_g(A,K,plat,plon,pazi,delx,dely,lat,lon,xm,ff)
end subroutine gtoxm_ak_dd_g

!> Given the ESG map specified by parameters (A,K) and geographical center and
!! orientation, plat,plon,pazi (radians), and a position, in map-space
!! coordinates given by the 2-vector, xm, return the geographical
!! coordinates, lat and lon (radians). If the transformation is
!! invalid for any reason, return instead with a raised failure flag,
!! FF= .true.
!!
!! @param[in] a parameter of an ESG mapping
!! @param[in] k parameter of an ESG mapping
!! @param[in] plat radians latitude of the projection center of the mapping
!! @param[in] plon radians longitude of the projection center of the mapping
!! @param[in] pazi Azimuth of orientation of the mapping at its center 
!! @param[in] xm center-relative 2-vector map space coordinates of a point 
!! @param[out] lat radians latitude of the point
!! @param[out] lon radians longitude of the point
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine xmtog_ak_rr_m(A,K,plat,plon,pazi,xm,lat,lon,ff)!      [xmtog_ak_rr]
use pmat5, only: ctogr
implicit none
real(dp),             intent(in ):: a,k,plat,plon,pazi
real(dp),dimension(2),intent(in ):: xm
real(dp),             intent(out):: lat,lon
logical,              intent(out):: ff
real(dp),dimension(3,2):: xcd
real(dp),dimension(3,3):: prot,azirot
real(dp)               :: clat,slat,clon,slon,cazi,sazi
real(dp),dimension(3)  :: xc
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xc=matmul(prot,xc)
call ctogr(xc,lat,lon)
end subroutine xmtog_ak_rr_m

!> For an ESG map with parameters, (A,K), and geographical
!! orientation, given by plon,plat,pazi (radians), and given a point
!! in grid-space units as the 2-vector, xm, return the geographical
!! coordinates, lat, lon, (radians) of this point. If instead the
!! transformation is invalid for any reason, then return the raised
!! failure flag, FF=.true.
!!
!! @param[in] a parameters of the ESG mapping
!! @param[in] k parameters of the ESG mapping
!! @param[in] plat radians latitude of the projection center of the mapping
!! @param[in] plon radians longitude of the projection center of the mapping
!! @param[in] pazi Azimuth of the orientation of the mapping at its center 
!! @param[in] delx central x-spacing of the grid in radians
!! @param[in] dely central y-spacing grid point in radians
!! @param[in] xm grid-space 2-vector coordinates of a point to be mapped
!! @param[out] lat radians latitude of the point
!! @param[out] lon radians longitude of the point
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine xmtog_ak_rr_g(A,K,plat,plon,pazi,delx,dely,xm,&!      [xmtog_ak_rr]
     lat,lon,ff)
implicit none
real(dp),             intent(in ):: a,k,plat,plon,pazi,delx,dely
real(dp),dimension(2),intent(in ):: xm
real(dp),             intent(out):: lat,lon
logical,              intent(out):: ff
real(dp),dimension(2):: xmt
xmt(1)=xm(1)*delx ! Convert from grid units to intrinsic map-space units
xmt(2)=xm(2)*dely !
call xmtog_ak_rr_m(A,K,plat,plon,pazi,xmt,lat,lon,ff)
end subroutine xmtog_ak_rr_g

!> Like xmtog_ak_rr_m, except lat, lon, azimuth, are expressed in degrees.
!!
!! @param[in] a parameters of the ESG mapping
!! @param[in] k parameters of the ESG mapping
!! @param[in] pdlat degrees latitude of the projection center of the mapping 
!! @param[in] pdlon degrees longitude of the projection center of the mapping
!! @param[in] pdazi Azimuth of the orientation of the mapping at its center
!! @param[in] xm map space 2-vector coordinates of a point 
!! @param[out] dlat degrees latitude of the point
!! @param[out] dlon degrees longitude of the point
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine xmtog_ak_dd_m(A,K,pdlat,pdlon,pdazi,xm,dlat,dlon,ff)! [xmtog_ak_dd]
use pmat5, only: ctogr
implicit none
real(dp),             intent(in ):: a,k,pdlat,pdlon,pdazi
real(dp),dimension(2),intent(in ):: xm
real(dp),             intent(out):: dlat,dlon
logical,              intent(out):: ff
real(dp):: plat,plon,pazi,lat,lon
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
call xmtog_ak_rr_m(A,K,plat,plon,pazi,xm,lat,lon,ff)
dlat=lat*rtod
dlon=lon*rtod
end subroutine xmtog_ak_dd_m

!> Like xmtog_ak_rr_g, except lat, lon, azimuth, are expressed in degrees.
!!
!! @param[in] a parameters of an ESG mapping
!! @param[in] k parameters of an ESG mapping
!! @param[in] pdlat degrees latitude of projection center of the mapping
!! @param[in] pdlon degrees longitude of projection center of the mapping
!! @param[in] pdazi Azimuth of the mapping orientation about its center 
!! @param[in] delx central x-spacing of the grid in radians
!! @param[in] dely central y-spacing of the grid in radians
!! @param[in] xm map coordinates, in grid units, of a point to be mapped
!! @param[out] dlat degrees latitude of the point
!! @param[out] dlon degrees longitude of the point
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine xmtog_ak_dd_g(A,K,pdlat,pdlon,pdazi,delx,dely,xm,&!   [xmtog_ak_dd]
     dlat,dlon,ff)
implicit none
real(dp),             intent(in ):: a,k,pdlat,pdlon,pdazi,delx,dely
real(dp),dimension(2),intent(in ):: xm
real(dp),             intent(out):: dlat,dlon
logical,              intent(out):: ff
real(dp),dimension(2):: xmt
real(dp)             :: plat,plon,pazi,lat,lon
xmt(1)=xm(1)*delx ! Convert from grid units to intrinsic map-space units
xmt(2)=xm(2)*dely !
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
call xmtog_ak_rr_m(A,K,plat,plon,pazi,xmt,lat,lon,ff)
dlat=lat*rtod
dlon=lon*rtod
end subroutine xmtog_ak_dd_g

end module pesg

