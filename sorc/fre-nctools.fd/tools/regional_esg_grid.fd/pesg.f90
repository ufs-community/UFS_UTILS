!
!                                              ***********************
!                                              *      pesg.f90       *
!                                              *    R. J. Purser     *
!                                              *   NOAA/NCEP/EMC     *
!                                              *     May 2020        * 
!                                              *                     *
!                                              * jim.purser@noaa.gov *                
!                                              ***********************
! Suite of routines to perform the 2-parameter family of Extended
! Schmidt Gnomonic (ESG) regional grid mappings, and to optimize the
! the two parameters, A and K, of those mappings for a given rectangular
! domain's principal (median) semi-arcs with respect to a domain-averaged
! measure of distortion. This criterion is itself endowed with a parameter, 
! lam (for "lambda" in [0,1) ) which gives weight to additional weight
! areal inhomogeneities instead of treating all distortion components
! equally.
!
! DEPENDENCIES
! Libraries: pmat, psym2, pfun
! Modules: pkind, pietc, pietc_s 
!=============================================================================
module pesg
!=============================================================================
use pkind, only: spi,dp
use pietc, only: F,T,u0,u1,u2,o2,rtod,dtor,pih
implicit none
private
public :: xctoxs,xstoxc,xstoxt,xttoxs,xttoxm,zttozm,zmtozt,xctoxm_ak,xmtoxc_ak,&
          getedges,get_qq,get_qofmap,get_bestesg,get_bestesg_inv, &
          hgrid_ak_rr,hgrid_ak_rc,hgrid_ak_dd,hgrid_ak_dc,hgrid_ak
interface xctoxs;         module procedure xctoxs;          end interface
interface xstoxc;         module procedure xstoxc;          end interface
interface xstoxt;         module procedure xstoxt;          end interface
interface xttoxs;         module procedure xttoxs;          end interface
interface xttoxm;         module procedure xttoxm;          end interface
interface zttozm;         module procedure zttozm;          end interface
interface zmtozt;         module procedure zmtozt;          end interface
interface xctoxm_ak;      module procedure xctoxm_ak;       end interface
interface xmtoxc_ak;      module procedure xmtoxc_ak;       end interface
interface getedges;       module procedure getedges;        end interface
interface get_wxy;        module procedure get_wxy;         end interface
interface get_qq;         module procedure get_qqw,get_qqt; end interface
interface get_qofmap;     module procedure get_qofmap;      end interface
interface get_bestesg;    module procedure get_bestesg;     end interface
interface get_bestesgt;   module procedure get_bestesgt;    end interface
interface get_bestesg_inv;module procedure get_bestesg_inv; end interface
interface hgrid_ak_rr
   module procedure hgrid_ak_rr,hgrid_ak_rr_c;              end interface
interface hgrid_ak_rc;    module procedure hgrid_ak_rc;     end interface
interface hgrid_ak_dd
   module procedure hgrid_ak_dd,hgrid_ak_dd_c;              end interface
interface hgrid_ak_dc;    module procedure hgrid_ak_dc;     end interface
interface hgrid_ak 
     module procedure hgrid_ak,hgrid_ak_c;                  end interface
contains

!=============================================================================
subroutine xctoxs(xc,xs)!                                             [xctoxs]
!=============================================================================
! Inverse of xstoxc. I.e., cartesians to stereographic
!=============================================================================
implicit none
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xs
!-----------------------------------------------------------------------------
real(dp):: zp
!=============================================================================
zp=u1+xc(3); xs=xc(1:2)/zp
end subroutine xctoxs

!=============================================================================
subroutine xstoxc(xs,xc,xcd)!                                         [xstoxc]
!=============================================================================
! Standard transformation from polar stereographic map coordinates, xs, to
! cartesian, xc, assuming the projection axis is polar.
! xcd=d(xc)/d(xs) is the jacobian matrix, encoding distortion and metric.
!=============================================================================
use pmat4, only: outer_product
implicit none
real(dp),dimension(2),  intent(in ):: xs
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
!-----------------------------------------------------------------------------
real(dp):: zp
!=============================================================================
zp=u2/(u1+dot_product(xs,xs)); xc(1:2)=xs*zp; xc(3)=zp
xcd=-outer_product(xc,xs)*zp; xcd(1,1)=xcd(1,1)+zp; xcd(2,2)=xcd(2,2)+zp
xc(3)=xc(3)-u1
end subroutine xstoxc

!=============================================================================
subroutine xstoxt(kappa,xs,xt,ff)!                                    [xstoxt]
!=============================================================================
! Inverse of xttoxs.
!=============================================================================
implicit none
real(dp),             intent(in ):: kappa
real(dp),dimension(2),intent(in ):: xs
real(dp),dimension(2),intent(out):: xt
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: s,sc
!=============================================================================
s=kappa*(xs(1)*xs(1)+xs(2)*xs(2)); sc=u1-s
ff=abs(s)>=u1; if(ff)return
xt=u2*xs/sc
end subroutine xstoxt

!==============================================================================
subroutine xttoxs(kappa,xt,xs,xsd,ff)!                                 [xttoxs]
!==============================================================================
! Scaled gnomonic plane xt to standard stereographic plane xs
!==============================================================================
implicit none
real(dp),               intent(in ):: kappa
real(dp),dimension(2),  intent(in ):: xt
real(dp),dimension(2),  intent(out):: xs
real(dp),dimension(2,2),intent(out):: xsd
logical,                intent(out):: ff
!------------------------------------------------------------------------------
real(dp):: s,sp,rsp,rspp,rspps,rspdx,rspdy
!==============================================================================
s=kappa*(xt(1)*xt(1) + xt(2)*xt(2)); sp=u1+s
ff=(sp<=u0); if(ff)return
rsp=sqrt(sp)
rspp=u1+rsp
rspps=rspp**2
xs=xt/rspp
rspdx=kappa*xt(1)/rsp
rspdy=kappa*xt(2)/rsp
xsd(1,1)=u1/rspp -xt(1)*rspdx/rspps
xsd(1,2)=        -xt(1)*rspdy/rspps
xsd(2,1)=        -xt(2)*rspdx/rspps
xsd(2,2)=u1/rspp -xt(2)*rspdy/rspps
end subroutine xttoxs

!=============================================================================
subroutine xttoxm(a,xt,xm,ff)!                                         [xttoxm]
!=============================================================================
! Inverse of xmtoxt
!============================================================================
implicit none
real(dp),             intent(in ):: a
real(dp),dimension(2),intent(in ):: xt
real(dp),dimension(2),intent(out):: xm
logical              ,intent(out):: ff
!-----------------------------------------------------------------------------
integer(spi):: i
!=============================================================================
do i=1,2; call zttozm(a,xt(i),xm(i),ff); if(ff)return; enddo
end subroutine xttoxm

!==============================================================================
subroutine xmtoxt(a,xm,xt,xtd,ff)!                                     [xmtoxt]
!==============================================================================
! Like zmtozt, but for 2-vector xm and xt, and 2*2 diagonal Jacobian xtd
!==============================================================================
implicit none
real(dp),               intent(in ):: a
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(2),  intent(out):: xt
real(dp),dimension(2,2),intent(out):: xtd
logical,                intent(out):: ff
!-----------------------------------------------------------------------------
integer(spi):: i
!==============================================================================
xtd=u0; do i=1,2; call zmtozt(a,xm(i),xt(i),xtd(i,i),ff); if(ff)return; enddo
end subroutine xmtoxt

!=============================================================================
subroutine zttozm(a,zt,zm,ff)!                                        [zttozm]
!=============================================================================
! Inverse of zmtozt
!=============================================================================
implicit none
real(dp),intent(in ):: a,zt
real(dp),intent(out):: zm
logical, intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: ra,razt
!=============================================================================
ff=F
if    (a>u0)then; ra=sqrt( a); razt=ra*zt; zm=atan (razt)/ra
elseif(a<u0)then; ra=sqrt(-a); razt=ra*zt; ff=abs(razt)>=u1; if(ff)return
                                           zm=atanh(razt)/ra
else                                     ; zm=zt
endif
end subroutine zttozm

!==============================================================================
subroutine zmtozt(a,zm,zt,ztd,ff)!                                     [zmtozt]
!==============================================================================
! Evaluate the function, zt = tan(sqrt(A)*z)/sqrt(A), and its derivative, ztd,
! for positive and negative A and for the limiting case, A --> 0
!==============================================================================
implicit none
real(dp),intent(in ):: a,zm
real(dp),intent(out):: zt,ztd
logical, intent(out):: ff
!------------------------------------------------------------------------------
real(dp):: ra
!==============================================================================
ff=f
if    (a>u0)then; ra=sqrt( a); zt=tan (ra*zm)/ra; ff=abs(ra*zm)>=pih
elseif(a<u0)then; ra=sqrt(-a); zt=tanh(ra*zm)/ra
else                         ; zt=zm
endif
ztd=u1+a*zt*zt
end subroutine zmtozt

!=============================================================================
subroutine xctoxm_ak(a,kappa,xc,xm,ff)!                             [xctoxm_ak]
!=============================================================================
! Inverse mapping of xmtoxc_ak. That is, go from given cartesian unit 3-vector,
! xc, to map coordinate 2-vector xm (or return a raised failure flag, FF, if
! the attempt fails).
!=============================================================================
implicit none
real(dp),             intent(in ):: a,kappa
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xm
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(2):: xs,xt
!=============================================================================
ff=F
call xctoxs(xc,xs)
call xstoxt(kappa,xs,xt,ff); if(ff)return
call xttoxm(a,xt,xm,ff)
end subroutine xctoxm_ak

!==============================================================================
subroutine xmtoxc_ak(a,kappa,xm,xc,xcd,ff)!                         [xmtoxc_ak]
!==============================================================================
! Assuming the A-Kappa parameterization of the Extended Schmidt-transformed
! Gnomonic (ESG) mapping, and given a map-space 2-vector, xm, find the
! corresponding cartesian unit 3-vector and its derivative wrt xm, jacobian
! matrix, xcd. If for any reason the mapping cannot be done, return a
! raised failure flag, FF.
!=============================================================================
implicit none
real(dp),               intent(in ):: a,kappa
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(3),  intent(out):: xc
real(dp),dimension(3,2),intent(out):: xcd
logical,                intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: xtd,xsd
real(dp),dimension(2)  :: xt,xs
!=============================================================================
call xmtoxt(a,xm,xt,xtd,ff);     if(ff)return
call xttoxs(kappa,xt,xs,xsd,ff); if(ff)return
xsd=matmul(xsd,xtd)
call xstoxc(xs,xc,xcd)
xcd=matmul(xcd,xsd)
end subroutine xmtoxc_ak

!=============================================================================
subroutine getedges(arcx,arcy,edgex,edgey)!                          [getedges]
!=============================================================================
! For angles (degrees) of the arcs spanning the halfwidths between the
! region's center and its x and y edges, get the two cartesian vectors
! that represent the locations of these edge midpoints in the positive x and y
! directions.
!=============================================================================
implicit none
real(dp),             intent(in ):: arcx,arcy
real(dp),dimension(3),intent(out):: edgex,edgey
!------------------------------------------------------------------------------
real(dp):: cx,sx,cy,sy
!==============================================================================
cx=cos(arcx*dtor); sx=sin(arcx*dtor)
cy=cos(arcy*dtor); sy=sin(arcy*dtor)
edgex=(/sx,u0,cx/); edgey=(/u0,sy,cy/)
end subroutine getedges

!=============================================================================
subroutine get_wxy(nxh,nyh,ncor,wxy)!                                 [get_wxy]
!=============================================================================
! Get the array of weights, wxy, for the positive quadrant of the rectangular
! grid having nxh*2 spaces in "x" and nyh spaces in "y" where it assumed the
! grid is uniform and so the extended trapezoidal integration scheme would
! be the appropriate source for these weights on the full grid. The extended
! scheme used is the one of order ncor. The weights are normalized, so they
! provide a definition of the grid average of the quantity they are applied to.
!=============================================================================
use pmat4, only: outer_product
implicit none
integer(spi),                   intent(in ):: nxh,nyh,ncor
real(dp),dimension(0:nxh,0:nyh),intent(out):: wxy
!-----------------------------------------------------------------------------
integer(spi),parameter     :: dencor0=2,dencor1=12,dencor2=24, &
                              dencor3=720,dencor4=1440! <-denominators
real(dp),dimension(0:nxh)  :: wx
real(dp),dimension(0:nyh)  :: wy
real(dp),dimension(0:ncor) :: cor     ! Becomes the full end-correction vector.
integer(spi),dimension(0:0):: numcor0 ! numerator for uncorrected trapezoidal.
integer(spi),dimension(0:1):: numcor1 ! numerators for common extended trap.
integer(spi),dimension(0:2):: numcor2 ! ..numerators for higher order schemes..
integer(spi),dimension(0:3):: numcor3 !
integer(spi),dimension(0:4):: numcor4 !
data numcor0/  1/
data numcor1/  5,  13/
data numcor2/  9,  28,  23/
data numcor3/251, 897, 633, 739/
data numcor4/475,1902,1104,1586,1413/
!=============================================================================
! Initialize the real end correction coefficients, cor, to make the
! "extended" trapezoidal integration schemes in both directions accurate
! to a higher order (probably ncor+2, formally). These corrections are, in
! some sense, like the discrete analogues of the Euler-Maclaurin formulae.
select case(ncor)
case(0); cor=numcor0; cor=cor/dencor0
case(1); cor=numcor1; cor=cor/dencor1
case(2); cor=numcor2; cor=cor/dencor2
case(3); cor=numcor3; cor=cor/dencor3
case(4); cor=numcor4; cor=cor/dencor4
case default; stop 'In get_wxy; this value of ncor is not valid'
end select
if(ncor<0 .or. ncor>4)stop 'In get_wxy; ncor is out of bounds'
if(ncor>=min(nxh,nyh))stop 'In get_wxy; ncor is too large for this small grid'
! the wx and wy are the weight coefficients for an unnormalized
! extended trapezoidal integration. The end correction coefficients can
! be found by staggering, then summing, the Adams-Moulton coefficients
! at both ends.
wx=u1; wx(0)=o2; wx(nxh:nxh-ncor:-1)=cor
wy=u1; wy(0)=o2; wy(nyh:nyh-ncor:-1)=cor
wxy=outer_product(wx,wy); wxy=wxy/sum(wxy)
end subroutine get_wxy

!=============================================================================
subroutine get_qqw(nxh,nyh,ncor,j0xy,tw,p,q)!                          [get_qq]
!=============================================================================
! Like get_qqt, except the square norm involved in the definition of Q is
! modified by including a "trace-weight" proportion, tw, of the squared-trace.
! (Elsewhere tw is also known as "lambda".)
! In the elasticity analogue, this extra degree of freedom is like being
! able to include a nontrivial Poisson ratio defining the elastic modulus.
!=============================================================================
use pmat4, only: outer_product
use psym2, only: logsym2,expsym2
implicit none
integer(spi),                       intent(in   ):: nxh,nyh,ncor
real(dp),dimension(3,2,0:nxh,0:nyh),intent(in   ):: j0xy
real(dp),                           intent(in   ):: tw
real(dp),dimension(2,2),            intent(inout):: p
real(dp),                           intent(  out):: q
!-----------------------------------------------------------------------------
integer(spi),parameter         :: nit=5
real(dp),parameter             :: acrit=1.e-8_dp,dpx=.0099e0_dp
real(dp),dimension(0:nxh,0:nyh):: wxy
real(dp),dimension(3,2)        :: j0,j
real(dp),dimension(2,2)        :: el,pf,elp,elmean,g,ppx,pmx,ppy,pmy
real(dp),dimension(2)          :: hess,grad
real(dp)                       :: anorm,q00,qpx,qmx,qpy,qmy,c,w,twc
integer(spi)                   :: ix,iy,it
!=============================================================================
call get_wxy(nxh,nyh,ncor,wxy)! <- get 2D extended trapezoidal averaging wts
twc=u1-tw
if(p(1,1)==u0)then; p=u0; p(1,1)=u1; p(2,2)=u1; endif
! Iteratively calibrate preconditioner, p, to make elmean vanish:
anorm=u1
do it=1,nit
   elmean=u0
   q=u0
   do iy=0,nyh; do ix=0,nxh
      j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
! Precondition the Jacobian using latest iteration of P:
      j=matmul(j0,p)
! Find the Gram matrix, G, implied by the column vectors of the new J:
      g=matmul(transpose(j),j)
! Find the matrix logarithm, L = 0.5*log(G), contributions to elmean and q:
      call logsym2(g,el); el=el*o2; elmean=elmean+w*el
      q=q+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   enddo;      ; enddo
   if(anorm<acrit)exit ! <- Convergence criterion was met at last iteration
! Use double extended trapezoidal integration to find the domain-mean of L:
   elmean(1,2)=u0; elmean(2,1)=u0 ! <-Symmetrize by zeroing out off-diagonals
   elp=-elmean; call expsym2(elp,pf); p=matmul(p,pf) ! <- update P
   anorm=maxval(abs(elmean))
enddo
if(it>nit)&
print'("WARNING: In get_qqw, relatively inferior convergence; anorm=",1x,e12.5)',anorm
q00=q
ppx=p; ppx(1,1)=ppx(1,1)*(u1+dpx);qpx=u0
pmx=p; pmx(1,1)=pmx(1,1)*(u1-dpx);qmx=u0
ppy=p; ppy(2,2)=ppy(2,2)*(u1+dpx);qpy=u0
pmy=p; pmy(2,2)=pmy(2,2)*(u1-dpx);qmy=u0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,ppx); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el*o2;qpx=qpx+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,pmx); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el*o2;qmx=qmx+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,ppy); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el*o2;qpy=qpy+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,pmy); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el*o2;qmy=qmy+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
enddo; enddo
! Estimate a (diagonal) Hessian matrix and a gradient vector:
hess=(/ (qpx-u2*q00+qmx)/dpx**2, (qpy-u2*q00+qmy)/dpx**2 /)
hess=(/8._dp,8._dp/)
grad=(/ (qpx-qmx)/(u2*dpx)     , (qpy-qmy)/(u2*dpx)      /)
! If the hessian is positive, polish p with a final Newton iteration:
if(hess(1)>u0 .and. hess(2)>u0)then
   c=u1-grad(1)/hess(1); p(:,1)=p(:,1)*c
   c=u1-grad(2)/hess(2); p(:,2)=p(:,2)*c
endif

! and calculate the new q. Keep it only if it's numerically smaller than before:
q00=0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,p); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el*o2;q00=q00+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
enddo; enddo
if(q00<q)q=q00
end subroutine get_qqw

!=============================================================================
subroutine get_qqt(nxh,nyh,ncor,j0xy,p,q)!                           [get_qq]
!=============================================================================
! Assume the grid to be mirror-symmetric across both medians, so that the
! computation of the quality diagnostic, Q, need only involve the positive
! quadrant of the grid. The norm associated with the definition of Q is the
! Frobenius norm (Q is the grid-mean of the squared-Frobenius norm of the
! log of the Gram matrix of the given distribution of jacobian matrices.)
!=============================================================================
use pmat4, only: outer_product
use psym2, only: logsym2,expsym2
implicit none
integer(spi),                       intent(in   ):: nxh,nyh,ncor
real(dp),dimension(3,2,0:nxh,0:nyh),intent(in   ):: j0xy
real(dp),dimension(2,2),            intent(inout):: p
real(dp),                           intent(  out):: q
!-----------------------------------------------------------------------------
integer(spi),parameter         :: nit=7
real(dp),parameter             :: acrit=1.e-8_dp,dpx=.0099_dp
real(dp),dimension(0:nxh,0:nyh):: wxy
real(dp),dimension(3,2)        :: j0,j
real(dp),dimension(2,2)        :: el,pf,elp,elmean,g,ppx,pmx,ppy,pmy
real(dp),dimension(2)          :: hess,grad
real(dp)                       :: anorm,q00,qpx,qmx,qpy,qmy,c,w
integer(spi)                   :: ix,iy,it
!=============================================================================
call get_wxy(nxh,nyh,ncor,wxy)! <- get 2D extended trapezoidal averaging wts
if(p(1,1)==u0)then; p=u0; p(1,1)=u1; p(2,2)=u1; endif
! Iteratively calibrate preconditioner, p, to make elmean vanish:
anorm=1
do it=1,nit
   elmean=u0
   q=u0
   do iy=0,nyh; do ix=0,nxh
      j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
! Precondition the Jacobian using latest iteration of P:
      j=matmul(j0,p)
! Find the Gram matrix, G, implied by the column vectors of the new J:
      g=matmul(transpose(j),j)
! Find the matrix logarithm, L = log(G), contrinutions to elmean and q:
      call logsym2(g,el); el=el*o2; elmean=elmean+w*el; q=q+w*sum(el**2)
   enddo       ; enddo
   if(anorm<acrit)exit ! <- Convergence criterion was met at last iteration
! Use double extended trapezoidal integration to find the domain-mean of L:
   elmean(1,2)=u0; elmean(2,1)=u0 ! <-Symmetrize by zeroing out off-diganonals
   elp=-elmean; call expsym2(elp,pf); p=matmul(p,pf) ! <- update P
   anorm=maxval(abs(elmean))
enddo
if(it>nit)print'(" WARNING: In get_qqt, relatively inferior convergence")'

q00=q
ppx=p; ppx(1,1)=ppx(1,1)*(1+dpx);qpx=0
pmx=p; pmx(1,1)=pmx(1,1)*(1-dpx);qmx=0
ppy=p; ppy(2,2)=ppy(2,2)*(1+dpx);qpy=0
pmy=p; pmy(2,2)=pmy(2,2)*(1-dpx);qmy=0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,ppx); g=matmul(transpose(j),j)
   call logsym2(g,el); el=el/2; qpx=qpx+w*sum(el**2)
   j=matmul(j0,pmx); g=matmul(transpose(j),j)
   call logsym2(g,el); el=el/2; qmx=qmx+w*sum(el**2)
   j=matmul(j0,ppy); g=matmul(transpose(j),j)
   call logsym2(g,el); el=el/2; qpy=qpy+w*sum(el**2)
   j=matmul(j0,pmy); g=matmul(transpose(j),j)
   call logsym2(g,el); el=el/2; qmy=qmy+w*sum(el**2)
enddo; enddo
! Estimate a (diagonal) Hessian matrix and a gradient vector:
hess=(/ (qpx-2*q00+qmx)/dpx**2, (qpy-2*q00+qmy)/dpx**2 /)
grad=(/ (qpx-qmx)/(2*dpx)     , (qpy-qmy)/(2*dpx)      /)

! If the hessian is positive, polish the final p with a final Newton iteration:
if(hess(1)>0 .and. hess(2)>0.)then
   c=u1-grad(1)/hess(1); p(:,1)=p(:,1)*c
   c=u1-grad(2)/hess(2); p(:,2)=p(:,2)*c
endif

! and calculate the new q. Keep it only if is numerically smaller than before:
q00=0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,p); g=matmul(transpose(j),j)
   call logsym2(g,el); el=el*o2; q00=q00+w*sum(el**2)
enddo; enddo
if(q00<q)q=q00
end subroutine get_qqt

!============================================================================
subroutine get_qofmap(nxh,nyh,a,k,lam,xcedgex,xcedgey,  & !      [get_qofmap]
     q,xmedgex,xmedgey,ff)
!============================================================================
! For the map distortion formula with parameter, lam ("lambda"), and the
! Extended Schmidt-Gnomonic mapping with parameters, (A,K), find the quality
! grid-averaged distortion criterion, q, for domain in standardized
! polar orientation whose right-edge and top-edge midpoints are the
! cartesian unit 3-vectors, xcedgex, xcedgey. Work-space arrays are provided
! for the positive quadrant of the domain, gridded (from (0,0) at the
! projection center) out to nxh and nyh in x and y; these arrays are the
! xcgrid and xcdgrid. If their is an failure of the computation at any of
! the grid points, the failure flag, ff, is raised upon return.
!============================================================================
implicit none
integer(spi),                       intent(in ):: nxh,nyh
real(dp),                           intent(in ):: A,K,lam
real(dp),dimension(3),              intent(in ):: xcedgex,xcedgey
real(dp),                           intent(out):: q
real(dp),dimension(2),              intent(out):: xmedgex,xmedgey
logical,                            intent(out):: ff
!----------------------------------------------------------------------------
integer(spi),parameter             :: iend=2! <- trapezoidal end correction index.
real(dp),dimension(3)              :: xc
real(dp),dimension(3,2,0:nxh,0:nyh):: xcdgrid
real(dp),dimension(2)              :: xm
real(dp),dimension(2,2)            :: p
real(dp)                           :: dx,dy
integer(spi)                       :: ix,iy
logical                            :: ffgrid
!============================================================================
call xctoxm_ak(a,k,xcedgex,xmedgex,ff)! <- get map (x,y) of "right edge" midpt.
if(ff)then
   print'(" In get_qofmap, at A; failure flag upon return from xctoxm_ak")'
   return
endif
call xctoxm_ak(a,k,xcedgey,xmedgey,ff)! <- get map (x,y) of "top edge" midpt.
if(ff)then
   print'(" In get_qofmap, at B; failure flag upon return from xctoxm_ak")'
   stop
endif
! Set up a uniform nxh by nyh grid in the positive quadrant:
dx=xmedgex(1)/nxh
dy=xmedgey(2)/nyh
ff=F ! "false" is the default; "true" will mean something in xmtoxc failed.
! Map each uniform-grid point of the positive quadrant of the domain to
! the unit-sphere and record the absolute matrix-jacobian at each point:
do ix=0,nxh
   xm(1)=dx*ix
   do iy=0,nyh
      xm(2)=dy*iy
      call xmtoxc_ak(a,k,xm,xc,xcdgrid(:,:,ix,iy),ffgrid)
      if(ffgrid)ff=T ! <- set failure flag to .true. if it wasn't already
   enddo
enddo
if(ff)then
   print'(" In get_qofmap, at C, failure flag was raised in a xmtoxc_ac call")'
   return
endif
p=0
call get_qqw(nxh,nyh,iend,xcdgrid,lam,p,q)! <- Find the average quality criterion, Q.
end subroutine get_qofmap

!=============================================================================
subroutine get_bestesg(lam,arcx,arcy, a,k,m_arcx,m_arcy,q,ff)!   [get_betsesg]
!=============================================================================
! With prescribed optimization criterion parameter, lam ("lambda", typically 0.8),
! semi-arcs, arcx and arcy, of the domain in degrees, return the following:
! The best Extended Schmidt-Gnomonic (ESG) mapping parameters, A and K;
! The nondimensional map-space coordinates, m_arcx, m_arcy, of
! the optimally-mapped domain edges in the positive x and y directions;
! the optimality diagnostic of the domain, Q(lam), which this routine
! minimizes.
! If the process fails for any reason, the logical failure flag, FF, is raised
! (.true.), usually accompanied by a terse explanation, and the other output 
! arguments are then meaningless.
!=============================================================================
implicit none
real(dp),intent(in ):: lam,arcx,arcy
real(dp),intent(out):: a,k,m_arcx,m_arcy,q
logical, intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),    parameter       :: eps=1.e-9_dp,u1eps=u1+eps
real(dp),    parameter       :: darc=2._dp+eps ! <- arcx spacing of asplims
real(dp),    parameter       :: dlam=.2_dp+eps  ! <- lam spacing of asplims
integer(spi),parameter       :: larc=51,marc=55 ! <- arc index range of asplims
integer(spi),parameter       :: mlam=5          ! <- lam index range of asplims
real(dp),dimension(0:mlam,larc:marc):: asplims ! <- table of aspect ratio limits
real(dp)                     :: tarcx,tarcy,tm_arcx,tm_arcy,slam,sarcx, &
                                wi0,wi1,wj0,wj1,asplim,asp
integer(spi)                 :: i0,i1,j0,j1
logical                      :: flip
data asplims/&
        u1eps,   u1eps,   u1eps,   u1eps,   u1eps,   u1eps, & ! <- arcx=102
     0.999_dp,   u1eps,   u1eps,   u1eps,   u1eps,   u1eps, & ! <- arcx=104
     0.881_dp,0.882_dp,0.883_dp,0.885_dp,0.888_dp,0.888_dp, & ! <- arcx=106
     0.651_dp,0.653_dp,0.655_dp,0.658_dp,0.662_dp,0.662_dp, & ! <- arcx=108
     0.327_dp,0.329_dp,0.330_dp,0.332_dp,0.332_dp,0.333_dp/   ! <- arcx=110
!============================================================================
ff= lam<u0 .or. lam>=u1
if(ff)then
   print'("In get_bestesg; invalid optimization criterion parameter, lam")'
   return
endif
ff= arcx<=u0 .or. arcy<=u0
if(ff)then
   print'("In get_bestesg; a nonpositive domain parameter, arcx or arcy")'
   return
endif
! Make tarcx the longer of the two semi-axes of the domain:
flip=arcy>arcx
if(flip)then; tarcx=arcy; tarcy=arcx ! <- switch
else        ; tarcx=arcx; tarcy=arcy ! <- don't switch
endif
asp=tarcy/tarcx ! <- Domain aspect ratio that does not exceed one
sarcx=tarcx/darc
i0=floor(sarcx); i1=i0+1; wi0=i1-sarcx; wi1=sarcx-i0
ff=i1>marc
if(ff)then
   print'("In get_bestesg; domain length too large")'; return
endif
asplim=u1 ! <- Default aspect ratio limit for small tarcx 
if(i0>=larc)then ! Interpolate aspect limit for this large tarcx from the table:
   slam=lam/dlam
   j0=floor(slam) ; j1=j0+1; wj0=j1-slam;  wj1=slam -j0
   asplim=(asplims(j0,i0)*wi0+asplims(j0,i1)*wi1)*wj0 + &
          (asplims(j1,i0)*wi0+asplims(j1,i1)*wi1)*wj1
endif
ff=asp>asplim
if(ff)then
   print'("In get_bestesg; domain width too large for given domain length")'
   return
endif
call get_bestesgt(lam,asp,tarcx, A,K,tm_arcx,tm_arcy,q,ff)
if(ff)return
if(flip)then; m_arcx=tm_arcy; m_arcy=tm_arcx
else        ; m_arcx=tm_arcx; m_arcy=tm_arcy
endif
end subroutine get_bestesg

!=============================================================================
subroutine get_bestesgt(lam,asp,arcx, a,k,m_arcx,m_arcy,q,ff)!   [get_betsesgt]
!=============================================================================
! With prescribed optimization criterion parameter, lam ("lambda"),
! aspect ratio, asp (0.1 < asp <= 1.), major semi-arc, in degrees, of arcx,
! return the best Extended Schmidt-Gnomonic (ESG) mapping parameters, A and K,
! for the rectangular domain whose edge midpoints are distant arcx and asp*arcx
! from its center.
! Also, return the map-space x and y coordinates, m_arcx, m_arcy, of
! the domain edges in the positive directions.
! The optimality criterion, Q(A,K), which this routine aims to minimize, 
! is provided upon return. However A and K are independenly scaled, the contours
! Q are highly elliptical; the differencing stencil is adjusted to mimic this
! stretching to preserve good numerical conditioning of the calculations and, we
! hope, tl thereby ensure an accurate and robust estimation of the location of
! the minimum.
! If the process fails for any reason (such as parameter combinations, lam, asp,
! arcx for which no minimum of Q in the proper interior of the valid space
! of these parameters exists), then the failure flag, FF, is raised (.true.)
! and the output parameters are then of course meaningless. 
!=============================================================================
use pietc, only: u5,o5,s18,s36,s54,s72,ms18,ms36,ms54,ms72,pi2
use psym2, only: chol2
implicit none
real(dp),intent(in ):: lam,arcx,asp
real(dp),intent(out):: a,k,m_arcx,m_arcy,q
logical, intent(out):: ff
!-----------------------------------------------------------------------------
integer(spi),parameter         :: narc=11,nasp=10,nit=8,nang=5
real(dp),parameter             :: eps=1.e-7_dp,u2o5=u2*o5,darc=10._dp+eps,&
                                  dasp=.1_dp+eps,dang=pi2*o5,r=0.0001_dp,rr=r*r, &
                                  urc=.4_dp, & ! <- Under-relaxation coefficient 
                                  enxyq=10000._dp! ~ pts in trial grid quadrant
real(dp),parameter             :: f18=u2o5*s18,f36=u2o5*s36,f54=u2o5*s54,f72=u2o5*s72,&
                                  mf18=-f18,mf36=-f36,mf54=-f54,mf72=-f72 !<- (Fourier)
real(dp),dimension(0:4,0:4)    :: em5 ! <- Fourier matrix for 5 points
real(dp),dimension(nasp,0:narc):: adarc,kdarc ! < 1st guess tables of A and K.
real(dp),dimension(0:nang-1)   :: qs
real(dp),dimension(2,0:nang-1) :: aks
real(dp),dimension(2)          :: grad,v2,ak
real(dp),dimension(2,2)        :: hess,el,basis0,basis
real(dp),dimension(3)          :: xcedgex,xcedgey
real(dp),dimension(2)          :: xmedgex,xmedgey
real(dp)                       :: ang,sarcx,sasp,wx0,wx1,wa0,wa1,qold
integer(spi)                   :: i,it,iarcx0,iasp0,iarcx1,iasp1,&
                                  nxh,nyh
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

data em5/o5, u2o5,   u0, u2o5,   u0, & ! <- The Fourier matrix for 5 points. Applied to
         o5,  f18,  f72, mf54,  f36, & ! the five 72-degree spaced values in a column-
         o5, mf54,  f36,  f18, mf72, & ! vector, the product vector has the components,
         o5, mf54, mf36,  f18,  f72, & ! wave-0, cos and sin wave-1, cos and sin wave-2.
         o5,  f18, mf72, mf54, mf36/

data basis0/0.1_dp,u0,  0.3_dp,0.3_dp/! <- initial basis for orienting (a,k) differencing
!============================================================================
sasp=asp/dasp
iasp0=floor(sasp); wa1=sasp-iasp0
iasp1=iasp0+1;     wa0=iasp1-sasp
sarcx=arcx/darc
iarcx0=floor(sarcx); wx1=sarcx-iarcx0
iarcx1=iarcx0+1;     wx0=iarcx1-sarcx
if(iasp0 <1 .or. iasp1 >nasp)stop 'Aspect ratio out of range'
if(iarcx0<0 .or. iarcx1>narc)stop 'Major semi-arc is out of range'
nxh=nint(sqrt(enxyq/asp))! < These nxh and nyh ensure nearly
nyh=nint(sqrt(enxyq*asp))! square trial grid cells, and about enxyq points in total.
call getedges(arcx,asp*arcx, xcedgex,xcedgey)

! Bilinearly interpolate A and K from crude table into a 2-vector:
ak=(/wx0*(wa0*adarc(iasp0,iarcx0)+wa1*adarc(iasp1,iarcx0))+ &
     wx1*(wa0*adarc(iasp0,iarcx1)+wa1*adarc(iasp1,iarcx1)), &
     wx0*(wa0*kdarc(iasp0,iarcx0)+wa1*kdarc(iasp1,iarcx0))+ &
     wx1*(wa0*kdarc(iasp0,iarcx1)+wa1*kdarc(iasp1,iarcx1))/)
basis=basis0 ! <- initial the basis to a representative guess.
qold=100._dp ! <- initialize qold to a meaningless large value to force at least one
             !    complete Newton iteration to occur.
do it=1,nit
   call get_qofmap(nxh,nyh,ak(1),ak(2),lam,xcedgex,xcedgey,&
        q,xmedgex,xmedgey,ff)
   if(ff)return
   if(q>=qold)exit ! <-Assume this condition indicates early convergence
   m_arcx=xmedgex(1)
   m_arcy=xmedgey(2)
   qold=q

! Place five additional sample points around the stencil-ellipse:
   do i=0,4
      ang=i*dang ! steps of 72 degrees 
      v2=(/cos(ang),sin(ang)/)*r ! points on a circle of radius r ...
      aks(:,i)=ak+matmul(basis,v2) ! ... become points on an ellipse.
! Get quality, qs(i)
      call get_qofmap(nxh,nyh,aks(1,i),aks(2,i),lam,xcedgex,xcedgey,&
           qs(i),xmedgex,xmedgey,ff)
   enddo
   if(ff)return
! Recover gradient and hessian estimates, normalized by q, from all 
! 6 samples, q, qs. These are wrt the basis, not (a,k) directly.
! Make qs the 5-pt discrete Fourier coefficients of the ellipse pts:
   qs=matmul(em5,qs)/q
   grad=qs(1:2)/r ! <- r is the finite differencing step size parameter
   qs(0)=qs(0)-u1 ! <- cos(2*ang) coefficient relative to the central value.
   hess(1,1)=qs(0)+qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
   hess(1,2)=qs(4)      ! <- sin(2*ang) coefficient
   hess(2,1)=qs(4)      !
   hess(2,2)=qs(0)-qs(3)! <- combine cos(0*ang) and cos(2*ang) coefficients
   hess=hess*u2/rr ! <- rr is r**2

! Perform a Cholesky decomposition of the hessian:
   call chol2(hess,el,ff)
   if(ff)then
      print'(" In get_bestESG, hessian is not positive; cholesky fails")'
      return
   endif
! Invert the cholesky factor in place:
   el(1,1)=u1/el(1,1)
   el(2,2)=u1/el(2,2)
   el(2,1)=-el(2,1)*el(1,1)*el(2,2)
! Estimate a Newton step towards the minimum of function Q(a,k):
   v2=-matmul(transpose(el),matmul(el,grad))
! qt=q+dot_product(grad,v2)*o2 ! <- Estimates what the new minimum will be
! Apply an under-relaxation in the first iteration:
   if(it<=1)v2=v2*urc
   ak=ak+matmul(basis,v2)! <- increment the parameter vector estimate
! Use the inverse cholesky factor to re-condition the basis. This is to make
! the next stencil-ellipse more closely share the shape of the elliptical
! contours of Q near its minumum -- essentially a preconditioning of the
! numerical optimization:
   basis=matmul(basis,transpose(el))
enddo
a=ak(1); k=ak(2)

end subroutine get_bestesgt

!=============================================================================
subroutine get_bestesg_inv(lam,m_arcx,m_arcy, a,k,arcx,arcy, q,ff)![get_bestesg_inv]
!=============================================================================
! A form of inverse of get_bestesg where the desired map-space edge
! coordinates, m_arcx and m_arcy, are input throught the argument
! list, and the A, K, arcx, and arcy (degrees) are output such that,
! if get_bestesg are called with these arcx and arcy as inputs, then
! m_arcx and m_arcx that would be returned are exactly the desired
! values. The A and K returned here are consistent also.
!=============================================================================
implicit none
real(dp),intent(in ):: lam,m_arcx,m_arcy
real(dp),intent(out):: a,k,arcx,arcy,q
logical, intent(out):: ff
!-----------------------------------------------------------------------------
integer(spi),parameter:: nit=40         ! <- Number of iterations allowed
real(dp),    parameter:: crit=1.e-12_dp ! <- Convergence criterion
real(dp)              :: rx,ry,m_arcxa,m_arcya
integer(spi)          :: it
!=============================================================================
arcx=m_arcx*rtod; arcy=m_arcy*rtod ! < 1st guess (degrees) of arcx and arcy
do it=1,nit
   call get_bestesg(lam,arcx,arcy, a,k,m_arcxa,m_arcya, q,ff)
   if(ff)then
      print'("In get_bestesg_inv; raised failure flag prevents completion")'
      return
   endif
   rx=m_arcx/m_arcxa; ry=m_arcy/m_arcya
   arcx=arcx*rx     ; arcy=arcy*ry
   if(abs(rx-u1)<=crit .and. abs(ry-u1)<=crit)return
enddo
print'("WARNING; in get_bestesg_inv;")'
print'("full convergence unattained after",i3," iterations")',nit
print'("Residual proportionate mismatch of arcx and of arcy:",2(1x,e12.5))',&
     rx-u1,ry-u1
end subroutine get_bestesg_inv

!=============================================================================
subroutine hgrid_ak_rr(lx,ly,nx,ny,A,K,plat,plon,pazi, & !       [hgrid_ak_rr]
     delx,dely,  glat,glon,garea, ff)
!=============================================================================
! Use a and k as the parameters of an Extended Schmidt-transformed
! Gnomonic (ESG) mapping centered at (plat,plon) and twisted about this center
! by an azimuth angle of pazi counterclockwise (these angles in radians).
!
! Assume the radius of the earth is unity, and using the central mapping
! point as the coordinate origin, set up the grid with central x-spacing delx
! and y-spacing dely. The grid index location of the left-lower
! corner of the domain is (lx,ly) (typically both NEGATIVE).
! The numbers of the grid spaces in x and y directions are nx and ny.
! (Note that, for a centered rectangular grid lx and ly are negative and, in
! magnitude, half the values of nx and ny respectively.)
! Return the latitude and longitude, in radians again, of the grid points thus
! defined in the arrays, glat and glon, and return a rectangular array, garea,
! of dimensions nx-1 by ny-1, that contains the areas of each of the grid cells
!
! if all goes well, return a lowered failure flag, ff=.false. . But if, for some
! reason, it is not possible to complete this task, return the raised failure
! flag, ff=.TRUE. .
!=============================================================================
use pmat4, only: sarea
use pmat5, only: ctogr
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(3,2):: xcd
real(dp),dimension(3)  :: xc
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,&
                          rlat,drlata,drlatb,drlatc,    &
                          rlon,drlona,drlonb,drlonc
integer(spi)           :: ix,iy,mx,my
!=============================================================================
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

!=============================================================================
subroutine hgrid_ak_rr_c(lx,ly,nx,ny,a,k,plat,plon,pazi, & !     [hgrid_ak_rr]
                    delx,dely,  glat,glon,garea,dx,dy,angle_dx,angle_dy, ff)
!=============================================================================
! Use a and k as the parameters of an extended Schmidt-transformed
! gnomonic (ESG) mapping centered at (plat,plon) and twisted about this center
! by an azimuth angle of pazi counterclockwise (these angles in radians).
!
! Using the central mapping point as the coordinate origin, set up the grid
! with central x-spacing delx and y-spacing dely in nondimensional units, (i.e.,
! as if the earth had unit radius) and with the location of the left-lower
! corner of the grid at center-relative grid index pair, (lx,ly) and with
! the number of the grid spaces in x and y directions given by nx and ny.
! (Note that, for a centered rectangular grid lx and ly are negative and, in
! magnitude, half the values of nx and ny respectively.)
! Return the latitude and longitude, again, in radians, of the grid points thus
! defined in the arrays, glat and glon; return a rectangular array, garea,
! of dimensions nx-1 by ny-1, that contains the areas of each of the grid cells
! in nondimensional "steradian" units.
! In this version, these grid cell areas are computed by 2D integrating the
! scalar jacobian of the transformation, using a 4th-order centered scheme.
! The estimated grid steps, dx and dy, are returned at the grid cell edges,
! using the same 4th-order scheme to integrate the 1D projected jacobian.
! The angles, relative to local east and north, are returned respectively
! as angle_dx and angle_dy at grid cell corners, in radians counterclockwise.
!
! if all goes well, return a .FALSE. failure flag, ff. If, for some
! reason, it is not possible to complete this task, return the failure flag
! as .TRUE.
!=============================================================================
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
!-----------------------------------------------------------------------------
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
!=============================================================================
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

!=============================================================================
subroutine hgrid_ak_rc(lx,ly,nx,ny,A,K,plat,plon,pazi, & !       [hgrid_ak_rc]
     delx,dely, xc,xcd,garea, ff)
!=============================================================================
! Use a and k as the parameters of an Extended Schmidt-transformed
! Gnomonic (ESG) mapping centered at (plat,plon) and twisted about this center
! by an azimuth angle of pazi counterclockwise (these angles in radians).
!
! Assume the radius of the earth is unity, and using the central mapping
! point as the coordinate origin, set up the grid with central x-spacing delx
! and y-spacing dely. The grid index location of the left-lower
! corner of the domain is (lx,ly) (typically both NEGATIVE).
! The numbers of the grid spaces in x and y directions are nx and ny.
! (Note that, for a centered rectangular grid lx and ly are negative and, in
! magnitude, half the values of nx and ny respectively.)
! Return the unit cartesian vectors xc of the grid points and their jacobian
! matrices xcd wrt the map coordinates, and return a rectangular array, garea,
! of dimensions nx-1 by ny-1, that contains the areas of each of the grid cells
!
! if all goes well, return a lowered failure flag, ff=.false. . But if, for some
! reason, it is not possible to complete this task, return the raised failure
! flag, ff=.TRUE. .
!=============================================================================
use pmat4, only: sarea
use pmat5, only: ctogr
implicit none
integer(spi),intent(in ):: lx,ly,nx,ny
real(dp),    intent(in ):: a,k,plat,plon,pazi,delx,dely
real(dp),dimension(3,  lx:lx+nx  ,ly:ly+ny  ),intent(out):: xc
real(dp),dimension(3,2,lx:lx+nx  ,ly:ly+ny  ),intent(out):: xcd
real(dp),dimension(    lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                      intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,                  &
                          rlat,rlata,rlatb,rlatc,drlata,drlatb,drlatc,    &
                          rlon,rlona,rlonb,rlonc,drlona,drlonb,drlonc
integer(spi)           :: ix,iy,mx,my
!=============================================================================
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

!=============================================================================
subroutine hgrid_ak_dd(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, & !    [hgrid_ak_dd]
     delx,dely,  gdlat,gdlon,garea, ff)
!=============================================================================
! Use a and k as the parameters of an Extended Schmidt-transformed
! Gnomonic (ESG) mapping centered at (pdlat,pdlon) and twisted about this center
! by an azimuth angle of pdazi counterclockwise (these angles in degrees).
! Like hgrid_ak_rr, return the grid points' lats and lons, except that here
! the angles are returned in degrees. Garea, the area of each grid cell, is
! returned as in hgrid_ak_rr, and a failure flag, ff, raised when a failure
! occurs anywhere in these calculations
!============================================================================
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: A,K,pdlat,pdlon,pdazi,&
                                                        delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: gdlat,gdlon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: plat,plon,pazi
!=============================================================================
plat=pdlat*dtor ! Convert these angles from degrees to radians
plon=pdlon*dtor !
pazi=pdazi*dtor !
call hgrid_ak_rr(lx,ly,nx,ny,A,K,plat,plon,pazi, &
    delx,dely,   gdlat,gdlon,garea, ff)
if(ff)return
gdlat=gdlat*rtod ! Convert these angles from radians to degrees
gdlon=gdlon*rtod !
end subroutine hgrid_ak_dd
!=============================================================================
subroutine hgrid_ak_dd_c(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, &!   [hgrid_ak_dd]
     delx,dely,  gdlat,gdlon,garea,dx,dy,dangle_dx,dangle_dy, ff)
!=============================================================================
! Like hgrid_ak_rr_c, except all the angle arguments (but not delx,dely)
! are in degrees instead of radians.
!=============================================================================
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,pdlat,pdlon,pdazi, &
                                                        delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: gdlat,gdlon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
real(dp),dimension(lx:lx+nx-1,ly:ly+ny  ),intent(out):: dx
real(dp),dimension(lx:lx+nx  ,ly:ly+ny-1),intent(out):: dy
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: dangle_dx,dangle_dy
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: plat,plon,pazi
!=============================================================================
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

!=============================================================================
subroutine hgrid_ak_dc(lx,ly,nx,ny,a,k,pdlat,pdlon,pdazi, & !    [hgrid_ak_dc]
     delx,dely, xc,xcd,garea, ff)
!=============================================================================
! Use a and k as the parameters of an Extended Schmidt-transformed
! Gnomonic (ESG) mapping centered at (pdlat,pdlon) and twisted about this center
! by an azimuth angle of pdazi counterclockwise (these angles in degrees).
! Like hgrid_ak_rx, return the grid points' cartesians xc and jacobian matrices,
! xcd. Garea, the area of each grid cell, is also
! returned as in hgrid_ak_rx, and a failure flag, ff, raised when a failure
! occurs anywhere in these calculations
!============================================================================
implicit none
integer(spi),intent(in ):: lx,ly,nx,ny
real(dp),    intent(in ):: A,K,pdlat,pdlon,pdazi,delx,dely
real(dp),dimension(3,  lx:lx+nx  ,ly:ly+ny  ),intent(out):: xc
real(dp),dimension(3,2,lx:lx+nx  ,ly:ly+ny  ),intent(out):: xcd
real(dp),dimension(    lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                      intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: plat,plon,pazi
!=============================================================================
plat=pdlat*dtor
plon=pdlon*dtor
pazi=pdazi*dtor
call hgrid_ak_rc(lx,ly,nx,ny,A,K,plat,plon,pazi, &
    delx,dely,   xc,xcd,garea, ff)
end subroutine hgrid_ak_dc

!=============================================================================
subroutine hgrid_ak(lx,ly,nx,ny,a,k,plat,plon,pazi, & !             [hgrid_ak]
     re,delxre,delyre,  glat,glon,garea, ff)
!=============================================================================
! Like hgrid_ak_rr_c except the argument list includes the earth radius, re,
! and this is used to express the map-space grid increments in the dimensional
! units, delxre, delyre on entry, and to express the grid cell areas, garea,
! in dimensional units upon return.
! The gridded lats and lons, glat and glon, remain in radians.
!============================================================================
implicit none
integer(spi),                             intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        re,delxre,delyre
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: delx,dely,rere
!=============================================================================
delx=delxre/re ! <- set nondimensional grid step delx
dely=delyre/re ! <- set nondimensional grid step dely
call hgrid_ak_rr(lx,ly,nx,ny,a,k,plat,plon,pazi, &
     delx,dely,  glat,glon,garea, ff)
if(ff)return
rere=re*re
garea=garea*rere ! <- Convert from steradians to physical area units.
end subroutine hgrid_ak

!=============================================================================
subroutine hgrid_ak_c(lx,ly,nx,ny,a,k,plat,plon,pazi, & !           [hgrid_ak]
     re,delxre,delyre,  glat,glon,garea,dx,dy,dangle_dx,dangle_dy, ff)
!=============================================================================
! Like hgrid_ak_rr_c except the argument list includes the earth radius, re,
! and this is used to express the map-space grid increments in the dimensional
! units, delxre, delyre on entry, and to express the grid cell areas, garea,
! and the x- and y- grid steps, dx and dy, in dimensional units upon return.
! The gridded lats and lons, glat and glon, remain in radians.
! Also, in order for the argument list to remain compatible with an earlier
! version of this routine, the relative rotations of the steps, dangle_dx
! and dangle_dy, are returned as degrees instead of radians (all other angles
! in the argument list, i.e., plat,plon,pazi,glat,glon, remain radians).
!============================================================================
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
!-----------------------------------------------------------------------------
real(dp):: delx,dely,rere
!=============================================================================
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

end module pesg

