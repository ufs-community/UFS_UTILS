!=============================================================================
subroutine get_qqt(nxh,nyh,ncor,j0xy,p,q)
!=============================================================================
! Assume the grid to be mirror-symmetric across both medians, so that the
! computation of the quality diagnostic, Q, need only involve the positive
! quadrant of the grid. The norm associated with the definition of Q is the
! Frobenius norm (Q is the grid-mean of the squared-Frobenius norm of the
! log of the Gram matrix of the given distribution of jacobian matrices.)
!=============================================================================
use pkind, only: dp
use pietc, only: u0,u1,o2
use pmat4, only: outer_product
use psym2
implicit none
integer,                            intent(in   ):: nxh,nyh,ncor
real(dp),dimension(3,2,0:nxh,0:nyh),intent(in   ):: j0xy
real(dp),dimension(2,2),            intent(inout):: p
real(dp),                           intent(  out):: q
!-----------------------------------------------------------------------------
integer,parameter                  :: nit=5
real(dp),parameter                 :: acrit=1.e-8,dpx=.0099
real(dp),dimension(0:nxh,0:nyh)    :: wxy
real(dp),dimension(3,2)            :: j0,j
real(dp),dimension(2,2)            :: el,pf,elp,elmean,g,ppx,pmx,ppy,pmy
real(dp),dimension(2)              :: hess,grad
real(dp)                           :: anorm,q00,qpx,qmx,qpy,qmy,c,w
integer                            :: ix,iy,lx,ly,it
!=============================================================================
call get_wxy(nxh,nyh,ncor,wxy)! <- get 2D extended trapezoidal averaging wts
if(p(1,1)==u0)then; p=0; p(1,1)=u1; p(2,2)=u1; endif
! Iteratively calibrate preconditioner, p, to make elmean vanish:
anorm=1
do it=1,nit
   elmean=0
   q=0
   do iy=0,nyh; do ix=0,nxh
      j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
! Precondition the Jacobian using latest iteration of P:
      j=matmul(j0,p)
! Find the Gram matrix, G, implied by the column vectors of the new J:
      g=matmul(transpose(j),j)
! Find the matrix logarithm, L = log(G), contrinutions to elmean and q:
      call logsym2(g,el); el=el/2; elmean=elmean+w*el; q=q+w*sum(el**2)
   enddo       ; enddo
   if(anorm<acrit)exit ! <- Convergence criterion was met at last iteration
! Use double extended trapezoidal integration to find the domain-mean of L:   
   elmean(1,2)=0; elmean(2,1)=0 ! <-Symmetrize by zeroing out off-diganonals
   elp=-elmean; call expsym2(elp,pf); p=matmul(p,pf) ! <- update P
   anorm=maxval(abs(elmean))
enddo
if(it>nit)then
   print'("WARNING: In get_qqt, apparent failure of iteration to converge")'
   read(*,*)
endif

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

!!print'('' hessian components:'',t30,2(1x,e20.14))',hess !!!!!!!
!!print'('' grad    components:'',t30,2(1x,e20.14))',grad !!!!!!!
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
   call logsym2(g,el); el=el/2; q00=q00+w*sum(el**2)
enddo; enddo
!!print'('' adjusted final q: '',e20.14)',q00
if(q00<q)q=q00
end subroutine get_qqt

!=============================================================================
subroutine get_qqw(nxh,nyh,ncor,j0xy,tw,p,q)
!=============================================================================
! Like get_qqt, except the square norm involved in the definition of Q is
! modified by including a "trace-weight" proportion, tw, of the squared-trace.
! (Elsewhere tw is also known as "lambda".)
! In the elasticity analogue, this extra degree of freedom is like being
! able to include a nontrivial Poisson ratio defining the elastic modulus.
!============================================================================= 
use pkind, only: dp
use pietc, only: u0,u1,o2
use pmat4, only: outer_product
use psym2
implicit none
integer,                            intent(in   ):: nxh,nyh,ncor
real(dp),dimension(3,2,0:nxh,0:nyh),intent(in   ):: j0xy
real(dp),                           intent(in   ):: tw
real(dp),dimension(2,2),            intent(inout):: p
real(dp),                           intent(  out):: q
!-----------------------------------------------------------------------------
integer,parameter              :: nit=5
real(dp),parameter             :: acrit=1.e-8,dpx=.0099
real(dp),dimension(0:nxh,0:nyh):: wxy
real(dp),dimension(3,2)        :: j0,j
real(dp),dimension(2,2)        :: el,pf,elp,elmean,g,ppx,pmx,ppy,pmy
real(dp),dimension(2)          :: hess,grad
real(dp)                       :: anorm,q00,qpx,qmx,qpy,qmy,c,w,twc
integer                        :: ix,iy,lx,ly,it
!=============================================================================
call get_wxy(nxh,nyh,ncor,wxy)! <- get 2D extended trapezoidal averaging wts
twc=u1-tw
if(p(1,1)==u0)then; p=0; p(1,1)=u1; p(2,2)=u1; endif
! Iteratively calibrate preconditioner, p, to make elmean vanish:
anorm=1
do it=1,nit
   elmean=0
   q=0
   do iy=0,nyh; do ix=0,nxh
      j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
! Precondition the Jacobian using latest iteration of P:
      j=matmul(j0,p)
! Find the Gram matrix, G, implied by the column vectors of the new J:
      g=matmul(transpose(j),j)
! Find the matrix logarithm, L = 0.5*log(G), contributions to elmean and q:
      call logsym2(g,el); el=el/2; elmean=elmean+w*el
      q=q+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   enddo;      ; enddo
   if(anorm<acrit)exit ! <- Convergence criterion was met at last iteration
! Use double extended trapezoidal integration to find the domain-mean of L:   
   elmean(1,2)=0; elmean(2,1)=0 ! <-Symmetrize by zeroing out off-diagonals
   elp=-elmean; call expsym2(elp,pf); p=matmul(p,pf) ! <- update P
   anorm=maxval(abs(elmean))
enddo
if(it>nit)then
   print'("WARNING: In get_qqt, apparent failure of iteration to converge")'
   read(*,*)
endif

q00=q
ppx=p; ppx(1,1)=ppx(1,1)*(1+dpx);qpx=0
pmx=p; pmx(1,1)=pmx(1,1)*(1-dpx);qmx=0
ppy=p; ppy(2,2)=ppy(2,2)*(1+dpx);qpy=0
pmy=p; pmy(2,2)=pmy(2,2)*(1-dpx);qmy=0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,ppx); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el/2;qpx=qpx+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,pmx); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el/2;qmx=qmx+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,ppy); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el/2;qpy=qpy+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
   j=matmul(j0,pmy); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el/2;qmy=qmy+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
enddo; enddo
! Estimate a (diagonal) Hessian matrix and a gradient vector:
hess=(/ (qpx-2*q00+qmx)/dpx**2, (qpy-2*q00+qmy)/dpx**2 /)
hess=(/8.,8./)
grad=(/ (qpx-qmx)/(2*dpx)     , (qpy-qmy)/(2*dpx)      /)
! If the hessian is positive, polish p with a final Newton iteration: 
if(hess(1)>0 .and. hess(2)>0.)then
   c=u1-grad(1)/hess(1); p(:,1)=p(:,1)*c
   c=u1-grad(2)/hess(2); p(:,2)=p(:,2)*c
endif

! and calculate the new q. Keep it only if it's numerically smaller than before:
q00=0
do iy=0,nyh; do ix=0,nxh
   j0=j0xy(:,:,ix,iy); w=wxy(ix,iy)
   j=matmul(j0,p); g=matmul(transpose(j),j)
   call logsym2(g,el);el=el/2;q00=q00+w*(twc*sum(el**2)+tw*(el(1,1)+el(2,2))**2)
enddo; enddo
if(q00<q)q=q00
end subroutine get_qqw

!=============================================================================
subroutine get_wxy(nxh,nyh,ncor,wxy)
!=============================================================================
use pkind, only: dp
use pietc, only: o2,u1
use pmat4, only: outer_product
implicit none
integer,                        intent(in ):: nxh,nyh,ncor
real(dp),dimension(0:nxh,0:nyh),intent(out):: wxy
!-----------------------------------------------------------------------------
integer,parameter         :: dencor0=2,dencor1=12,dencor2=24, &
                             dencor3=720,dencor4=1440! <-denominators
real(dp),dimension(0:nxh) :: wx
real(dp),dimension(0:nyh) :: wy
real(dp),dimension(0:ncor):: cor ! Becomes the full end-correction vector
integer,dimension(0:0)    :: numcor0 ! numerator for uncorrected trapezoidal
integer,dimension(0:1)    :: numcor1 ! numerators for common extended trap.
integer,dimension(0:2)    :: numcor2 ! ..numerators for higher order schemes ..
integer,dimension(0:3)    :: numcor3 !
integer,dimension(0:4)    :: numcor4 ! 
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
subroutine getedges(arcx,arcy,edgex,edgey)
!=============================================================================
! For angles (degrees) of the arcs spanning the halfwidths between the
! region's center and its x and y edges, get the two cartesian vectors
! that represent the locations of these edge midpoints in the positive x and y
! directions.
!=============================================================================
use pkind, only: dp
use pietc, only: u0,dtor
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

!==============================================================================
subroutine xmtoxc_ak(a,kappa,xm,xc,xcd,ff)
!==============================================================================
! Assuming the A-kappa parameterization of the generalized schmidt-transformed
! gnomonic mapping, and given a map-space 2-vector, xm, find the corresponding
! cartesian unit 3-vector and its derivative wrt xm, jacobian matrix, xcd.
! If for any reason the mapping cannot be done, return a raised failure
! flag, FF.
!=============================================================================
use pkind, only: dp
use pietc, only: T,F
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
subroutine xctoxm_ak(a,kappa,xc,xm,ff)
!=============================================================================
! Inverse mapping of xmtoxc_ak. That is, go from given cartesian unit 3-vector,
! xc, to map coordinate 2-vector xm (or return a raised failure flag, FF, if
! the attempt fails).
!=============================================================================
use pkind, only: dp
use pietc, only: F,T,u0,u1
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
subroutine zmtozt(a,zm,zt,ztd,ff)
!==============================================================================
! Evaluate the function, zt = tan(sqrt(A)*z)/sqrt(A), and its deivative, ztd,
! for positive and negative A and for the limiting case, A --> 0
!==============================================================================
use pkind, only: dp
use pietc, only: F,T,u1,pih
implicit none
real(dp),intent(in ):: a,zm
real(dp),intent(out):: zt,ztd
logical, intent(out):: ff
!------------------------------------------------------------------------------
real(dp):: ra
!==============================================================================
ff=f
if    (a>0)then; ra=sqrt( a); zt=tan (ra*zm)/ra; ff=abs(ra*zm)>=pih
elseif(a<0)then; ra=sqrt(-a); zt=tanh(ra*zm)/ra
else                        ; zt=zm
endif
ztd=u1+a*zt*zt
end subroutine zmtozt

!=============================================================================
subroutine zttozm(a,zt,zm,ff)
!=============================================================================
! Inverse of zmtozt
!=============================================================================
use pkind, only: dp
use pietc, only: F,u1
implicit none
real(dp),intent(in ):: a,zt
real(dp),intent(out):: zm
logical, intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: ra,razt
!=============================================================================
ff=F
if    (a>0)then; ra=sqrt( a); razt=ra*zt; zm=atan (razt)/ra
elseif(a<0)then; ra=sqrt(-a); razt=ra*zt; ff=abs(razt)>=u1; if(ff)return
                                          zm=atanh(razt)/ra
else                                    ; zm=zt
endif
end subroutine zttozm

!==============================================================================
subroutine xmtoxt(a,xm,xt,xtd,ff)
!==============================================================================
! Like zmtozt, but for 2-vector xm and xt, and 2*2 diagonal Jacobian xtd
!==============================================================================
use pkind, only: dp
implicit none
real(dp),               intent(in ):: a
real(dp),dimension(2),  intent(in ):: xm
real(dp),dimension(2),  intent(out):: xt
real(dp),dimension(2,2),intent(out):: xtd
logical,                intent(out):: ff
!-----------------------------------------------------------------------------
integer:: i
!==============================================================================
xtd=0; do i=1,2; call zmtozt(a,xm(i),xt(i),xtd(i,i),ff); if(ff)return; enddo
end subroutine xmtoxt

!=============================================================================
subroutine xttoxm(a,xt,xm,ff)
!=============================================================================
! Inverse of xmtoxt
!============================================================================
use pkind, only: dp
use pietc, only: F
implicit none
real(dp),             intent(in ):: a
real(dp),dimension(2),intent(in ):: xt
real(dp),dimension(2),intent(out):: xm
logical              ,intent(out):: ff
!-----------------------------------------------------------------------------
integer:: i
!=============================================================================
do i=1,2; call zttozm(a,xt(i),xm(i),ff); if(ff)return; enddo
end subroutine xttoxm

!==============================================================================
subroutine xttoxs(kappa,xt,xs,xsd,ff)
!==============================================================================
! Scaled gnomonic plane xt to standard stereographic plane xs
!==============================================================================
use pkind, only: dp
use pietc, only: u0,u1
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
subroutine xstoxt(kappa,xs,xt,ff)
!=============================================================================
! Inverse of xttoxs.
!=============================================================================
use pkind, only: dp
use pietc, only: u1
implicit none
real(dp),             intent(in ):: kappa
real(dp),dimension(2),intent(in ):: xs
real(dp),dimension(2),intent(out):: xt
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(dp):: s,sc
!=============================================================================
s=kappa*(xs(1)*xs(1)+xs(2)*xs(2)); sc=u1-s
ff=(sc<=0); if(ff)return
xt=2*xs/sc
end subroutine xstoxt

!=============================================================================
subroutine xstoxc(xs,xc,xcd)
!=============================================================================
! Standard transformation from polar stereographic map coordinates, xs, to
! cartesian, xc, assuming the projection axis is polar.
! xcd=d(xc)/d(xs) is the jacobian matrix
!=============================================================================
use pkind, only: dp
use pietc, only: u1,u2
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
subroutine xctoxs(xc,xs)
!=============================================================================
! Inverse of xstoxc. I.e., cartesians to stereographic
!=============================================================================
use pkind, only: dp
use pietc, only: u1
implicit none
real(dp),dimension(3),intent(in ):: xc
real(dp),dimension(2),intent(out):: xs
!-----------------------------------------------------------------------------
real(dp):: zp
!=============================================================================
zp=u1+xc(3); xs=xc(1:2)/zp
end subroutine xctoxs
