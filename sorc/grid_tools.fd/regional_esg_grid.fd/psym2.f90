!> @file
!! @brief Matrix routines.
!! @author R. J. Purser @date September 2018

!> A suite of routines to perform the eigen-decomposition of symmetric 2*2
!! matrices and to deliver basic analytic functions, and the
!! derivatives of these functions, of such matrices. In addition, we
!! include a simple cholesky routine.
!!
!! @author R. J. Purser
module psym2
use pkind, only: spi,dp
use pietc, only: u0,u1,o2
implicit none
private
public:: eigensym2,invsym2,sqrtsym2,expsym2,logsym2,id2222,chol2

real(dp),dimension(2,2,2,2):: id !< ID.
data id/u1,u0,u0,u0, u0,o2,o2,u0,  u0,o2,o2,u0, u0,u0,u0,u1/! Effective identity

interface eigensym2;   module procedure eigensym2,eigensym2d; end interface
interface invsym2;     module procedure invsym2,invsym2d;     end interface
interface sqrtsym2;    module procedure sqrtsym2,sqrtsym2d;   end interface
interface sqrtsym2d_e; module procedure sqrtsym2d_e;          end interface
interface sqrtsym2d_t; module procedure sqrtsym2d_t;          end interface
interface expsym2;     module procedure expsym2,expsym2d;     end interface
interface expsym2d_e;  module procedure expsym2d_e;           end interface
interface expsym2d_t;  module procedure expsym2d_t;           end interface
interface logsym2;     module procedure logsym2,logsym2d;     end interface
interface id2222;      module procedure id2222;               end interface
interface chol2;       module procedure chol2;                end interface

contains

!> Get the orthogonal eigenvectors, vv, and diagonal matrix of
!! eigenvalues, oo, of the symmetric 2*2 matrix, em.
!!
!! @param[in] em symmetric 2*2 matrix
!! @param[out] vv orthogonal eigenvectors
!! @param[out] oo diagonal matrix of eigenvalues
!! @author R. J. Purser
subroutine eigensym2(em,vv,oo)!                                    [eigensym2]
implicit none
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: vv,oo
real(dp):: a,b,c,d,e,f,g,h
a=em(1,1); b=em(1,2); c=em(2,2)
d=a*c-b*b! <- det(em)
e=(a+c)*o2; f=(a-c)*o2
h=sqrt(f**2+b**2)
g=sqrt(b**2+(h+abs(f))**2)
if    (g==u0)then; vv(:,1)=(/u1,u0/)
elseif(f> u0)then; vv(:,1)=(/h+f,b/)/g
else             ; vv(:,1)=(/b,h-f/)/g
endif      
vv(:,2)=(/-vv(2,1),vv(1,1)/)
oo=matmul(transpose(vv),matmul(em,vv))
oo(1,2)=u0; oo(2,1)=u0
end subroutine eigensym2

!> For a symmetric 2*2 matrix, em, return the normalized eigenvectors,
!! vv, and the diagonal matrix of eigenvalues, oo. If the two
!! eigenvalues are equal, proceed no further and raise the logical
!! failure flag, ff, to .true.; otherwise, return with vvd=d(vv)/d(em)
!! and ood=d(oo)/d(em) and ff=.false., and maintain the symmetries
!! between the last two of the indices of these derivatives.
!!
!! @param[in] em symmetric 2*2 matrix
!! @param[out] vv normalized eigenvectors
!! @param[out] oo diagonal matrix of eigenvalues
!! @param[out] vvd vvd=d(vv)/d(em)
!! @param[out] ood ood=d(oo)/d(em)
!! @param[out] ff logical failure flag
!! @author R. J. Purser
subroutine eigensym2d(em,vv,oo,vvd,ood,ff)!                        [eigensym2]
implicit none
real(dp),dimension(2,2),    intent(in ):: em
real(dp),dimension(2,2),    intent(out):: vv,oo
real(dp),dimension(2,2,2,2),intent(out):: vvd,ood
logical,                    intent(out):: ff
real(dp),dimension(2,2):: emd,tt,vvr
real(dp)               :: oodif,dtheta
integer(spi)           :: i,j
call eigensym2(em,vv,oo); vvr(1,:)=-vv(2,:); vvr(2,:)=vv(1,:)
oodif=oo(1,1)-oo(2,2); ff=oodif==u0; if(ff)return
ood=0
vvd=0
do j=1,2
   do i=1,2
      emd=0
      if(i==j)then
         emd(i,j)=u1
      else
         emd(i,j)=o2; emd(j,i)=o2
      endif
      tt=matmul(transpose(vv),matmul(emd,vv))
      ood(1,1,i,j)=tt(1,1)
      ood(2,2,i,j)=tt(2,2)
      dtheta=tt(1,2)/oodif
      vvd(:,:,i,j)=vvr*dtheta
   enddo
enddo
end subroutine eigensym2d

!> Get the inverse of a 2*2 matrix (need not be symmetric in this
!! case).
!!
!! @param[in] em 2*2 matrix
!! @param[out] z inverse of a 2*2 matrix
!! @author R. J. Purser
subroutine invsym2(em,z)!                                            [invsym2]
implicit none
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: z
real(dp):: detem
z(1,1)=em(2,2); z(2,1)=-em(2,1); z(1,2)=-em(1,2); z(2,2)=em(1,1)
detem=em(1,1)*em(2,2)-em(2,1)*em(1,2)
z=z/detem
end subroutine invsym2

!> Get the inverse, z,of a 2*2 symmetric matrix, em, and its
!! derivative, zd, with respect to symmetric variations of its
!! components. I.e., for a symmetric infinitesimal change, delta_em,
!! in em, the resulting infinitesimal change in z would be:
!! <pre>delta_z(i,j) = matmul(zd(i,j,:,:),delta_em)</pre>
!!
!! @param[in] em 2*2 symmetric matrix
!! @param[out] z inverse of a 2*2 symmetric matrix
!! @param[out] zd derivative of the 2*2 symmetric matrix
!! @author R. J. Purser
subroutine invsym2d(em,z,zd)!                                        [invsym2]
implicit none
real(dp),dimension(2,2)    ,intent(in ):: em
real(dp),dimension(2,2)    ,intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
integer(spi):: k,l
call invsym2(em,z)
call id2222(zd)
do l=1,2; do k=1,2
   zd(:,:,k,l)=-matmul(matmul(z,zd(:,:,k,l)),z)
enddo;    enddo
end subroutine invsym2d

!> Get the sqrt of a symmetric positive-definite 2*2 matrix.
!!
!! @param[in] em 2*2 symmetric matrix
!! @param[out] z sqrt of a symmetric positive-definite 2*2 matrix
!! @author R. J. Purser
subroutine sqrtsym2(em,z)!                                          [sqrtsym2]
implicit none
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: z
real(dp),dimension(2,2):: vv,oo
integer(spi)           :: i
call eigensym2(em,vv,oo)
do i=1,2
if(oo(i,i)<0)stop 'In sqrtsym2; matrix em is not non-negative'
oo(i,i)=sqrt(oo(i,i)); enddo
z=matmul(vv,matmul(oo,transpose(vv)))
end subroutine sqrtsym2

!> General routine to evaluate z=sqrt(x), and the symmetric
!! derivative, zd = dz/dx, where x is a symmetric 2*2
!! positive-definite matrix. If the eigenvalues are very close
!! together, extract their geometric mean for "preconditioning" a
!! scaled version, px, of x, whose sqrt, and hence its derivative, can
!! be easily obtained by the series expansion method. Otherwise, use
!! the eigen-method (which entails dividing by the difference in the
!! eignevalues to get zd, and which therefore fails when the
!! eigenvalues become too similar).
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z sqrt(x) result
!! @param[out] zd symmetric derivative 
!! @author R. J. Purser
subroutine sqrtsym2d(x,z,zd)!                                        [sqrtsym2]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
real(dp),dimension(2,2):: px
real(dp)               :: rdetx,lrdetx,htrpxs,q
rdetx=sqrt(x(1,1)*x(2,2)-x(1,2)*x(2,1)) ! <- sqrt(determinant of x)
lrdetx=sqrt(rdetx)
px=x/rdetx                  ! - preconditioned x (has unit determinant)
htrpxs= ((px(1,1)+px(2,2))/2)**2 ! - {half-trace-px}-squared
q=htrpxs-u1
if(q<.05_dp)then               ! - Taylor expansion method
   call sqrtsym2d_t(px,z,zd)
   z=z*lrdetx; zd=zd/lrdetx
else
   call sqrtsym2d_e(x,z,zd) ! - Eigen-method
endif
end subroutine sqrtsym2d

!> Eigen-method.
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z sqrt(x) result
!! @param[out] zd symmetric derivative
!! @author R. J. Purser
subroutine sqrtsym2d_e(x,z,zd)!                                  [sqrtsym2d_e]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
real(dp),dimension(2,2,2,2):: vvd,ood
real(dp),dimension(2,2)    :: vv,oo,oori,tt
integer(spi)               :: i,j
logical                    :: ff
call eigensym2(x,vv,oo,vvd,ood,ff)
z=u0; z(1,1)=sqrt(oo(1,1)); z(2,2)=sqrt(oo(2,2))
z=matmul(matmul(vv,z),transpose(vv))
oori=u0; oori(1,1)=u1/sqrt(oo(1,1)); oori(2,2)=u1/sqrt(oo(2,2))
do j=1,2
do i=1,2
   tt=matmul(vvd(:,:,i,j),transpose(vv))
   zd(:,:,i,j)=o2*matmul(matmul(matmul(vv,ood(:,:,i,j)),oori),transpose(vv))&
        +matmul(tt,z)-matmul(z,tt)
enddo
enddo
end subroutine sqrtsym2d_e

!> Use the Taylor-series method (eigenvalues both fairly close to unity).
!! For a 2*2 positive definite symmetric matrix x, try to get both the
!! z=sqrt(x) and dz/dx using the binomial-expansion method applied to
!! the intermediate matrix,
!! <pre>r = (x-1). ie z=sqrt(x) = (1+r)^{1/2} = I + (1/2)*r -(1/8)*r^2 ...
!!  + [(-)^n *(2n)!/{(n+1)! * n! *2^{2*n-1}} ]*r^{n+1}</pre>
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z sqrt(x) result
!! @param[out] zd symmetric derivative
!! @author R. J. Purser
subroutine sqrtsym2d_t(x,z,zd)!                                  [sqrtsym2d_t]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
integer(spi),parameter     :: nit=300 ! number of iterative increments allowed
real(dp),parameter         :: crit=1.e-17
real(dp),dimension(2,2)    :: r,rp,rd,rpd
real(dp)                   :: c
integer(spi)               :: i,j,n
r=x; r(1,1)=x(1,1)-1; r(2,2)=x(2,2)-1
z=u0; z(1,1)=u1; z(2,2)=u1
rp=r
c=o2
do n=0,nit
   z=z+c*rp
   rp=matmul(rp,r)
   if(sum(abs(rp))<crit)exit
   c=-c*(n*2+1)/(2*(n+2))
enddo
do j=1,2; do i=1,2
   rd=id(:,:,i,j)
   rpd=rd
   zd(:,:,i,j)=0
   rp=r
   c=o2
   do n=0,nit
      zd(:,:,i,j)=zd(:,:,i,j)+c*rpd
      rpd=matmul(rd,rp)+matmul(r,rpd)
      rp=matmul(r,rp)
      if(sum(abs(rp))<crit)exit
      c=-c*(n*2+1)/(2*(n+2))
   enddo
enddo; enddo
end subroutine sqrtsym2d_t

!> Get the exp of a symmetric 2*2 matrix.
!!
!! @param[in] em symmetric 2*2 matrix
!! @param[out] expem exp of a symmetric 2*2 matrix
!! @author R. J. Purser
subroutine expsym2(em,expem)!                                        [expsym2]
implicit none
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: expem
real(dp),dimension(2,2):: vv,oo
integer(spi)           :: i
call eigensym2(em,vv,oo)
do i=1,2; oo(i,i)=exp(oo(i,i)); enddo
expem=matmul(vv,matmul(oo,transpose(vv)))
end subroutine expsym2

!> Get the exp of a symmetric 2*2 matrix, and its symmetric derivative.
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z exp of symmetric 2*2 matrix x
!! @param[out] zd symmetric derivative wrt x of exp of x
!! @author R. J. Purser
subroutine expsym2d(x,z,zd)!                                         [expsym2]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
real(dp),dimension(2,2):: px
real(dp)               :: trxh,detpx
trxh=(x(1,1)+x(2,2))*o2
px=x;px(1,1)=x(1,1)-trxh;px(2,2)=x(2,2)-trxh
detpx=abs(px(1,1)*px(2,2)-px(1,2)*px(2,1))
if(detpx<.1_dp)then; call expsym2d_t(px,z,zd)
else               ; call expsym2d_e(px,z,zd)
endif
z=z*exp(trxh)
end subroutine expsym2d

!> Get the exponential and its symmetric derivative for a symmetric 2*2 matrix
!! using eigen-decomposition.
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z exp of symmetrix matrix x
!! @param[out] zd symmetric derivative of z wrt x
!! @author R. J. Purser
subroutine expsym2d_e(x,z,zd)!                                    [expsym2d_e]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
real(dp),dimension(2,2,2,2):: vvd,ood
real(dp),dimension(2,2)    :: vv,oo,ooe,tt
integer(spi)               :: i,j
logical                    :: ff
call eigensym2(x,vv,oo,vvd,ood,ff)
z=u0; z(1,1)=exp(oo(1,1)); z(2,2)=exp(oo(2,2))
z=matmul(matmul(vv,z),transpose(vv))
ooe=u0; ooe(1,1)=exp(oo(1,1)); ooe(2,2)=exp(oo(2,2))
do j=1,2
do i=1,2
   tt=matmul(vvd(:,:,i,j),transpose(vv))
   zd(:,:,i,j)=matmul(matmul(matmul(vv,ood(:,:,i,j)),ooe),transpose(vv))&
        +matmul(tt,z)-matmul(z,tt)
enddo
enddo
end subroutine expsym2d_e

!> Use the Taylor-series method (eigenvalues both fairly close to
!! zero). For a 2*2 symmetric matrix x, try to get both the z=exp(x)
!! and dz/dx using the Taylor series expansion method.
!!
!! @param[in] x symmetric 2*2 positive-definite matrix
!! @param[out] z Taylor series expansion method exp(x)
!! @param[out] zd symmetric derivative
!! @author R. J. Purser
subroutine expsym2d_t(x,z,zd)!                                    [expsym2d_t]
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
integer(spi),parameter     :: nit=100 ! number of iterative increments allowed
real(dp),parameter         :: crit=1.e-17_dp
real(dp),dimension(2,2)    :: xp,xd,xpd
real(dp)                   :: c
integer(spi)               :: i,j,n
z=0; z(1,1)=u1; z(2,2)=u1
xp=x
c=u1
do n=1,nit
   z=z+c*xp
   xp=matmul(xp,x)
   if(sum(abs(xp))<crit)exit
   c=c/(n+1)
enddo
do j=1,2; do i=1,2
   xd=id(:,:,i,j)
   xpd=xd
   zd(:,:,i,j)=0
   xp=x
   c=u1
   do n=1,nit
      zd(:,:,i,j)=zd(:,:,i,j)+c*xpd
      xpd=matmul(xd,xp)+matmul(x,xpd)
      xp=matmul(x,xp)
      if(sum(abs(xpd))<crit)exit
      c=c/(n+1)
   enddo
enddo; enddo
end subroutine expsym2d_t

!> Get the log of a symmetric positive-definite 2*2 matrix.
!!
!! @param[in] em symmetric 2*2 matrix
!! @param[out] logem log of a symmetric positive-definite 2*2 matrix
!! @author R. J. Purser
subroutine logsym2(em,logem)!                                        [logsym2]
implicit none
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: logem
real(dp),dimension(2,2):: vv,oo
integer(spi)           :: i
call eigensym2(em,vv,oo)
do i=1,2
   if(oo(i,i)<=u0)stop 'In logsym2; matrix em is not positive definite'
   oo(i,i)=log(oo(i,i))
enddo
logem=matmul(vv,matmul(oo,transpose(vv)))
end subroutine logsym2

!> General routine to evaluate the logarithm, z=log(x), and the
!! symmetric derivative, zd = dz/dx, where x is a symmetric 2*2
!! positive-definite matrix.
!!
!! @param[in] zd the symmetric derivative
!! @param[out] x a symmetric 2*2 positive-definite matrix
!! @param[out] z evaluate the logarithm log(x)
!! @author R. J. Purser
subroutine logsym2d(x,z,zd)!                                         [logsym2]
use pfun, only: sinhox
implicit none
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
real(dp),dimension(2,2):: vv,oo,d11,d12,d22,pqr
real(dp)               :: c,s,cc,cs,ss,c2h,p,q,r,lp,lq,L
integer(spi)           :: i
call eigensym2(x,vv,oo)
if(oo(1,1)<=u0 .or. oo(2,2)<=u0)stop 'In logsym2; matrix x is not positive definite'
c=vv(1,1); s=vv(1,2); cc=c*c; cs=c*s; ss=s*s; c2h=(cc-ss)*o2
p=u1/oo(1,1); q=u1/oo(2,2); lp=log(p); lq=log(q); oo(1,1)=-lp; oo(2,2)=-lq
z=matmul(vv,matmul(oo,transpose(vv)))
L=(lp-lq)*o2; r=sqrt(p*q)/sinhox(L)
d11(1,:)=(/ cc,cs /); d11(2,:)=(/ cs,ss/)
d12(1,:)=(/-cs,c2h/); d12(2,:)=(/c2h,cs/)
d22(1,:)=(/ ss,-cs/); d22(2,:)=(/-cs,cc/)
pqr(1,:)=(/p,r/)    ; pqr(2,:)=(/r,q/)
zd(:,:,1,1)=matmul(vv,matmul(d11*pqr,transpose(vv)))
zd(:,:,1,2)=matmul(vv,matmul(d12*pqr,transpose(vv)))
zd(:,:,2,2)=matmul(vv,matmul(d22*pqr,transpose(vv)))
zd(:,:,2,1)=zd(:,:,1,2)
end subroutine logsym2d

!> General routine for a symmetrized 4th-rank tensor that acts as
!! an effective identity for operations on symmetric matrices.
!!
!! @param[out] em symmetrized effective identity in space of symmetrix matrices.
!! @author R. J. Purser
subroutine id2222(em)!                                                [id2222]
implicit none
real(dp),dimension(2,2,2,2),intent(out):: em
real(dp),dimension(2,2,2,2)            :: id
!data id/u1,u0,u0,u0, u0,o2,o2,u0,  u0,o2,o2,u0, u0,u0,u0,u1/! Effective identity
em=id
end subroutine id2222

!> Return the cholesky lower triangular factor, C, of the 2X2
!! symmetric matrix, S, or raise the failure flag, FF, if S is not
!! positive-definite.
!!
!! @param[in] s 2X2 symmetric matrix
!! @param[out] c cholesky lower triangular factor
!! @param[out] ff raise the failure flag
!! @author R. J. Purser
subroutine chol2(s,c,ff)!                                            [chol2]
use pietc, only: u0
implicit none
real(dp),dimension(2,2),intent(in ):: s
real(dp),dimension(2,2),intent(out):: c
logical                ,intent(out):: ff
real(dp):: r
ff=s(1,1)<=u0; if(ff)return
c(1,1)=sqrt(s(1,1))
c(1,2)=u0
c(2,1)=s(2,1)/c(1,1)
r=s(2,2)-c(2,1)**2
ff=r<=u0;      if(ff)return
c(2,2)=sqrt(r)
end subroutine chol2

end module psym2

