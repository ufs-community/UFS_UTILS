!                                        ***********************************
!                                        *    module psym2                 *
!                                        *    R. J. Purser                 *
!                                        *    NOAA/NCEP/EMC September 2018 *
!                                        *    jim.purser@noaa.gov          *
!                                        ***********************************
!
! A suite of routines to perform the eigen-decomposition of symmetric 2*2
! matrices and to deliver basic analytic functions, and the derivatives
! of these functions, of such matrices.
!
! DIRECT DEPENDENCIES
! Module: pkind, pietc
!
!=============================================================================
module psym2
!=============================================================================
use pkind, only: dp
use pietc, only: u0,u1,o2
implicit none
private
public:: eigensym2,invsym2,sqrtsym2,expsym2,logsym2,id2222

real(dp),dimension(2,2,2,2):: id
data id/u1,u0,u0,u0, u0,o2,o2,u0,  u0,o2,o2,u0, u0,u0,u0,u1/! Effective identity

interface eigensym2;   module procedure eigensym2,eigensym2d; end interface
interface invsym2;     module procedure invsym2,invsym2d;     end interface
interface sqrtsym2;    module procedure sqrtsym2,sqrtsym2d;   end interface
interface sqrtsym2d_e; module procedure sqrtsym2d_e;          end interface
interface sqrtsym2d_t; module procedure sqrtsym2d_t;          end interface
interface expsym2;     module procedure expsym2,expsym2d;     end interface
interface expsym2d_e;  module procedure expsym2d_e;           end interface
interface expsym2d_t;  module procedure expsym2d_t;           end interface
interface logsym2;     module procedure logsym2,logsym2d ;    end interface
interface logsym2d_e;  module procedure logsym2d_e;           end interface
interface logsym2d_t;  module procedure logsym2d_t;           end interface
interface id2222;      module procedure id2222;               end interface

contains

!=============================================================================
subroutine eigensym2(em,vv,oo)!                                    [eigensym2]
!=============================================================================
! Get the orthogonal eigenvectors, vv, and diagonal matrix of eigenvalues, oo, 
! of the symmetric 2*2 matrix, em.
!=============================================================================
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: vv,oo
!-----------------------------------------------------------------------------
real(dp):: a,b,c,d,e,f,g,h
!=============================================================================
a=em(1,1); b=em(1,2); c=em(2,2)
d=a*c-b*b! <- det(em)
e=(a+c)/2; f=(a-c)/2
h=sqrt(f**2+b**2)
g=sqrt(b**2+(h+abs(f))**2)
if    (g==0)then; vv(:,1)=(/u1,u0/)
elseif(f> 0)then; vv(:,1)=(/h+f,b/)/g
else            ; vv(:,1)=(/b,h-f/)/g
endif      
vv(:,2)=(/-vv(2,1),vv(1,1)/)
oo=matmul(transpose(vv),matmul(em,vv))
oo(1,2)=u0; oo(2,1)=u0
end subroutine eigensym2
!=============================================================================
subroutine eigensym2d(em,vv,oo,vvd,ood,ff)!                        [eigensym2]
!=============================================================================
! For a symmetric 2*2 matrix, em, return the normalized eigenvectors, vv, and
! the diagonal matrix of eigenvalues, oo. If the two eigenvalues are equal,
! proceed no further and raise the logical failure flagg, ff, to .true.;
! otherwise, return with vvd=d(vv)/d(em) and ood=d(oo)/d(em) and ff=.false.,
! and maintain the symmetries between the last two of the indices of 
! these derivatives.
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: em
real(dp),dimension(2,2),    intent(out):: vv,oo
real(dp),dimension(2,2,2,2),intent(out):: vvd,ood
logical,                    intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: emd,tt,vvr
real(dp)               :: oodif,dtheta
integer                :: i,j
!=============================================================================
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

!=============================================================================
subroutine invsym2(em,z)!                                            [invsym2]
!=============================================================================
! Get the inverse of a 2*2 matrix (need not be symmetric in this case).
!=============================================================================
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: z
!-----------------------------------------------------------------------------
real(dp):: detem
!=============================================================================
z(1,1)=em(2,2); z(2,1)=-em(2,1); z(1,2)=-em(1,2); z(2,2)=em(1,1)
detem=em(1,1)*em(2,2)-em(2,1)*em(1,2)
z=z/detem
end subroutine invsym2
!=============================================================================
subroutine invsym2d(em,z,zd)!                                        [invsym2]
!=============================================================================
! Get the inverse, z,of a 2*2 symmetric matrix, em, and its derivative, zd, 
! with respect to symmetric variations of its components. I.e., for a 
! symmetric infinitesimal change, delta_em, in em, the resulting
! infinitesimal change in z would be:
! delta_z(i,j) = matmul(zd(i,j,:,:),delta_em)
!=============================================================================
real(dp),dimension(2,2)    ,intent(in ):: em
real(dp),dimension(2,2)    ,intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
integer:: k,l
!=============================================================================
call invsym2(em,z)
call id2222(zd)
do l=1,2; do k=1,2
   zd(:,:,k,l)=-matmul(matmul(z,zd(:,:,k,l)),z)
enddo;    enddo
end subroutine invsym2d

!=============================================================================
subroutine sqrtsym2(em,z)!                                          [sqrtsym2]
!=============================================================================
! Get the sqrt of a symmetric positive-definite 2*2 matrix
!=============================================================================
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: z
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: vv,oo
integer                :: i
!=============================================================================
call eigensym2(em,vv,oo)
do i=1,2
if(oo(i,i)<0)stop 'In sqrtsym2; matrix em is not non-negative'
oo(i,i)=sqrt(oo(i,i)); enddo
z=matmul(vv,matmul(oo,transpose(vv)))
end subroutine sqrtsym2

!=============================================================================
subroutine sqrtsym2d(x,z,zd)!                                        [sqrtsym2]
!=============================================================================
! General routine to evaluate z=sqrt(x),  and the symmetric
! derivative, zd = dz/dx, where x is a symmetric 2*2 positive-definite
! matrix. If the eigenvalues are very close together, extract their
! geometric mean for "preconditioning" a scaled version, px, of x, whose
! sqrt, and hence its derivative, can be easily obtained by the series
! expansion method. Otherwise, use the eigen-method (which entails dividing
! by the difference in the eignevalues to get zd, and which therefore
! fails when the eigenvalues become too similar).
!============================================================================= 
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: px
real(dp)               :: rdetx,lrdetx,htrpxs,q
!=============================================================================
rdetx=sqrt(x(1,1)*x(2,2)-x(1,2)*x(2,1)) ! <- sqrt(determinant of x)
lrdetx=sqrt(rdetx)
px=x/rdetx                  ! <- preconditioned x (has unit determinant)
htrpxs= ((px(1,1)+px(2,2))/2)**2 ! <- {half-trace-px}-squared
q=htrpxs-u1
if(q<.05)then               ! <- Taylor expansion method
   call sqrtsym2d_t(px,z,zd)
   z=z*lrdetx; zd=zd/lrdetx
else
   call sqrtsym2d_e(x,z,zd) ! <- Eigen-method
endif
end subroutine sqrtsym2d

!=============================================================================
subroutine sqrtsym2d_e(x,z,zd)!                                  [sqrtsym2d_e]
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),dimension(2,2,2,2):: vvd,ood
real(dp),dimension(2,2)    :: vv,oo,oori,tt
integer                    :: i,j
logical                    :: ff
!=============================================================================
call eigensym2(x,vv,oo,vvd,ood,ff)
z=u0; z(1,1)=sqrt(oo(1,1)); z(2,2)=sqrt(oo(2,2))
z=matmul(matmul(vv,z),transpose(vv))
oori=0; oori(1,1)=u1/sqrt(oo(1,1)); oori(2,2)=u1/sqrt(oo(2,2))
do j=1,2
do i=1,2
   tt=matmul(vvd(:,:,i,j),transpose(vv))
   zd(:,:,i,j)=o2*matmul(matmul(matmul(vv,ood(:,:,i,j)),oori),transpose(vv))&
        +matmul(tt,z)-matmul(z,tt)
enddo
enddo
end subroutine sqrtsym2d_e

!=============================================================================
subroutine sqrtsym2d_t(x,z,zd)!                                  [sqrtsym2d_t]
!=============================================================================
! Use the Taylor-series method (eigenvalues both fairly close to unity).
! For a 2*2 positive definite symmetric matrix x, try to get both the z=sqrt(x)
! and dz/dx using the binomial-expansion method applied to the intermediate
! matrix, r = (x-1). ie z=sqrt(x) = (1+r)^{1/2} = I + (1/2)*r -(1/8)*r^2 ...
!  + [(-)^n *(2n)!/{(n+1)! * n! *2^{2*n-1}} ]*r^{n+1}
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
integer,parameter          :: nit=300 ! number of iterative increments allowed
real(dp),parameter         :: crit=1.e-17
real(dp),dimension(2,2)    :: r,rp,rd,rpd
real(dp)                   :: c
integer                    :: i,j,n
!=============================================================================
r=x; r(1,1)=x(1,1)-1; r(2,2)=x(2,2)-1
z=0; z(1,1)=u1; z(2,2)=u1
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

!=============================================================================
subroutine expsym2(em,expem)!                                        [expsym2]
!=============================================================================
! Get the exp of a symmetric 2*2 matrix
!=============================================================================
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: expem
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: vv,oo
integer                :: i
!=============================================================================
call eigensym2(em,vv,oo)
do i=1,2; oo(i,i)=exp(oo(i,i)); enddo
expem=matmul(vv,matmul(oo,transpose(vv)))
end subroutine expsym2
!=============================================================================
subroutine expsym2d(x,z,zd)!                                         [expsym2]
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: px
real(dp)               :: trxh,detpx
!=============================================================================
trxh=(x(1,1)+x(2,2))/2
px=x;px(1,1)=x(1,1)-trxh;px(2,2)=x(2,2)-trxh
detpx=abs(px(1,1)*px(2,2)-px(1,2)*px(2,1))
if(detpx<.1)then; call expsym2d_t(px,z,zd)
else            ; call expsym2d_e(px,z,zd)
endif
z=z*exp(trxh)
end subroutine expsym2d

!=============================================================================
subroutine expsym2d_e(x,z,zd)!                                    [expsym2d_e]
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),dimension(2,2,2,2):: vvd,ood
real(dp),dimension(2,2)    :: vv,oo,ooe,tt
integer                    :: i,j
logical                    :: ff
!=============================================================================
call eigensym2(x,vv,oo,vvd,ood,ff)
z=u0; z(1,1)=exp(oo(1,1)); z(2,2)=exp(oo(2,2))
z=matmul(matmul(vv,z),transpose(vv))
ooe=0; ooe(1,1)=exp(oo(1,1)); ooe(2,2)=exp(oo(2,2))
do j=1,2
do i=1,2
   tt=matmul(vvd(:,:,i,j),transpose(vv))
   zd(:,:,i,j)=matmul(matmul(matmul(vv,ood(:,:,i,j)),ooe),transpose(vv))&
        +matmul(tt,z)-matmul(z,tt)
enddo
enddo
end subroutine expsym2d_e

!=============================================================================
subroutine expsym2d_t(x,z,zd)!                                    [expsym2d_t]
!=============================================================================
! Use the Taylor-series method (eigenvalues both fairly close to zero).
! For a 2*2 symmetric matrix x, try to get both the z=exp(x)
! and dz/dx using the Taylor series expansion method.
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
integer,parameter          :: nit=100 ! number of iterative increments allowed
real(dp),parameter         :: crit=1.e-17
real(dp),dimension(2,2)    :: xp,xd,xpd
real(dp)                   :: c
integer                    :: i,j,n
!=============================================================================
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

!=============================================================================
subroutine logsym2(em,logem)!                                        [logsym2]
!=============================================================================
! Get the log of a symmetric positive-definite 2*2 matrix
!=============================================================================
real(dp),dimension(2,2),intent(in ):: em
real(dp),dimension(2,2),intent(out):: logem
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: vv,oo
integer                :: i
!=============================================================================
call eigensym2(em,vv,oo)
do i=1,2
   if(oo(i,i)<=0)stop 'In logsym2; matrix em is not positive definite'
   oo(i,i)=log(oo(i,i))
enddo
logem=matmul(vv,matmul(oo,transpose(vv)))
end subroutine logsym2
!=============================================================================
subroutine logsym2d(x,z,zd)!                                         [logsym2]
!=============================================================================
! General routine to evaluate the logarithm, z=log(x),  and the symmetric
! derivative, zd = dz/dx, where x is a symmetric 2*2 positive-definite
! matrix. If the eigenvalues are very close together, extract their
! geometric mean for "preconditioning" a scaled version, px, of x, whose
! log, and hence its derivative, can be easily obtained by the series
! expansion method. Otherwise, use the eigen-method (which entails dividing
! by the difference in the eignevalues to get zd, and which therefore
! fails when the eigenvalues become too similar).
!=============================================================================  
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),parameter     :: o8=u1/8
real(dp),dimension(2,2):: px
real(dp)               :: rdetx,lrdetx,htrpxs,q
!=============================================================================
rdetx=sqrt(x(1,1)*x(2,2)-x(1,2)*x(2,1)) ! <- sqrt(determinant of x)
lrdetx=log(rdetx)
px=x/rdetx                  ! <- preconditioned x (has unit determinant)
htrpxs= ((px(1,1)+px(2,2))/2)**2 ! <- {half-trace-px}-squared
q=htrpxs-u1
if(q<o8)then               ! <- Taylor expansion method
   call logsym2d_t(px,z,zd)
   z(1,1)=z(1,1)+lrdetx; z(2,2)=z(2,2)+lrdetx; zd=zd/rdetx
else
   call logsym2d_e(x,z,zd) ! <- Eigen-method
endif
end subroutine logsym2d

!=============================================================================
subroutine logsym2d_e(x,z,zd)!                                    [logsym2d_e]
!=============================================================================
! Use the Eigen-decomposition method (eigenvalues not close together)
! For a 2*2 positive definite symmetric matrix x, try to get both the z=log(x)
! and dz/dx using the eigen-decomposition method.
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
real(dp),dimension(2,2,2,2):: vvd,ood
real(dp),dimension(2,2)    :: vv,oo,ooi,tt
integer                    :: i,j
logical                    :: ff
!==============================================================================
call eigensym2(x,vv,oo,vvd,ood,ff)
z=u0; z(1,1)=log(oo(1,1)); z(2,2)=log(oo(2,2))
z=matmul(matmul(vv,z),transpose(vv))
ooi=0; ooi(1,1)=u1/oo(1,1); ooi(2,2)=u1/oo(2,2)
do j=1,2
do i=1,2
   tt=matmul(vvd(:,:,i,j),transpose(vv))
   zd(:,:,i,j)=matmul(matmul(matmul(vv,ood(:,:,i,j)),ooi),transpose(vv))&
        +matmul(tt,z)-matmul(z,tt)
enddo
enddo
end subroutine logsym2d_e

!=============================================================================
subroutine logsym2d_t(x,z,zd)!                                    [logsym2d_t]
!=============================================================================
! Use the Taylor-series method (eigenvalues both fairly close to unity).
! For a 2*2 positive definite symmetric matrix x, try to get both the z=log(x)
! and dz/dx using the power-expansion method applied to the intermediate
! matrix, y = (x-1)*(x+1)^{-1}. ie z=log(x) = log(1+Y)-log(1-Y) = 
! 2*{Y + Y**3/3 + Y**5/5 + ....}, so, if Y' = dY/dX, then dZ/dX =
! 2*{Y'+(Y**2*Y' + Y*Y'*Y + Y'*Y**2)/3 + (Y**4*Y' + ...+Y'*Y**4)/5 + ... }.
! (Note, we are using the atanh function of Y.)
!=============================================================================
real(dp),dimension(2,2),    intent(in ):: x
real(dp),dimension(2,2),    intent(out):: z
real(dp),dimension(2,2,2,2),intent(out):: zd
!-----------------------------------------------------------------------------
integer,parameter          :: nit=30000 ! number of iterative increments allowed
real(dp),parameter         :: crit=1.e-17
real(dp),dimension(2,2)    :: xp,xm,xpi,y,yy,xd,yd,yyd,yp,ypd
real(dp)                   :: den
integer                    :: i,j,it
!=============================================================================
xp=x; xp(1,1)=x(1,1)+1; xp(2,2)=x(2,2)+1
xm=x; xm(1,1)=x(1,1)-1; xm(2,2)=x(2,2)-1
call invsym2(xp,xpi); y=matmul(xm,xpi); yy=matmul(y,y)
z=0
zd=0
do j=1,2; do i=1,2
   xd=id(:,:,i,j)
   yd=2*matmul(xpi,matmul(xd,xpi)); yp=y; ypd=yd 
   yyd=matmul(y,yd); yyd=yyd+transpose(yyd)
   if(j*i==1)z=z+y
   zd(:,:,i,j)=zd(:,:,i,j)+yd
   den=1
   do it=1,nit
      den=den+2
      ypd=matmul(yyd,yp)+matmul(yy,ypd)
      yp =matmul(yy,yp)
      if(i*j==1)z=z+yp/den
      zd(:,:,i,j)=zd(:,:,i,j)+ypd/den
      if(sum(abs(yp))<crit)exit
   enddo
enddo; enddo
z=2*z
zd=2*zd
end subroutine logsym2d_t

!=============================================================================
subroutine id2222(em)!                                                [id2222]
!=============================================================================
real(dp),dimension(2,2,2,2),intent(out):: em
real(dp),dimension(2,2,2,2)            :: id
data id/u1,u0,u0,u0, u0,o2,o2,u0,  u0,o2,o2,u0, u0,u0,u0,u1/! Effective identity
em=id
end subroutine id2222

end module psym2
