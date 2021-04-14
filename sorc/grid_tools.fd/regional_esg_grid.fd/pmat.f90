!> @file
!! @brief Utility routines for various linear inversions and Cholesky.
!!
!! @author R. J. Purser, NOAA/NCEP/EMC, Tsukasa Fujita, JMA.

!! Utility routines for various linear inversions and Cholesky.
!!
!! Originally, these routines were copies of the purely "inversion" members
!! of pmat1.f90 (a most extensive collection of matrix routines -- not just
!! inversions). As well as having both single and double precision versions
!! of each routine, these versions also make provision for a more graceful
!! termination in cases where the system matrix is detected to be 
!! essentially singular (and therefore noninvertible). This provision takes
!! the form of an optional "failure flag", FF, which is normally returned
!! as .FALSE., but is returned as .TRUE. when inversion fails.
!! In Sep 2012, these routines were collected together into pmat.f90 so
!! that all the main matrix routines could be in the same library, pmat.a.
!! 
!! @author R. J. Purser
module pmat
use pkind, only: spi,sp,dp,spc,dpc
use pietc, only: t,f
implicit none
private
public:: ldum,udlmm,inv,L1Lm,LdLm,invl,invu
interface swpvv;  module procedure sswpvv,dswpvv,cswpvv;        end interface
interface ldum
   module procedure sldum,dldum,cldum,sldumf,dldumf,cldumf;     end interface
interface udlmm
   module procedure sudlmm,dudlmm,cudlmm,sudlmv,dudlmv,cudlmv;  end interface
interface inv
   module procedure                                                           &
sinvmt, dinvmt, cinvmt, slinmmt, dlinmmt, clinmmt, slinmvt, dlinmvt, clinmvt, &
sinvmtf,dinvmtf,cinvmtf,slinmmtf,dlinmmtf,clinmmtf,slinmvtf,dlinmvtf,clinmvtf,&
iinvf
                                                               end interface
interface L1Lm;   module procedure sL1Lm,dL1Lm,sL1Lmf,dL1Lmf;  end interface
interface LdLm;   module procedure sLdLm,dLdLm,sLdLmf,dLdLmf;  end interface
interface invl;   module procedure sinvl,dinvl,slinlv,dlinlv;  end interface
interface invu;   module procedure sinvu,dinvu,slinuv,dlinuv;  end interface

contains

!> Swap a pair of single precision vectors
!!
!! @param[inout] d vector
!! @param[inout] e vector
!! @author R. J. Purser
subroutine sswpvv(d,e)!                                                [swpvv]
real(sp),    intent(inout) :: d(:), e(:)
real(sp)                   :: tv(size(d))
tv = d; d = e; e = tv
end subroutine sswpvv

!> Swap a pair of double precision vectors
!!
!! @param[inout] d vector
!! @param[inout] e vector
!! @author R. J. Purser
subroutine dswpvv(d,e)!                                                [swpvv]
real(dp), intent(inout) :: d(:), e(:)
real(dp)                :: tv(size(d))
tv = d; d = e; e = tv
end subroutine dswpvv

!> Swap a pair of complex vectors
!!
!! @param[inout] d vector
!! @param[inout] e vector
!! @author R. J. Purser
subroutine cswpvv(d,e)!                                                [swpvv]
complex(dpc),intent(inout) :: d(:), e(:)
complex(dpc)               :: tv(size(d))
tv = d; d = e; e = tv
end subroutine cswpvv

!> Invert single precision matrix in place
!!
!! @param[inout] a matrix
!! @author R. J. Purser
subroutine sinvmt(a)!                                                    [inv]
real(sp),dimension(:,:),intent(INOUT):: a
logical                              :: ff
call sinvmtf(a,ff)
if(ff)stop 'In sinvmt; Unable to invert matrix'
end subroutine sinvmt

!> Invert double precision matrix in place.
!!
!! @param[inout] a matrix
!! @author R. J. Purser
subroutine dinvmt(a)!                                                    [inv]
real(dp),dimension(:,:),intent(inout):: a
logical                              :: ff
call dinvmtf(a,ff)
if(ff)stop 'In dinvmt; Unable to invert matrix'
end subroutine dinvmt

!> Invert complex matrix in place.
!!
!! @param[inout] a matrix
!! @author R. J. Purser
subroutine cinvmt(a)!                                                    [inv]
complex(dpc),dimension(:,:),intent(inout):: a
logical                                  :: ff
call cinvmtf(a,ff)
if(ff)stop 'In cinvmt; Unable to invert matrix'
end subroutine cinvmt

!> Invert a single precision matrix in place, or flag if process fails.
!!
!! @param[inout] a matrix
!! @param[out] ff flag for error condition
!! @author R. J. Purser
subroutine sinvmtf(a,ff)!                                                [inv]
use pietc_s, only: u1
real(sp),dimension(:,:),intent(inout):: a
logical,                intent(  out):: ff 
integer(spi)                     :: m,i,j,jp,l
real(sp)                         :: d
integer(spi),dimension(size(a,1)):: ipiv
m=size(a,1)
if(m /= size(a,2))stop 'In sinvmtf; matrix passed to sinvmtf is not square'
! Perform a pivoted L-D-U decomposition on matrix a:
call sldumf(a,ipiv,d,ff)
if(ff)then
   print '(" In sinvmtf; failed call to sldumf")'
   return
endif
! Invert upper triangular portion U in place:
do i=1,m; a(i,i)=u1/a(i,i); enddo
do i=1,m-1
   do j=i+1,m; a(i,j)=-a(j,j)*dot_product(a(i:j-1,j),a(i,i:j-1)); enddo
enddo
! Invert lower triangular portion L in place:
do j=1,m-1; jp=j+1
   do i=jp,m; a(i,j)=-a(i,j)-dot_product(a(jp:i-1,j),a(i,jp:i-1)); enddo
enddo
!  Form the product of U**-1 and L**-1 in place
do j=1,m-1; jp=j+1
   do i=1,j; a(i,j)=a(i,j)+dot_product(a(jp:m,j),a(i,jp:m)); enddo
   do i=jp,m; a(i,j)=dot_product(a(i:m,j),a(i,i:m));         enddo
enddo
!  Permute columns according to ipiv
do j=m-1,1,-1; l=ipiv(j); call sswpvv(a(:,j),a(:,l)); enddo
end subroutine sinvmtf

!> Invert a double precision matrix in place, or flag if process fails
!!
!! @param[inout] a matrix
!! @param[out] ff flag for error condition
!! @author R. J. Purser
subroutine dinvmtf(a,ff)!                                                [inv]
real(dp),dimension(:,:),intent(inout):: a
logical,                intent(  out):: ff
integer(spi)                         :: m,i,j,jp,l
real(dp)                             :: d
integer(spi), dimension(size(a,1))   :: ipiv
m=size(a,1)
if(m /= size(a,2))stop 'In inv; matrix passed to dinvmtf is not square'
! Perform a pivoted L-D-U decomposition on matrix a:
call dldumf(a,ipiv,d,ff)
if(ff)then
   print '(" In dinvmtf; failed call to dldumf")'
   return
endif
! Invert upper triangular portion U in place:
do i=1,m; a(i,i)=1_dp/a(i,i); enddo
do i=1,m-1
   do j=i+1,m; a(i,j)=-a(j,j)*dot_product(a(i:j-1,j),a(i,i:j-1)); enddo
enddo
! Invert lower triangular portion L in place:
do j=1,m-1; jp=j+1
   do i=jp,m; a(i,j)=-a(i,j)-dot_product(a(jp:i-1,j),a(i,jp:i-1)); enddo
enddo
!  Form the product of U**-1 and L**-1 in place
do j=1,m-1; jp=j+1
   do i=1,j; a(i,j)=a(i,j)+dot_product(a(jp:m,j),a(i,jp:m)); enddo
   do i=jp,m; a(i,j)=dot_product(a(i:m,j),a(i,i:m));         enddo
enddo
!  Permute columns according to ipiv
do j=m-1,1,-1; l=ipiv(j); call dswpvv(a(:,j),a(:,l)); enddo
end subroutine dinvmtf

!> Invert a complex matrix in place, or flag if process fails.
!!
!! @param[inout] a matrix
!! @param[out] ff flag for error condition
!! @author R. J. Purser
subroutine cinvmtf(a,ff)!                                                [inv]
use pietc, only: c1
complex(dpc),dimension(:,:),intent(INOUT):: a
logical,                    intent(  OUT):: ff
integer(spi)                     :: m,i,j,jp,l
complex(dpc)                     :: d
integer(spi),dimension(size(a,1)):: ipiv
m=size(a,1)
if(m /= size(a,2))stop 'In inv; matrix passed to cinvmtf is not square'
! Perform a pivoted L-D-U decomposition on matrix a:
call cldumf(a,ipiv,d,ff)
if(ff)then
   print '(" In cinvmtf; failed call to cldumf")'
   return
endif
! Invert upper triangular portion U in place:
do i=1,m; a(i,i)=c1/a(i,i); enddo
do i=1,m-1
   do j=i+1,m; a(i,j)=-a(j,j)*sum(a(i:j-1,j)*a(i,i:j-1)); enddo
enddo
! Invert lower triangular portion L in place:
do j=1,m-1; jp=j+1
   do i=jp,m; a(i,j)=-a(i,j)-sum(a(jp:i-1,j)*a(i,jp:i-1)); enddo
enddo
!  Form the product of U**-1 and L**-1 in place
do j=1,m-1; jp=j+1
   do i=1,j; a(i,j)=a(i,j)+sum(a(jp:m,j)*a(i,jp:m)); enddo
   do i=jp,m; a(i,j)=sum(a(i:m,j)*a(i,i:m));         enddo
enddo
!  Permute columns according to ipiv
do j=m-1,1,-1; l=ipiv(j); call cswpvv(a(:,j),a(:,l)); enddo
end subroutine cinvmtf

!> Invert linear system with multiple right-hand side vectors.
!! Single precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors 
!! @author R. J. Purser
subroutine slinmmt(a,b)!                                                 [inv]
real(sp),dimension(:,:),intent(inout):: a,b
logical                              :: ff
call slinmmtf(a,b,ff)
if(ff)stop 'In slinmmt; unable to invert linear system'
end subroutine slinmmt

!> Invert linear system with multiple right-hand side vectors.
!! Double precision version
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors
!! @author R. J. Purser
subroutine dlinmmt(a,b)!                                                 [inv]
real(dp),dimension(:,:),intent(inout):: a,b
logical                              :: ff
call dlinmmtf(a,b,ff)
if(ff)stop 'In dlinmmt; unable to invert linear system'
end subroutine dlinmmt

!> Invert complex linear system with multiple right-hand side vectors.
!! Complex double precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors
!! @author R. J. Purser
subroutine clinmmt(a,b)!                                                 [inv]
complex(dpc),dimension(:,:),intent(inout):: a,b
logical                                  :: ff
call clinmmtf(a,b,ff)
if(ff)stop 'In clinmmt; unable to invert linear system'
end subroutine clinmmt

!> Invert linear system with multiple right-hand side vectors, or flag failure.
!! Single precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine slinmmtf(a,b,ff)!                                             [inv]
real(sp),   dimension(:,:),intent(inout):: a,b
logical,                   intent(  out):: ff
integer(spi),dimension(size(a,1))       :: ipiv
integer(spi)                            :: m
real(sp)                                :: d
m=size(a,1)
if(m /= size(a,2))stop 'In inv; matrix passed to slinmmtf is not square'
if(m /= size(b,1))&
     stop 'In inv; matrix and vectors in slinmmtf have unmatched sizes'
call sldumf(a,ipiv,d,ff)
if(ff)then
   print '("In slinmmtf; failed call to sldumf")'
   return
endif
call sudlmm(a,b,ipiv)
end subroutine slinmmtf

!> Invert linear system with multiple right-hand side vectors, or flag failure.
!! Double precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine dlinmmtf(a,b,ff)!                                             [inv]
real(dp),dimension(:,:),   intent(inout):: a,b
logical,                   intent(  out):: ff
integer(spi),dimension(size(a,1)):: ipiv
integer(spi):: m 
real(dp)    :: d
m=size(a,1)
if(m /= size(a,2))stop 'In inv; matrix passed to dlinmmtf is not square'
if(m /= size(b,1))&
     stop 'In inv; matrix and vectors in dlinmmtf have unmatched sizes'
call dldumf(a,ipiv,d,ff)
if(ff)then
   print '("In dlinmmtf; failed call to dldumf")'
   return
endif
call dudlmm(a,b,ipiv)
end subroutine dlinmmtf

!> Invert linear system with multiple right-hand side vectors, or flag failure.
!! Complex double precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vectors, output solution vectors
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine clinmmtf(a,b,ff)!                                             [inv]
complex(dpc),dimension(:,:),intent(INOUT):: a,b
logical,                    intent(  OUT):: ff
integer(spi),dimension(size(a,1)):: ipiv
integer(spi)                     :: m 
complex(dpc)                     :: d
m=size(a,1)
if(m /= size(a,2))stop 'In inv; matrix passed to dlinmmtf is not square'
if(m /= size(b,1))&
     stop 'In inv; matrix and vectors in dlinmmtf have unmatched sizes'
call cldumf(a,ipiv,d,ff)
if(ff)then
   print '("In clinmmtf; failed call to cldumf")'
   return
endif
call cudlmm(a,b,ipiv)
end subroutine clinmmtf

!> Invert linear system with single right-hand side vector.
!! Single precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @author R. J. Purser
subroutine slinmvt(a,b)!                                                 [inv]
real(sp),dimension(:,:),intent(inout):: a
real(sp),dimension(:),  intent(inout):: b
logical:: ff
call slinmvtf(a,b,ff)
if(ff)stop 'In slinmvt; matrix singular, unable to continue'
end subroutine slinmvt

!> Invert linear system with single right-hand side vector.
!! Double precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @author R. J. Purser
subroutine dlinmvt(a,b)!                                                 [inv]
real(dp),dimension(:,:),intent(inout):: a
real(dp),dimension(:),  intent(inout):: b
logical                              :: ff
call dlinmvtf(a,b,ff)
if(ff)stop 'In dlinmvt; matrix singular, unable to continue'
end subroutine dlinmvt

!> Invert linear system with single right-hand side vector.
!! Complex double precision version.
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @author R. J. Purser
subroutine clinmvt(a,b)!                                                 [inv]
complex(dpc),   dimension(:,:),intent(inout):: a
complex(dpc),   dimension(:),  intent(inout):: b
logical                                     :: ff
call clinmvtf(a,b,ff)
if(ff)stop 'In clinmvt; matrix singular, unable to continue'
end subroutine clinmvt

!> Invert linear system with single right-hand side vector
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine slinmvtf(a,b,ff)!                                             [inv]
real(sp),dimension(:,:),intent(inout):: a
real(sp),dimension(:),  intent(inout):: b
logical,                intent(  out):: ff
integer(spi),dimension(size(a,1))    :: ipiv
real(sp)                             :: d
if(size(a,1) /= size(a,2).or. size(a,1) /= size(b))&
     stop 'In inv; In slinmvtf; incompatible array dimensions'
call sldumf(a,ipiv,d,ff)
if(ff)then
   print '("In slinmvtf; failed call to sldumf")'
   return
endif
call sudlmv(a,b,ipiv) 
end subroutine slinmvtf

!> Invert linear system with single right-hand side vector
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine dlinmvtf(a,b,ff)!                                             [inv]
real(dp),dimension(:,:),intent(inout):: a
real(dp),dimension(:),  intent(inout):: b
logical,                intent(  out):: ff
integer(spi), dimension(size(a,1))   :: ipiv
real(dp)                             :: d
if(size(a,1) /= size(a,2).or. size(a,1) /= size(b))&
     stop 'In inv; incompatible array dimensions passed to dlinmvtf'
call dldumf(a,ipiv,d,ff)
if(ff)then
   print '("In dlinmvtf; failed call to dldumf")'
   return
endif
call dudlmv(a,b,ipiv)
end subroutine dlinmvtf

!> Invert complex linear system with single right-hand side vector
!!
!! @param[inout] a Invertible system matrix, destroyed on output
!! @param[inout] b input RHS vector, output solution vector
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine clinmvtf(a,b,ff)!                                             [inv]
complex(dpc),dimension(:,:),intent(inout):: a
complex(dpc),dimension(:),  intent(inout):: b
logical,                    intent(  out):: ff
integer, dimension(size(a,1))            :: ipiv
complex(dpc)                             :: d
if(size(a,1) /= size(a,2).or. size(a,1) /= size(b))&
     stop 'In inv; incompatible array dimensions passed to clinmvtf'
call cldumf(a,ipiv,d,ff)
if(ff)then
   print '("In clinmvtf; failed call to cldumf")'
   return
endif
call cudlmv(a,b,ipiv)
end subroutine clinmvtf

!> Invert integer square matrix, imat, if possible, but flag ff=.true.
!! if not possible. (Determinant of imat must be +1 or -1)
!!
!! @param[inout] imat integer square matrix
!! @param[out] ff error flag
!! @author R. J. Purser
subroutine iinvf(imat,ff)!                                               [inv]
integer(spi),dimension(:,:),intent(INOUT):: imat
logical,                    intent(  OUT):: ff
real(dp),parameter                           :: eps=1.e-6_dp
real(dp),dimension(size(imat,1),size(imat,1)):: dmat
integer(spi)                                 :: m,i,j
m=size(imat,1)
if(m /= size(imat,2))stop 'In inv; matrix passed to iinvf is not square'
dmat=imat; call inv(dmat,ff)
if(.not.ff)then
   do j=1,m
      do i=1,m
         imat(i,j)=nint(dmat(i,j)); if(abs(dmat(i,j)-imat(i,j))>eps)ff=t
      enddo
   enddo
endif
end subroutine iinvf

!> Perform L*D*U decomposition, with pivoting, of square matrix.
!! Single precision version.
!!
!! @param[inout] a input square matrix, output L,D,U factors
!! @param[out] d determinant sign change indicator (+1 or -1)
!! @param[out] ipiv vector of pivots
!! @author R. J. Purser
subroutine sldum(a,ipiv,d)!                                             [ldum]
real(sp),    intent(inout) :: a(:,:) 
real(sp),    intent(  out) :: d
integer(spi),intent(  out) :: ipiv(:)
logical:: ff
call sldumf(a,ipiv,d,ff)
if(ff)stop 'In sldum; matrix singular, unable to continue'
end subroutine sldum

!> Perform L*D*U decomposition, with pivoting, of square matrix.
!! Double precision version.
!!
!! @param[inout] a input square matrix, output L,D,U factors
!! @param[out] d determinant sign change indicator (+1 or -1)
!! @param[out] ipiv vector of pivots
!! @author R. J. Purser
subroutine dldum(a,ipiv,d)!                                             [ldum]
real(dp),    intent(inout) :: a(:,:) 
real(dp),    intent(  out) :: d
integer(spi),intent(  out) :: ipiv(:)
logical:: ff
call dldumf(a,ipiv,d,ff)
if(ff)stop 'In dldum; matrix singular, unable to continue'
end subroutine dldum

!> Perform L*D*U decomposition, with pivoting, of square matrix.
!! Complex double precision version.
!!
!! @param[inout] a input square matrix, output L,D,U factors
!! @param[out] d determinant sign change indicator (+1 or -1)
!! @param[out] ipiv vector of pivots
!! @author R. J. Purser
subroutine cldum(a,ipiv,d)!                                             [ldum]
complex(dpc),intent(inout) :: a(:,:) 
complex(dpc),intent(out  ) :: d
integer(spi),intent(out  ) :: ipiv(:)
logical:: ff
call cldumf(a,ipiv,d,ff)
if(ff)stop 'In cldum; matrix singular, unable to continue'
end subroutine cldum

!> Perform l-d-u decomposition of square matrix a in place with pivoting.
!! Single precision version.
!!
!! @param[inout] a  square matrix to be factorized
!! @param[out] ipiv vector encoding the pivoting sequence
!! @param[out] d    indicator for possible sign change of determinant
!! @param[out] ff:  failure flag, set to .true. when determinant of a vanishes.
!! @author R. J. Purser
subroutine sldumf(a,ipiv,d,ff)!                                         [ldum]
use pietc_s,only: u0,u1
real(sp),    intent(inout) :: a(:,:) 
real(sp),    intent(  out) :: d
integer(spi),intent(  out) :: ipiv(:)
logical,     intent(  out) :: ff
integer(spi):: m,i, j, jp, ibig, jm
real(sp)    :: s(size(a,1)),  aam, aa, abig,  ajj, ajji, aij
ff=f
m=size(a,1)
do i=1,m
  aam=u0
  do j=1,m
    aa=abs(a(i,j))
    if(aa > aam)aam=aa
  enddo
  if(aam == u0)then
    print '("In sldumf; row ",i6," of matrix vanishes")',i
    ff=t
    return
  endif
  s(i)=u1/aam
enddo
d=1_sp
ipiv(m)=m
do j=1,m-1
   jp=j+1
   abig=s(j)*abs(a(j,j))
   ibig=j
   do i=jp,m
      aa=s(i)*abs(a(i,j))
      if(aa > abig)then
         ibig=i
         abig=aa
      endif
   enddo
!  swap rows, recording changed sign of determinant
   ipiv(j)=ibig
   if(ibig /= j)then
      d=-d
      call sswpvv(a(j,:),a(ibig,:))
      s(ibig)=s(j)
   endif
   ajj=a(j,j)
   if(ajj == u0)then
      jm=j-1
      print '(" failure in sldumf:"/" matrix singular, rank=",i3)',jm
      ff=t
      return
   endif
   ajji=u1/ajj
   do i=jp,m
      aij=ajji*a(i,j)
      a(i,j)=aij
      a(i,jp:m) = a(i,jp:m) - aij*a(j,jp:m)
   enddo
enddo
end subroutine sldumf

!> Perform l-d-u decomposition of square matrix a in place with pivoting.
!! Double precision version.
!!
!! @param[inout] a  square matrix to be factorized
!! @param[out] ipiv vector encoding the pivoting sequence
!! @param[out] d    indicator for possible sign change of determinant
!! @param[out] ff:  failure flag, set to .true. when determinant of a vanishes.
!! @author R. J. Purser
subroutine dldumf(a,ipiv,d,ff)!                                         [ldum]
use pietc, only: u0,u1
real(dp),    intent(inout) :: a(:,:) 
real(dp),    intent(  out) :: d
integer,     intent(  out) :: ipiv(:)
logical(spi),intent(  out) :: ff
integer(spi)               :: m,i, j, jp, ibig, jm
real(dp)                   :: s(size(a,1)),  aam, aa, abig,  ajj, ajji, aij
ff=f
m=size(a,1)
do i=1,m
   aam=u0
   do j=1,m
      aa=abs(a(i,j))
      if(aa > aam)aam=aa
   enddo
   if(aam == u0)then
      print '("In dldumf;  row ",i6," of matrix vanishes")',i
      ff=t
      return
   endif
   s(i)=u1/aam
enddo
d=u1
ipiv(m)=m
do j=1,m-1
   jp=j+1
   abig=s(j)*abs(a(j,j))
   ibig=j
   do i=jp,m
      aa=s(i)*abs(a(i,j))
      if(aa > abig)then
         ibig=i
         abig=aa
      endif
   enddo
   !  swap rows, recording changed sign of determinant
   ipiv(j)=ibig
   if(ibig /= j)then
      d=-d
      call dswpvv(a(j,:),a(ibig,:))
      s(ibig)=s(j)
   endif
   ajj=a(j,j)
   if(ajj == u0)then
      jm=j-1
      print '(" Failure in dldumf:"/" matrix singular, rank=",i3)',jm
      ff=t
      return
   endif
   ajji=u1/ajj
   do i=jp,m
      aij=ajji*a(i,j)
      a(i,j)=aij
      a(i,jp:m) = a(i,jp:m) - aij*a(j,jp:m)
   enddo
enddo
end subroutine dldumf

!> Perform l-d-u decomposition of square matrix a in place with pivoting.
!! Complex double precision version.
!!
!! @param[inout] a  square matrix to be factorized
!! @param[out] ipiv vector encoding the pivoting sequence
!! @param[out] d    indicator for possible sign change of determinant
!! @param[out] ff:  failure flag, set to .true. when determinant of a vanishes.
!! @author R. J. Purser
subroutine cldumf(a,ipiv,d,ff)!                                         [ldum]
use pietc, only: u0,u1,c0,c1
complex(dpc), intent(inout)  :: a(:,:) 
complex(dpc), intent(  out)  :: d
integer(spi), intent(  out)  :: ipiv(:)
logical,      intent(  out)  :: ff
integer(spi)                 :: m,i, j, jp, ibig, jm
complex(dpc)                 :: ajj, ajji, aij
real(dp)                     :: aam,aa,abig
real(dp),dimension(size(a,1)):: s
ff=f
m=size(a,1)
do i=1,m
   aam=u0
   do j=1,m
      aa=abs(a(i,j))
      if(aa > aam)aam=aa
   enddo
   if(aam == u0)then
      print '("In cldumf;  row ",i6," of matrix vanishes")',i
      ff=t
      return
   endif
   s(i)=u1/aam
enddo
d=c1
ipiv(m)=m
do j=1,m-1
   jp=j+1
   abig=s(j)*abs(a(j,j))
   ibig=j
   do i=jp,m
      aa=s(i)*abs(a(i,j))
      if(aa > abig)then
         ibig=i
         abig=aa
      endif
   enddo
   !  swap rows, recording changed sign of determinant
   ipiv(j)=ibig
   if(ibig /= j)then
      d=-d
      call cswpvv(a(j,:),a(ibig,:))
      s(ibig)=s(j)
   endif
   ajj=a(j,j)
   if(ajj == c0)then
      jm=j-1
      print '(" Failure in cldumf:"/" matrix singular, rank=",i3)',jm
      ff=t
      return
   endif
   ajji=c1/ajj
   do i=jp,m
      aij=ajji*a(i,j)
      a(i,j)=aij
      a(i,jp:m) = a(i,jp:m) - aij*a(j,jp:m)
   enddo
enddo
end subroutine cldumf

!> Use l-u factors in A to back-substitute for several rhs in B, using
!! ipiv to define the pivoting permutation used in the l-u
!! decomposition.
!!
!! @param[in] a L-D-U factorization of linear system matrux
!! @param[inout] b rt-hand-sides vectors on input, corresponding
!! solutions on return
!! @param[in] ipiv vector encoding the pivoting sequence
!! @author R. J. Purser
subroutine sudlmm(a,b,ipiv)!                                           [udlmm]
use pietc_s, only: u1
integer(spi),dimension(:),  intent(in)    :: ipiv 
real(sp),    dimension(:,:),intent(in)    :: a 
real(sp),    dimension(:,:),intent(inout) :: b 
integer(spi):: m,i, k, l
real(sp)    :: s,aiii
m=size(a,1)
do k=1,size(b,2) !loop over columns of b
  do i=1,m
    l=ipiv(i)
    s=b(l,k)
    b(l,k)=b(i,k)
    s = s - sum(b(1:i-1,k)*a(i,1:i-1))
    b(i,k)=s
  enddo
  b(m,k)=b(m,k)/a(m,m)
  do i=m-1,1,-1
    aiii=u1/a(i,i)
    b(i,k) = b(i,k) - sum(b(i+1:m,k)*a(i,i+1:m))
    b(i,k)=b(i,k)*aiii
  enddo
enddo
end subroutine sudlmm

!> Use l-u factors in A to back-substitute for several rhs in B, using
!! ipiv to define the pivoting permutation used in the l-u
!! decomposition.
!!
!! @param[in] a square matrix to be factorized
!! @param[inout] b rt-hand-sides vectors on input, corresponding
!! solutions on return
!! @param[in] ipiv vector encoding the pivoting sequence
!! @author R. J. Purser
subroutine dudlmm(a,b,ipiv)!                                           [udlmm]
use pietc, only: u1
integer(spi),dimension(:),  intent(in   ) :: ipiv 
real(dp),    dimension(:,:),intent(in   ) :: a 
real(dp),    dimension(:,:),intent(inout) :: b 
integer(spi):: m,i, k, l
real(dp)    :: s,aiii
m=size(a,1)
do k=1, size(b,2)!loop over columns of b
   do i=1,m
      l=ipiv(i)
      s=b(l,k)
      b(l,k)=b(i,k)
      s = s - sum(b(1:i-1,k)*a(i,1:i-1))
      b(i,k)=s
   enddo
   b(m,k)=b(m,k)/a(m,m)
   do i=m-1,1,-1
      aiii=u1/a(i,i)
      b(i,k) = b(i,k) - sum(b(i+1:m,k)*a(i,i+1:m))
      b(i,k)=b(i,k)*aiii
   enddo
enddo
end subroutine dudlmm

!> Use l-u factors in A to back-substitute for several rhs in B, using
!! ipiv to define the pivoting permutation used in the l-u
!! decomposition.
!!
!! @param[in] a square matrix to be factorized
!! @param[inout] b rt-hand-sides vectors on input, corresponding
!! solutions on return
!! @param[in] ipiv vector encoding the pivoting sequence
!! @author R. J. Purser
subroutine cudlmm(a,b,ipiv)!                                           [udlmm]
use pietc, only: c1
integer(spi),dimension(:),  intent(in   ) :: ipiv 
complex(dpc),dimension(:,:),intent(in   ) :: a 
complex(dpc),dimension(:,:),intent(inout) :: b 
integer(spi):: m,i, k, l
complex(dpc):: s,aiii
m=size(a,1)
do k=1, size(b,2)!loop over columns of b
   do i=1,m
      l=ipiv(i)
      s=b(l,k)
      b(l,k)=b(i,k)
      s = s - sum(b(1:i-1,k)*a(i,1:i-1))
      b(i,k)=s
   enddo
   b(m,k)=b(m,k)/a(m,m)
   do i=m-1,1,-1
      aiii=c1/a(i,i)
      b(i,k) = b(i,k) - sum(b(i+1:m,k)*a(i,i+1:m))
      b(i,k)=b(i,k)*aiii
   enddo
enddo
end subroutine cudlmm

!> Use l-u factors in A to back-substitute for 1 rhs in B, using ipiv to
!! define the pivoting permutation used in the l-u decomposition.
!!
!! @param[in] a L-D-U factorization of linear system matrix
!! @param[inout] b right-hand-side vector on input, corresponding
!! solution on return
!! @param[in] ipiv vector encoding the pivoting sequence
!! @author R. J. Purser
subroutine sudlmv(a,b,ipiv)!                                           [udlmv]
use pietc_s, only: u1
integer(spi),dimension(:),  intent(in   ):: ipiv 
real(sp),    dimension(:,:),intent(in   ):: a 
real(sp),    dimension(:),  intent(inout):: b 
integer(spi):: m,i, l
real(sp)    :: s,aiii
m=size(a,1)
do i=1,m
   l=ipiv(i)
   s=b(l)
   b(l)=b(i)
   s = s - sum(b(1:i-1)*a(i,1:i-1))
   b(i)=s
enddo
b(m)=b(m)/a(m,m)
do i=m-1,1,-1
   aiii=u1/a(i,i)
   b(i) = b(i) - sum(b(i+1:m)*a(i,i+1:m))
   b(i)=b(i)*aiii
enddo
end subroutine sudlmv

!> Use l-u factors in A to back-substitute for 1 rhs in B, using ipiv to
!! define the pivoting permutation used in the l-u decomposition.
!!
!! @param[in] a square matrix to be factorized
!! @param[inout] b right-hand side vector on input, corresponding
!! solution on return
!! @param[in] ipiv array encoding the pivoting sequence
!! @author R. J. Purser
subroutine dudlmv(a,b,ipiv)!                                           [udlmv]
use pietc, only: u1
integer(spi),dimension(:),  intent(in   ) :: ipiv(:) 
real(dp),    dimension(:,:),intent(in   ) :: a(:,:) 
real(dp),    dimension(:),  intent(inout) :: b(:) 
integer(spi):: m,i, l
real(dp)    :: s,aiii
m=size(a,1)
do i=1,m
   l=ipiv(i)
   s=b(l)
   b(l)=b(i)
   s = s - sum(b(1:i-1)*a(i,1:i-1))
   b(i)=s
enddo
b(m)=b(m)/a(m,m)
do i=m-1,1,-1
   aiii=u1/a(i,i)
   b(i) = b(i) - sum(b(i+1:m)*a(i,i+1:m))
   b(i)=b(i)*aiii
enddo
end subroutine dudlmv

!> Use l-u factors in A to back-substitute for 1 rhs in B, using ipiv to
!! define the pivoting permutation used in the l-u decomposition.
!!
!! @param[in] a square matrix to be factorized
!! @param[inout] b right-hand side vector on input, corresponding 
!! solution on return
!! @param[in] ipiv array encoding the pivoting sequence
!! @author R. J. Purser
subroutine cudlmv(a,b,ipiv)!                                           [udlmv]
use pietc, only: c1
integer(spi),dimension(:),  intent(in   ) :: ipiv(:) 
complex(dpc),dimension(:,:),intent(in   ) :: a(:,:) 
complex(dpc),dimension(:),  intent(inout) :: b(:) 
integer(spi):: m,i, l
complex(dpc):: s,aiii
m=size(a,1)
do i=1,m
   l=ipiv(i)
   s=b(l)
   b(l)=b(i)
   s = s - sum(b(1:i-1)*a(i,1:i-1))
   b(i)=s
enddo
b(m)=b(m)/a(m,m)
do i=m-1,1,-1
   aiii=c1/a(i,i)
   b(i)= b(i) - sum(b(i+1:m)*a(i,i+1:m))
   b(i)=b(i)*aiii
enddo
end subroutine cudlmv

!> Cholesky, M -> L*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b Cholesky factor matrix.
!! @author R. J. Purser
subroutine sl1lm(a,b) !                                                 [l1lm]
real(sp),intent(in   ):: a(:,:)
real(sp),intent(inout):: b(:,:)
logical:: ff
call sl1lmf(a,b,ff)
if(ff)stop 'In sl1lm; matrix singular, unable to continue'
end subroutine sl1lm

!> Cholesky, M -> L*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b Cholesky factor matrix.
!! @author R. J. Purser
subroutine dl1lm(a,b) !                                                 [l1lm]
real(dp),intent(in   ):: a(:,:)
real(dp),intent(inout):: b(:,:)
logical:: ff
call dl1lmf(a,b,ff)
if(ff)stop 'In dl1lm; matrix singular, unable to continue'
end subroutine dl1lm

!> Cholesky, M -> L*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b Cholesky factor matrix.
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine sl1lmf(a,b,ff)!                                              [L1Lm] 
use pietc_s, only: u0
real(sp),intent(in   ):: a(:,:)
real(sp),intent(inout):: b(:,:)
logical, intent(  out):: ff
integer(spi):: m,j, jm, jp, i
real(sp)    :: s, bjji
m=size(a,1)
ff=f
do j=1,m
   jm=j-1
   jp=j+1
   s = a(j,j) - sum(b(j,1:jm)*b(j,1:jm))
   ff=(S <= u0)
   if(ff)then
      print '("sL1Lmf detects nonpositive a, rank=",i6)',jm
      return
   endif
   b(j,j)=sqrt(s)
   bjji=1_sp/b(j,j)
   do i=jp,m
      s = a(i,j) - sum(b(i,1:jm)*b(j,1:jm))
      b(i,j)=s*bjji
   enddo
   b(1:jm,j) = u0
enddo
end subroutine sl1lmf

!> Cholesky, M -> L*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b Cholesky factor matrix.
!! @param[out] ff failure flag
!! @author R. J. Purser
subroutine dl1lmf(a,b,ff) !                                             [L1Lm]
use pietc, only: u0,u1
real(dp),intent(in   ) :: a(:,:) 
real(dp),intent(inout) :: b(:,:) 
logical, intent(  out) :: ff
integer(spi):: m,j, jm, jp, i
real(dp)    :: s, bjji
m=size(a,1)
ff=f
do j=1,m
  jm=j-1
  jp=j+1
  s = a(j,j) - sum(b(j,1:jm)*b(j,1:jm))
  ff=(s <= u0)
  if(ff)then
     print '("dL1LMF detects nonpositive A, rank=",i6)',jm
     return
  endif
  b(j,j)=sqrt(s)
  bjji=u1/b(j,j)
  do i=jp,m
    s = a(i,j) - sum(b(i,1:jm)*b(j,1:jm))
    b(i,j)=s*bjji
  enddo
  b(1:jm,j) = u0
enddo
end subroutine dl1lmf

!> Modified Cholesky decompose Q --> L*D*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b output modified cholesky factor, L.
!! @param[out] d diagonal matrix, D.
!! @author R. J. Purser
subroutine sldlm(a,b,d)!                                                [LdLm]
real(sp),intent(in   ):: a(:,:)
real(sp),intent(inout):: b(:,:)
real(sp),intent(  out):: d(:)
logical:: ff
call sldlmf(a,b,d,ff)
if(ff)stop 'In sldlm; matrix singular, unable to continue'
end subroutine sldlm

!> Modified Cholesky decompose Q --> L*D*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b output modified cholesky factor, L.
!! @param[out] d diagonal matrix, D.
!! @author R. J. Purser
subroutine dldlm(a,b,d)!                                                [LdLm]
real(dp),intent(in   ):: a(:,:)
real(dp),intent(inout):: b(:,:)
real(dp),intent(  out):: d(:)
logical:: ff
call dldlmf(a,b,d,ff)
if(ff)stop 'In dldlm; matrix singular, unable to continue'
end subroutine dldlm

!> Modified Cholesky decompose Q --> L*D*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix
!! @param[inout] b modified cholesky factor, L.
!! @param[out] d diagonal matrix, D.
!! @param[out] ff error flag
!! @author R. J. Purser
subroutine sldlmf(a,b,d,ff) !                                           [LDLM]
use pietc_s, only: u0,u1
real(sp), intent(in   ):: a(:,:)
real(sp), intent(inout):: b(:,:)
real(sp), intent(  out):: d(:)
logical,  intent(  out):: ff
integer(spi):: m,j, jm, jp, i
real(sp)    :: bjji
m=size(a,1)
ff=f
do j=1,m
  jm=j-1
  jp=j+1
  d(j)=a(j,j) - sum(b(1:jm,j)*b(j,1:jm))
  b(j,j) = u1
  ff=(d(j) == u0)
  if(ff)then
     print '("In sldlmf; singularity of matrix detected")'
     print '("Rank of matrix: ",i6)',jm
     return
  endif
  bjji=u1/d(j)
  do i=jp,m
     b(j,i)=a(i,j) - dot_product(b(1:jm,j),b(i,1:jm))
     b(i,j)=b(j,i)*bjji
  enddo
  b(1:jm,j)=u0
enddo
end subroutine sldlmf

!> Modified Cholesky  Q --> L*D*U, U(i,j)=L(j,i)
!!
!! @param[in] a symmetric matrix.
!! @param[inout] b modified Cholesky factor, L.
!! @param[out] d diagonal matrix, D.
!! @param[out] ff error flag
!! @author R. J. Purser 
subroutine dldlmf(a,b,d,ff) !                                           [LDLM]
use pietc, only: u0,u1
real(dp), intent(IN   ) :: a(:,:)
real(dp), intent(INOUT) :: b(:,:)
real(dp), intent(  OUT) :: d(:)
logical,  intent(  OUT) :: ff
integer(spi):: m,j, jm, jp, i
real(dp)    :: bjji
m=size(a,1)
ff=f
do j=1,m; jm=j-1; jp=j+1
  d(j)=a(j,j) - sum(b(1:jm,j)*b(j,1:jm))
  b(j,j) = 1
  ff=(d(j) == u0)
  if(ff)then
     print '("In dldlmf; singularity of matrix detected")'
     print '("Rank of matrix: ",i6)',jm
     return
  endif
  bjji=u1/d(j)
  do i=jp,m
     b(j,i)=a(i,j) - dot_product(b(1:jm,j),b(i,1:jm))
     b(i,j)=b(j,i)*bjji
  enddo
  b(1:jm,j)=u0
enddo
end subroutine dldlmf

!> Invert the upper triangular matrix in place by transposing, calling
!! invl, and transposing again. Single precision version.
!!
!! @param[inout] a  upper triangular matrix.
!! @author R. J. Purser
subroutine sinvu(a)!                                                     [invu]
real(sp),dimension(:,:),intent(inout):: a
a=transpose(a); call sinvl(a); a=transpose(a)
end subroutine sinvu

!> Invert the upper triangular matrix in place by transposing, calling
!! invl, and transposing again. Double precision version.
!!
!! @param[inout] a  upper triangular matrix.
!! @author R. J. Purser
subroutine dinvu(a)!                                                     [invu]
real(dp),dimension(:,:),intent(inout):: a
a=transpose(a); call dinvl(a); a=transpose(a)
end subroutine dinvu

!> Invert lower triangular matrix in place. Single precision.
!!
!! @param[inout] a  lower triangular matrix.
!! @author R. J. Purser
subroutine sinvl(a)!                                                     [invl]
use pietc_s, only: u0,u1
real(sp), intent(inout) :: a(:,:) 
integer(spi):: m,j, i
m=size(a,1)
do j=m,1,-1
   a(1:j-1,j) = u0
   a(j,j)=u1/a(j,j)
   do i=j+1,m
      a(i,j)=-a(i,i)*sum(a(j:i-1,j)*a(i,j:i-1))
   enddo
enddo
end subroutine sinvl

!> Invert lower triangular matrix in place. Double precision.
!!
!! @param[inout] a  lower triangular matrix.
!! @author R. J. Purser
subroutine dinvl(a)!                                                     [invl]
use pietc, only: u0,u1
real(dp), intent(inout) :: a(:,:) 
integer(spi):: m,j, i
m=size(a,1)
do j=m,1,-1
   a(1:j-1,j) = u0
   a(j,j)=u1/a(j,j)
   do i=j+1,m
      a(i,j)=-a(i,i)*sum(a(j:i-1,j)*a(i,j:i-1))
   enddo
enddo
end subroutine dinvl

!> Solve linear system involving lower triangular system matrix.
!! Single precision version.
!!
!! @param[in] a lower triangular matrix.
!! @param[inout] u input RHS vector, output solution vector.
!! @author R. J. Purser
subroutine slinlv(a,u)!                                                  [invl]
real(sp),intent(in   ) :: a(:,:)
real(sp),intent(inout) :: u(:)
integer(spi):: i
if(size(a,1) /= size(a,2) .or. size(a,1) /= size(u))&
     stop 'In slinlv; incompatible array dimensions'
do i=1,size(u); u(i)=(u(i) - sum(u(:i-1)*a(i,:i-1)))/a(i,i); enddo
end subroutine slinlv

!> Solve linear system involving lower triangular system matrix.
!! Double precision version.
!!
!! @param[in] a lower triangular matrix.
!! @param[inout] u input RHS vector, output solution vector.
!! @author R. J. Purser
subroutine dlinlv(a,u)!                                                  [invl]
real(dp),intent(in   ) :: a(:,:)
real(dp),intent(inout) :: u(:)
integer(spi):: i
if(size(a,1) /= size(a,2) .or. size(a,1) /= size(u))&
     stop 'In dlinlv; incompatible array dimensions'
do i=1,size(u); u(i)=(u(i) - sum(u(:i-1)*a(i,:i-1)))/a(i,i); enddo
end subroutine dlinlv

!> Solve linear system involving upper triangular system matrix.
!! Single precision version.
!!
!! @param[in] a upper triangular matrix.
!! @param[inout] u input RHS vector, output solution vector.
!! @author R. J. Purser
subroutine slinuv(a,u)!                                                  [invu]
real(sp),intent(in   ) :: a(:,:)
real(sp),intent(inout) :: u(:)
integer(spi):: i
if(size(a,1) /= size(a,2) .or. size(a,1) /= size(u))&
     stop 'In linuv; incompatible array dimensions'
do i=size(u),1,-1; u(i)=(u(i) - sum(a(i+1:,i)*u(i+1:)))/a(i,i); enddo
end subroutine slinuv

!> Solve linear system involving upper triangular system matrix.
!! Double precision version.
!!
!! @param[in] a upper triangular matrix.
!! @param[inout] u input RHS vector, output solution vector.
!! @author R. J. Purser
subroutine dlinuv(a,u)!                                                  [invu]
real(dp), intent(in   ) :: a(:,:)
real(dp), intent(inout) :: u(:)
integer(spi)            :: i
if(size(a,1) /= size(a,2) .or. size(a,1) /= size(u))&
     stop 'In dlinuv; incompatible array dimensions'
do i=size(u),1,-1; u(i)=(u(i) - sum(a(i+1:,i)*u(i+1:)))/a(i,i); enddo
end subroutine dlinuv

end module pmat

