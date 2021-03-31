!> @file
!! @brief Routines dealing with the operations of banded matrices.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994/1999

!> Routines dealing with the operations of banded matrices.
!! The three special routines allow the construction of compact or
!! conventional interpolation and differencing stencils to a general
!! order of accuracy. These are:
!! - avco()  Averaging, or interpolating;
!! - dfco()  Differentiating (once);
!! - dfco2() Differentiating (twice).
!!
!! Other routines provide the tools for applying compact schemes, and for
!! the construction and application of recursive filters.
!!
!! Last modified (Purser):                              January 6th 2005
!!  added nonredundant ldltb and ltdlbv routines for symmetric matrices,
!!  and remove obsolescent routines.
!!
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994/1999
module pmat2
!============================================================================
use    pkind, only: spi,sp,dp,dpc
implicit none
private
public:: avco !< compact averaging coefficients for mid-interval interpolation
public:: dfco !< compact 1st derivative coefficients (or quadrature)
public:: dfco2 !< compact 2nd derivative coefficients
public:: clipb !< clip the unused values in a banded matrix representation
public:: cad1b !< Incorporate symmetry, and clip near end-1 of a band matrix.
public:: csb1b !< Incorporate antisymmetry, and clip near end-1 of a band matrix.
public:: cad2b !< Incorporate symmetry, and clip near end-2 of a band matrix.
public:: csb2b !< Incorporate antisymmetry, and clip near end-2 of a band matrix.
public:: ldub !< L*D*U factoring of possibly asymmetric band matrix
public:: ldltb !< L*D*L^T factors of symmetric band matrix (root-free Cholesky)
public:: udlb !< U*D*L factoring of possibly asymmetric band matrix
public:: l1ubb !< L*D*U factoring, and modification of an associated matrix
public:: l1ueb !< Special case of L1ueb adapted for compact quadratures
public:: ltdlbv !< Back-substitution, symmetric matrix, root-free Cholesky case. 
public:: udlbv !< Back-substitution for col. vector, LDU-factored matrix case.
public:: udlbx !< Like udlbv, but for x-vectors of an xy array
public:: udlby !< Like udlby, but for y-vectors of an xy array
public:: udlvb !< Like udlbv, but for row vector (instead of column vector)
public:: udlxb !< Like udlvb, but for x-vectors of an xy array
public:: udlyb !< Like udlvb, but for y-vectors of an xy array
public:: u1lbv !< Like udlbv, but for LDU representation provided by l1ubb
public:: u1lbx !< Like u1lbv, but for x-vectors of an xy array
public:: u1lby !< Like u1lbv, but for y-vectors of an xy array
public:: u1lvb !< Like u1lbv, but for row vector (instead of column vector)
public:: u1lxb !< Like u1lvb, but for x-vectors of an xy array
public:: u1lyb !< Like u1lvb, but for y-vectors of an xy array
public:: linbv !< Solve linear system with square banded matrix and a vector
public:: wrtb !< Write out, interactively, the contents of a banded matrix
real(dp),parameter:: zero=0  !< Double precision real zero

interface AVCO;   module procedure AVCO,  DAVCO,  TAVCO;   end interface
interface DFCO;   module procedure DFCO,  DDFCO,  TDFCO;   end interface
interface DFCO2;  module procedure DFCO2, DDFCO2, TDFCO2;  end interface
interface CLIPB;  module procedure clib,  clib_d, clib_c;  end interface
interface CAD1B;  module procedure CAD1B;                  end interface
interface CSB1B;  module procedure CSB1B;                  end interface
interface CAD2B;  module procedure CAD2B;                  end interface
interface CSB2B;  module procedure CSB2B;                  end interface
interface LDUB;   module procedure LDUB,  DLDUB;           end interface
interface LDLTB;  module procedure LDLTB, DLDLTB;          end interface
interface L1UBB;  module procedure L1UBB, DL1UBB;          end interface
interface L1UEB;  module procedure L1UEB, DL1UEB;          end interface
interface ltDLBV; module procedure ltdlbv,dltdlbv;         end interface
interface UDLB;   module procedure UDLB,  DUDLB;           end interface
interface UDLBV;  module procedure UDLBV, dudlbv;          end interface
interface UDLBX;  module procedure UDLBX;                  end interface
interface UDLBY;  module procedure UDLBY;                  end interface
interface UDLVB;  module procedure UDLVB;                  end interface
interface UDLXB;  module procedure UDLXB;                  end interface
interface UDLYB;  module procedure UDLYB;                  end interface
interface U1LBV;  module procedure U1LBV;                  end interface
interface U1LBX;  module procedure U1LBX;                  end interface
interface U1LBY;  module procedure U1LBY;                  end interface
interface U1LVB;  module procedure U1LVB;                  end interface
interface U1LXB;  module procedure U1LXB;                  end interface
interface U1LYB;  module procedure U1LYB;                  end interface
interface LINBV;  module procedure LINBV;                  end interface
interface WRTB;   module procedure WRTB;                   end interface
contains

!>  Compute one row of the coefficients for the compact mid-interval
!!  interpolation scheme characterized by matrix equation of the form,
!!			 A.t = B.s			       (*)
!!  Where s is the vector of "source" values, t the staggered "target" values.
!!
!! @param[in] na number of t-points operated on by this row of the A of (*)
!! @param[in] nb number of s-points operated on by this row of the B of (*)
!! @param[in] za coordinates of t-points used in this row of (*)
!! @param[in] zb coordinates of s-points used in this row of (*)
!! @param[in] z0 nominal point of application of this row of (*)
!! @param[out] a the NA coefficients A for this scheme
!! @param[out] b the NB coefficients B for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine AVCO(na,nb,za,zb,z0,a,b) !                                   [AVCO]
use pietc, only: u0,u1
use pmat,  only: inv
implicit none
integer(spi),intent(in ):: na,nb
real(sp),    intent(in ):: za(na),zb(nb),z0
real(sp),    intent(out):: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi)                    :: na1,nab,i
real(sp), dimension(na+nb,na+nb):: w
real(sp), dimension(na)         :: za0,pa
real(sp), dimension(nb)         :: zb0,pb
real(sp), dimension(na+nb)      :: ab
!=============================================================================
na1=na+1;     nab=na+nb
za0=za-z0;    zb0=zb-z0
pa=u1;        pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=2,nab; w(i,1:na)=pa;    pa=pa*za0; w(i,na1:nab)=pb; pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine AVCO 

!> Double precision version of subroutine avco for midpoint interpolation.
!!
!! @param[in] na number of target points involved in formula
!! @param[in] nb number of source points involved in formula
!! @param[in] za coordinates of target points
!! @param[in] zb coordinates of source points
!! @param[in] z0 nominal point of application of compact interpolation formula
!! @param[out] a the coefficients A for this scheme
!! @param[out] b the coefficients B for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DAVCO(na,nb,za,zb,z0,a,b) !                                  [AVCO]
use pietc, only: u0,u1
use pmat, only: inv
implicit none
integer(spi),intent(IN ):: na,nb
real(dp),    intent(IN ):: za(na),zb(nb),z0
real(dp),    intent(OUT):: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi)                   :: na1,nab,i
real(dp),dimension(na+nb,na+nb):: w
real(dp),dimension(na)         :: za0,pa
real(dp),dimension(nb)         :: zb0,pb
real(dp),dimension(na+nb)      :: ab
!=============================================================================
na1=na+1;     nab=na+nb
za0=za-z0;    zb0=zb-z0
pa=u1;        pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=2,nab; w(i,1:na)=pa;    pa=pa*za0; w(i,na1:nab)=pb; pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine DAVCO

!> Simplified computation of compact midpoint interpolation coefficients.
!!
!! @param[in] xa coordinates of target points relative to point of application
!! @param[in] xb coordinates of source points relative to point of application
!! @param[out] a the coefficients A for this scheme
!! @param[out] b the coefficients B for this scheme 
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine TAVCO(xa,xb,a,b)!                                            [AVCO]
implicit none
real(dp),dimension(:),intent(IN ):: xa,xb
real(dp),dimension(:),intent(OUT):: a,b
!-----------------------------------------------------------------------------
integer(spi):: na,nb
!=============================================================================
na=size(xa); if(na /= size(a))stop 'In tavco; sizes of a and xa different'
nb=size(xb); if(nb /= size(b))stop 'In tavco; sizes of b and xb different'
call DAVCO(na,nb,xa,xb,zero,a,b)
end subroutine TAVCO

!>  Compute one row of the coefficients for either the compact differencing or
!!  quadrature scheme characterized by matrix equation of the form,
!!			 A.d = B.c			       (*)
!!  In either case, d is the derivative (or density) of cumulative c.
!!
!! @param[in] na number of d-points operated on by this row of the A of (*)
!! @param[in] nb number of c-points operated on by this row of the B of (*)
!! @param[in] za coordinates of d-points used in this row of (*)
!! @param[in] zb coordinates of c-points used in this row of (*)
!! @param[in] z0 nominal point of application of this row of (*)
!! @param[in] A the A-coefficients for this scheme
!! @param[in] B the B-coefficients for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DFCO(na,nb,za,zb,z0,a,b)!                                    [DFCO]
use pietc_s, only: u0,u1
use pmat,    only: inv
implicit none
integer(spi),intent(IN )        :: na,nb
real(sp),    intent(IN )        :: za(na),zb(nb),z0
real(sp),    intent(OUT)        :: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi):: na1,nab,i
real(sp), dimension(na+nb,na+nb):: w
real(sp), dimension(na)         :: za0,pa
real(sp), dimension(nb)         :: zb0,pb
real(sp), dimension(na+nb)      :: ab
!=============================================================================
na1=na+1; nab=na+nb
za0=za-z0; zb0=zb-z0
pa=u1;     pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=3,nab; w(i,1:na)   =pa*(i-2); pa=pa*za0; enddo
do i=2,nab; w(i,na1:nab)=pb;       pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine DFCO 

!> Double precision version of dfco for compact differentiation coefficients
!!
!! @param[in] na number of A coefficients multiplying derivatives
!! @param[in] nb number of B coefficients multiplying cumulatives
!! @param[in] za coordinates of the density (d) points
!! @param[in] zb coordinates of the cumulative (c) points
!! @param[in] z0 nominal point of application of the compact formula
!! @param[in] a the A-coefficients for this scheme
!! @param[in] b the B-coefficients for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DDFCO(na,nb,za,zb,z0,a,b) ! Real(dp) version of              [DFCO]
use pietc, only: u0,u1
use pmat,  only: inv
implicit none
integer(spi),intent(in) :: na,nb
real(dp),    intent(in) :: za(na),zb(nb),z0
real(dp),    intent(out):: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi)                    :: na1,nab,i
real(dp), dimension(na+nb,na+nb):: w
real(dp), dimension(na)         :: za0,pa
real(dp), dimension(nb)         :: zb0,pb
real(dp), dimension(na+nb)      :: ab
!=============================================================================
na1=na+1; nab=na+nb
za0=za-z0; zb0=zb-z0
pa=u1;     pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=3,nab; w(i,1:na)   =pa*(i-2); pa=pa*za0; enddo
do i=2,nab; w(i,na1:nab)=pb;       pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine DDFCO 

!> Simplified computation of compact differencing coefficients to get 
!! derivatives d from cumulatives c, or vice-versa.
!!
!! @param[in] xa coordinates, relative to point of application, of d values
!! @param[in] xb coordinates, relatuve to point of application, of c values
!! @param[in] a the coefficients, A, for this scheme
!! @param[in] b the coefficients, B, for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine TDFCO(xa,xb,a,b)!                                            [DFCO]
implicit none
real(dp),dimension(:),intent(IN ):: xa,xb
real(dp),dimension(:),intent(OUT):: a,b
!-----------------------------------------------------------------------------
integer(spi):: na,nb
!=============================================================================
na=size(xa); if(na /= size(a))stop 'In tdfco; sizes of a and xa different'
nb=size(xb); if(nb /= size(b))stop 'In tdfco; sizes of b and xb different'
call DDFCO(na,nb,xa,xb,zero,a,b)
end subroutine TDFCO

!>  Compute one row of the coefficients for either the compact second-
!!  differencing scheme characterized by matrix equation of the form,
!!			 A.d = B.c			       (*)
!!  Where d is the second-derivative of c.
!!
!! @param[in] na number of d-points operated on by this row of the A of (*)
!! @param[in] nb number of c-points operated on by this row of the B of (*)
!! @param[in] za coordinates of d-points used in this row of (*)
!! @param[in] zb coordinates of c-points used in this row of (*)
!! @param[in] z0 nominal point of application of this row of (*)
!! @param[in] a the NA coefficients A for this scheme
!! @param[in] b the NB coefficients B for this scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DFCO2(na,nb,za,zb,z0,a,b)!                                  [DFCO2] 
use pietc_s, only: u0,u1
use pmat, only: inv
implicit none
integer(spi), intent(IN ):: na,nb
real(sp),     intent(IN ):: za(na),zb(nb),z0
real(sp),     intent(OUT):: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi)                    :: na1,nab,i
real(sp), dimension(na+nb,na+nb):: w
real(sp), dimension(na)         :: za0,pa
real(sp), dimension(nb)         :: zb0,pb
real(sp), dimension(na+nb)      :: ab
!=============================================================================
na1=na+1; nab=na+nb
za0=za-z0; zb0=zb-z0
pa=u1;     pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=4,nab; w(i,1:na)   =pa*(i-2)*(i-3); pa=pa*za0; enddo
do i=2,nab; w(i,na1:nab)=pb;             pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine DFCO2 

!> Double precision version of DFCO2 to get 2nd-derivative coefficients 
!!
!! @param[in] na number of 2nd-derivative (d) points in compact formula
!! @param[in] nb number of source points (c)
!! @param[in] za coordinates of 2nd-derivative points
!! @param[in] zb coordinates of source points
!! @param[in] z0 nominal coordinate of application
!! @param[in] a coefficients A for derivative points
!! @param[in] b coefficients B for source points
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DDFCO2(na,nb,za,zb,z0,a,b) ! Real(dp) version of            [DFCO2]
use pietc, only: u0,u1
use pmat, only: inv
implicit none
integer(spi),intent(IN )           :: na,nb
real(dp),    intent(IN )           :: za(na),zb(nb),z0
real(dp),    intent(OUT)           :: a(na),b(nb)
!-----------------------------------------------------------------------------
integer(spi)                    :: na1,nab,i
real(dp), dimension(na+nb,na+nb):: w
real(dp), dimension(na)         :: za0,pa
real(dp), dimension(nb)         :: zb0,pb
real(dp), dimension(na+nb)      :: ab
!=============================================================================
na1=na+1; nab=na+nb
za0=za-z0; zb0=zb-z0
pa=u1;     pb=-u1
w=u0;         ab=u0
w(1,1:na)=u1; ab(1)=u1
do i=4,nab; w(i,1:na)   =pa*(i-2)*(i-3); pa=pa*za0; enddo
do i=2,nab; w(i,na1:nab)=pb;             pb=pb*zb0; enddo
call INV(w,ab)
a=ab(1:na); b=ab(na1:nab)
end subroutine ddfco2 

!> Simplified computation of compact 2nd-derivative coefficients
!!
!! @param[in] xa Relative coordinates of derivatives
!! @param[in] xb Relative coordinates of source points
!! @param[out] a coefficients A for the derivatives in compact scheme
!! @param[out] b coefficients B for source values in the compact scheme
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine TDFCO2(xa,xb,a,b)!                                          [DFCO2]
!=============================================================================
real(dp),dimension(:),intent(IN ):: xa,xb
real(dp),dimension(:),intent(OUT):: a,b
!-----------------------------------------------------------------------------
integer(spi):: na,nb
!=============================================================================
na=size(xa); if(na /= size(a))stop 'In tdfco2; sizes of a and xa different'
nb=size(xb); if(nb /= size(b))stop 'In tdfco2; sizes of b and xb different'
call DDFCO2(na,nb,xa,xb,zero,a,b)
end subroutine TDFCO2

!> Clip (set to zero) the unused values in a banded matrix representation
!!
!! @param[in] m1 number of matrix rows
!! @param[in] m2 number of matrix columns
!! @param[in] mah1 number of subdiagonals
!! @param[in] mah2 number of superdiagonals
!! @param[inout] a single precision matrix elements, stored compactly as rows
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
pure subroutine CLIB(m1,m2,mah1,mah2,a)!                               [CLIPB]
use pietc_s, only: u0
implicit none
integer(spi), intent(IN   ) :: m1, m2, mah1, mah2
real(sp),     intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=u0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=u0; enddo
end subroutine CLIB

!> Clip (set to zero) the unused values in a banded matrix representation
!!
!! @param[in] m1 number of matrix rows
!! @param[in] m2 number of matrix columns
!! @param[in] mah1 number of subdiagonals
!! @param[in] mah2 number of superdiagonals
!! @param[in] a double precision matrix elements, stored compactly as rows
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
pure subroutine clib_d(m1,m2,mah1,mah2,a)!                             [CLIPB]
use pietc, only: u0
implicit none
integer(spi),intent(IN   ) :: m1, m2, mah1, mah2
real(dp),    intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=u0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=u0; enddo
end subroutine clib_d

!> Clip (set to zero) the unused values in a banded matrix representation
!!
!! @param[in] m1 number of matrix rows
!! @param[in] m2 number of matrix columns
!! @param[in] mah1 number of subdiagonals
!! @param[in] mah2 number of superdiagonals
!! @param[in] a complex matrix elements, stored compactly as rows
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
pure subroutine clib_c(m1,m2,mah1,mah2,a)!                              [CLIPB]
use pietc, only: c0
implicit none
integer(spi), intent(IN   ) :: m1, m2, mah1, mah2
complex(dpc), intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=c0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=c0; enddo
end subroutine clib_c

!> Incorporate operand symmetry and clip near end-1 of a band matrix operator.
!!
!! @param[in] m1 Size of implied full matrix
!! @param[in] mah1 Left semi-bandwidths (subdiagonals) of A.
!! @param[in] mah2 Right semi-bandwidths (superdiagonals) of A.
!! @param[in] mirror2 2*location of symmetry axis relative to end-1 operand element.
!! @param[inout] a Input operator, output as symmetrized and clipped at end-1
!!      Note: although m2 is not used here, it IS used in companion routines
!!            cad2b and csb2b; it is retained in the interests of uniformity.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine CAD1B(m1,mah1,mah2,mirror2,a)!                              [CAD1B]
use pietc_s, only: u0
implicit none
integer(spi),intent(IN   ):: m1,mah1,mah2,mirror2
real(sp),    intent(INOUT):: a(0:m1-1,-mah1:mah2)
!-----------------------------------------------------------------------------
integer(spi):: i,i2,jm,jp,jpmax
!=============================================================================
if(mirror2+mah1 > mah2)stop 'In CAD1B; mah2 insufficient'
do i=0,m1-1; i2=i*2; jpmax=mirror2+mah1-i2; if(jpmax <= -mah1)exit
   do jm=-mah1,mah2; jp=mirror2-jm-i2; if(jp <= jm)exit
      a(i,jp)=a(i,jp)+a(i,jm) ! Reflect and add
      a(i,jm)=u0              ! zero the exterior part
   enddo
enddo
end subroutine CAD1B

!> Like cad1b, but for antisymmetric operand.
!!
!! @param[in] m1 Size of implied full matrix
!! @param[in] mah1 Left semi-bandwidths (subdiagonals) of A.
!! @param[in] mah2 Right semi-bandwidths (superdiagonals) of A.
!! @param[in] mirror2 2*location of symmetry axis relative to end-1 operand element.
!! @param[inout] a Input operator, output as clipped antisymmetric at end-1.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine CSB1B(m1,mah1,mah2,mirror2,a)!                              [CSB1B]
use pietc_s, only: u0
implicit none
integer(spi),intent(IN   ):: m1,mah1,mah2,mirror2
real(sp),    intent(INOUT):: a(0:m1-1,-mah1:mah2)
!-----------------------------------------------------------------------------
integer(spi):: i,i2,jm,jp,jpmax
!=============================================================================
if(mirror2+mah1 > mah2)stop 'In CSB1B; mah2 insufficient'
do i=0,m1-1; i2=i*2; jpmax=mirror2+mah1-i2; if(jpmax < -mah1)exit
   do jm=-mah1,mah2; jp=mirror2-jm-i2; if(jp < jm)exit
      a(i,jp)=a(i,jp)-a(i,jm) ! Reflect and subtract
      a(i,jm)=u0              ! zero the exterior part
   enddo
enddo
end subroutine CSB1B

!> Incorporate symmetry and clip near end-2 of a band matrix.
!!
!! @param[in] m1 Number of rows of full matrix A
!! @param[in] m2 Number of columns of implied full matrix A
!! @param[in] mah1 Left semi-bandwidths (subdiagonals) of A.
!! @param[in] mah2 Right semi-bandwidths (superdiagonals) of A.
!! @param[in] mirror2 2*location of symmetry axis relative to end-2 operand element.
!! @param[inout] a Input operator, output as symmetrized and clipped at end-2.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine CAD2B(m1,m2,mah1,mah2,mirror2,a)!                           [CAD2B]
use pietc_s, only: u0
implicit none
integer(spi),intent(IN   ):: m1,m2,mah1,mah2,mirror2
real(sp),    intent(INOUT):: a(1-m1:0,m1-m2-mah1:m1-m2+mah2)
!-----------------------------------------------------------------------------
integer(spi):: i,i2,jm,jp,jmmin,nah1,nah2
!=============================================================================
nah1=mah1+m2-m1; nah2=mah2+m1-m2 ! Effective 2nd-index bounds of A
if(mirror2-nah1 > -nah2)stop 'In CAD2B; mah1 insufficient'
do i=0,1-m1,-1; i2=i*2; jmmin=mirror2-nah2-i2; if(jmmin >= nah2)exit
   do jp=nah2,nah1,-1; jm=mirror2-jp-i2; if(jm >= jp)exit
      a(i,jm)=a(i,jm)+a(i,jp) ! Reflect and add
      a(i,jp)=u0              ! zero the exterior part
   enddo
enddo
end subroutine CAD2B

!> Incorporate operand antisymmetry and clip near end-2 of a band matrix.
!!
!! @param[in] m1 Number of rows of matrix A
!! @param[in] m2 Number of columns of matrix A
!! @param[in] mah1 Left semi-bandwidths (subdiagonals) of A.
!! @param[in] mah2 Right semi-bandwidths (superdiagonals) of A.
!! @param[in] mirror2 2*location of symmetry axis relative to end-2 operand element.
!! @param[inout] a Input operator, output as antisymmetrized and clipped at end-2.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine CSB2B(m1,m2,mah1,mah2,mirror2,a)!                           [CSB2B]
use pietc_s, only: u0
implicit none
integer(spi),intent(IN   ):: m1,m2,mah1,mah2,mirror2
real(sp),    intent(INOUT):: a(1-m1:0,m1-m2-mah1:m1-m2+mah2)
!-----------------------------------------------------------------------------
integer(spi):: i,i2,jm,jp,jmmin,nah1,nah2
!=============================================================================
nah1=mah1+m2-m1; nah2=mah2+m1-m2 ! Effective 2nd-index bounds of A
if(mirror2-nah1 > -nah2)stop 'In CSB2B; mah1 insufficient'
do i=0,1-m1,-1; i2=i*2; jmmin=mirror2-nah2-i2; if(jmmin > nah2)exit
   do jp=nah2,nah1,-1; jm=mirror2-jp-i2; if(jm > jp)exit
      a(i,jm)=a(i,jm)-a(i,jp) ! Reflect and subtract
      a(i,jp)=u0              ! zero the exterior part
   enddo
enddo
end subroutine CSB2B

!> [L]*[D]*[U] factoring of single precision band-matrix.
!! [L] is lower triangular with unit main diagonal
!! [D] is a diagonal matrix
!! [U] is upper triangular with unit main diagonal
!! @param[in] m The number of rows of array A
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[inout] a input as the asymmetric band matrix, [A]. On output, it 
!! contains the factors in the encoding, [L-I]+[D**-1]+[U-I], I=identity.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine LDUB(m,mah1,mah2,a)!                                         [LDUB]
use pietc_s, only: u0,u1
implicit none
integer(spi),intent(IN   ):: m,mah1, mah2 
real(sp),    intent(INOUT):: a(m,-mah1:mah2) 
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jp, i
real(sp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jmost=min(m,j+mah2)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)then
      print '(" Failure in LDUB:"/" Matrix requires pivoting or is singular")'
      stop
   endif
   ajji=u1/ajj
   a(j,0)=ajji
   do i=jp,imost
      aij=ajji*a(i,j-i)
      a(i,j-i)=aij
      a(i,jp-i:jmost-i)=a(i,jp-i:jmost-i)-aij*a(j,1:jmost-j)
   enddo
   a(j,1:jmost-j)=ajji*a(j,1:jmost-j)
enddo
end subroutine LDUB

!> [L]*[D]*[U] factoring of double precision band-matrix.
!! [L] is lower triangular with unit main diagonal
!! [D] is a diagonal matrix
!! [U] is upper triangular with unit main diagonal
!! @param[in] m The number of rows of array A
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[inout] a input as the asymmetric band matrix, [A]. On output, it 
!! contains the factors in the encoding, [L-I]+[D**-1]+[U-I], I=identity.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DLDUB(m,mah1,mah2,a) ! Real(dp) version of                   [LDUB]
use pietc, only: u0,u1
implicit none
integer(spi),intent(IN   ):: m,mah1, mah2 
real(dp),    intent(INOUT):: a(m,-mah1:mah2) 
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jp, i
real(dp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
  imost=min(m,j+mah1)
  jmost=min(m,j+mah2)
  jp=j+1
  ajj=a(j,0)
  if(ajj == u0)then
    print '(" Fails in LDUB_d:"/" Matrix requires pivoting or is singular")'
    stop
  endif
  ajji=u1/ajj
  a(j,0)=ajji
  do i=jp,imost
    aij=ajji*a(i,j-i)
    a(i,j-i)=aij
    a(i,jp-i:jmost-i)=a(i,jp-i:jmost-i)-aij*a(j,1:jmost-j)
  enddo
  a(j,1:jmost-j)=ajji*a(j,1:jmost-j)
enddo
end subroutine DLDUB

!> [L]*[D]*[L^T] factoring of symmetric band matrix A (root-free Cholesky).
!! [L] is lower triangular with unit main diagonal
!! [D] is a diagonal matrix.
!!
!! @param[in] m size of symmetric matrix A
!! @param[in] mah1 semi-bandwidth of matrix A
!! @param[inout] a input lower (left) part of symmetric A; output its factors
!! encoded as [L-I]+[D**-1]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine LDLTB(m,mah1,a) ! Real(sp) version of                       [LDLTB]
use pietc_s, only: u0,u1
integer(spi),intent(IN   ):: m,mah1
real(sp),    intent(INOUT):: a(m,-mah1:0) 
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jp, i,k
real(sp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
  imost=min(m,j+mah1)
  jp=j+1
  ajj=a(j,0)
  if(ajj == u0)then
    print '(" Fails in LDLTB:"/" Matrix requires pivoting or is singular")'
    stop
  endif
  ajji=u1/ajj
  a(j,0)=ajji
  do i=jp,imost
    aij=a(i,j-i)
    a(i,j-i)=ajji*aij
    do k=jp,i
       a(i,k-i)=a(i,k-i)-aij*a(k,j-k)
    enddo
  enddo
enddo
end subroutine LDLTB

!> [L]*[D]*[L^T] factoring of symmetric matrix A (root-free Cholesky).
!! [L] is lower triangular with unit main diagonal
!! [D] is a diagonal matrix
!! @param[in] m size of symmetric matrix A
!! @param[in] mah1 semi-bandwidth of matrix A
!! @param[inout] a input lower (left) part of symmetric A; output its factors
!! encoded as [L-I]+[D**-1]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DLDLTB(m,mah1,a) ! Real(dp) version of                   [LDLTB]
use pietc, only: u0,u1
integer(spi),intent(IN   ) :: m,mah1
real(dp),    intent(INOUT) :: a(m,-mah1:0) 
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jp, i,k
real(dp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)then
      print '(" Fails in LDLTB_d:"/" Matrix requires pivoting or is singular")'
      stop
   endif
   ajji=u1/ajj
   a(j,0)=ajji
   do i=jp,imost
      aij=a(i,j-i)
      a(i,j-i)=ajji*aij
      do k=jp,i
         a(i,k-i)=a(i,k-i)-aij*a(k,j-k)
      enddo
   enddo
enddo
end subroutine DLDLTB

!> [U]*[D]*[L] factoring of single precision band matrix A
!! [U] is upper triangular with unit main diagonal
!! [D] is a diagonal matrix
!! [L] is lower triangular with unit main diagonal
!! @param[in] m number of rows of A
!! @param[in] mah1 number of subdiagonals of A
!! @param[in] mah2 number of superdiagonals of A
!! @param[inout] a Input single precision band matrix A; output its factors
!! encoded as [U-I]+[D**-1]+[L-I]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine UDLB(m,mah1,mah2,a) ! Reversed-index version of ldub        [UDLB]
implicit none
integer(spi),                    intent(IN   ) :: m,mah1,mah2
real(sp),dimension(m,-mah1:mah2),intent(INOUT) :: a(m,-mah1:mah2)
!-----------------------------------------------------------------------------
real(sp),dimension(m,-mah2:mah1):: at
!=============================================================================
at=a(m:1:-1,mah2:-mah1:-1); call LDUB(m,mah2,mah1,at)
a=at(m:1:-1,mah1:-mah2:-1)
end subroutine UDLB 

!> [U]*[D]*[L] factoring of double precision band matrix A
!! [U] is upper triangular with unit main diagonal
!! [D] is a diagonal matrix
!! [L] is lower triangular with unit main diagonal
!! @param[in] m number of rows of A
!! @param[in] mah1 number of subdiagonals of A
!! @param[in] mah2 number of superdiagonals of A
!! @param[inout] a Input double precision band matrix A; output its factors
!! encoded as [U-I]+[D**-1]+[L-I]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DUDLB(m,mah1,mah2,a) ! real(dp) version of udlb              [UDLB]
implicit none
integer(spi),                    intent(IN   ) :: m,mah1,mah2
real(dp),dimension(m,-mah1:mah2),intent(INOUT) :: a(m,-mah1:mah2)
!-----------------------------------------------------------------------------
real(dp),dimension(m,-mah2:mah1):: at
!=============================================================================
at=a(m:1:-1,mah2:-mah1:-1); call DLDUB(m,mah2,mah1,at)
a=at(m:1:-1,mah1:-mah2:-1)
end subroutine DUDLB 

!> [L]*[D]*[U] factoring of band-matrix  [A], modify [B] --> [D**-1]*[B] 
!! [L] lower triangular with unit main diagonal
!! [D] diagonal matrix
!! [U] upper triangular with unit main diagonal
!! [B] associated band matrix with same number of rows as [A] 
!!  lower triangular elements of [A] by [D**-1]*[L]*[D], the upper by [U],
!!  replace matrix [B] by [D**-1]*[B].
!!
!! @param[in] m Number of rows of A and B
!! @param[in] mah1 number of subdiagonals of A
!! @param[in] mah2 number of superdiagonals of A
!! @param[in] mbh1 number of subdiagonals of B
!! @param[in] mbh2 number of superdiagonals of B
!! @param[inout] a input as band matrix, output as lower and upper triangulars with 1s
!! implicitly assumed to lie on the main diagonal. The product of these
!! triangular matrices is [D**-1]*[A], where [D] is a diagonal matrix.
!! @param[inout] b Input single precision band matrix B; output [D**-1 B]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine L1UBB(m,mah1,mah2,mbh1,mbh2,a,b)!                           [L1UBB] 
use pietc_s, only: u0,u1
implicit none
integer(spi), intent(IN   ) ::  m,mah1, mah2, mbh1, mbh2 
real(sp),     intent(INOUT) :: a(m,-mah1:mah2), b(m,-mbh1:mbh2)
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jleast, jp, i
real(sp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jmost=min(m,j+mah2)
   jleast=max(1,j-mah1)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)stop 'In L1UBB; zero element found in diagonal factor'
   ajji=u1/ajj
   a(j,jleast-j:jmost-j) = ajji * a(j,jleast-j:jmost-j)
   do i=jp,imost
      aij=a(i,j-i)
      a(i,jp-i:jmost-i) = a(i,jp-i:jmost-i) - aij*a(j,jp-j:jmost-j)
   enddo
   a(j,0)=u1
   b(j,-mbh1:mbh2) = ajji * b(j,-mbh1:mbh2)
enddo
end subroutine L1UBB

!> Double precision version of L1UBB
!!
!! @param[in] m Number of rows of A and B
!! @param[in] mah1 left half-width of fortran array A
!! @param[in] mah2 right half-width of fortran array A
!! @param[in] mbh1 left half-width of fortran array B
!! @param[in] mbh2 left half-width of fortran array B
!! @param[inout] a Input double precision band matrix A; output factors encoded
!! as [D**-1 * L * D]+[U-I]
!! @param[inout] b Input double precision band matrix B; output [D**-1 B]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DL1UBB(m,mah1,mah2,mbh1,mbh2,a,b) ! Real(dp) version of     [L1UBB]
use pietc, only: u0,u1
implicit none
integer(spi),intent(IN   ) :: m,mah1, mah2, mbh1, mbh2 
real(dp),    intent(INOUT) :: a(m,-mah1:mah2), b(m,-mbh1:mbh2)
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jleast, jp, i
real(dp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jmost=min(m,j+mah2)
   jleast=max(1,j-mah1)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)stop 'In L1UBB_d; zero element found in diagonal factor'
   ajji=u1/ajj
   a(j,jleast-j:jmost-j) = ajji * a(j,jleast-j:jmost-j)
   do i=jp,imost
      aij=a(i,j-i)
      a(i,jp-i:jmost-i) = a(i,jp-i:jmost-i) - aij*a(j,jp-j:jmost-j)
   enddo
   a(j,0)=u1
   b(j,-mbh1:mbh2) = ajji * b(j,-mbh1:mbh2)
enddo
end subroutine DL1UBB

!>  Form the [L]*[D]*[U] decomposition of asymmetric band-matrix [A]
!!  replace all but row zero of the lower triangular elements of [A]
!!  by [D**-1]*[L]*[D], the upper by [U], replace matrix [B] by
!!  [D**-1]*[B].
!!
!!  This is a special adaptation of L1UBB used to process quadrature
!!  weights for QEDBV etc in which the initial quadrature value is
!!  provided as input instead of being implicitly assumed zero (which
!!  is the case for QZDBV etc).
!!
!! @param[in] m number of rows of B, one less than the rows of A (which has "row 0")
!! @param[in] mah1 left half-width of fortran array A
!! @param[in] mah2 right half-width of fortran array A
!! @param[in] mbh1 left half-width of fortran array B
!! @param[in] mbh2 right half-width of fortran array B
!! @param[inout] a input as band matrix, output as lower and upper triangulars with 1s
!! implicitly assumed to lie on the main diagonal. The product of these
!! triangular matrices is [D**-1]*[A], where [D] is a diagonal matrix.
!! @param[inout] b in as band matrix, out as same but premultiplied by diagonal [D**-1]
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1998
subroutine L1UEB(m,mah1,mah2,mbh1,mbh2,a,b)!                           [L1UEB]
use pietc_s, only: u0,u1
implicit none
integer(spi),intent(IN   ) :: m,mah1, mah2, mbh1, mbh2 
real(sp),    intent(INOUT) :: a(0:m,-mah1:mah2), b(m,-mbh1:mbh2)
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jleast, jp, i
real(sp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jmost=min(m,j+mah2)
   jleast=max(0,j-mah1)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)stop 'In L1UEB; zero element found in diagonal factor'
   ajji=u1/ajj
   a(j,jleast-j:jmost-j) = ajji * a(j,jleast-j:jmost-j)
   do i=jp,imost
      aij=a(i,j-i)
      a(i,jp-i:jmost-i) = a(i,jp-i:jmost-i) - aij*a(j,jp-j:jmost-j)
   enddo
   a(j,0)=u1
   b(j,-mbh1:mbh2) = ajji * b(j,-mbh1:mbh2)
enddo
end subroutine L1UEB

!> Double precision version of L1UEB
!!
!! @param[in] m number of rows of band matrices A and B
!! @param[in] mah1 number of subdiagonals of A
!! @param[in] mah2 number of superdiagonals of A
!! @param[in] mbh1 number of subdiagonals of B
!! @param[in] mbh2 number of superdiagonals of B
!! @param[inout] a Input matrix A; output encoded L, D, U, factors 
!! @param[inout] b Input matrix B; output D**-1 * B
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine DL1UEB(m,mah1,mah2,mbh1,mbh2,a,b) ! Real(dp) version of     [L1UEB]
use pietc, only: u0,u1
implicit none
integer(spi),intent(IN   ):: m,mah1, mah2, mbh1, mbh2 
real(dp),    intent(INOUT):: a(0:,-mah1:), b(:,-mbh1:)
!-----------------------------------------------------------------------------
integer(spi):: j, imost, jmost, jleast, jp, i
real(dp)    :: ajj, ajji, aij
!=============================================================================
do j=1,m
   imost=min(m,j+mah1)
   jmost=min(m,j+mah2)
   jleast=max(0,j-mah1)
   jp=j+1
   ajj=a(j,0)
   if(ajj == u0)stop 'In L1UEB_D; zero element found in diagonal factor'
   ajji=u1/ajj
   a(j,jleast-j:jmost-j) = ajji * a(j,jleast-j:jmost-j)
   do i=jp,imost
      aij=a(i,j-i)
      a(i,jp-i:jmost-i) = a(i,jp-i:jmost-i) - aij*a(j,jp-j:jmost-j)
   enddo
   a(j,0)=u1
   b(j,-mbh1:mbh2) = ajji * b(j,-mbh1:mbh2)
enddo
end subroutine DL1UEB

!>  Back-substitution step of linear inversion involving banded LDU factored 
!!  matrix and input vector, v. Output vector is [U**-1]*[D**-1]*[L**-1]*v
!!
!! @param[in] m the number of rows assumed for A and for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] a encodes the L*D*U factorization of the linear-system
!! matrix, as supplied by subroutine LDUB, in single precision
!! @param[inout] v input as right-hand-side vector, output as solution vector
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLBV(m,mah1,mah2,a,v)!                                    [UDLBV]
implicit none
integer(spi),intent(IN   ):: m, mah1, mah2
real(sp),    intent(IN   ):: a(m,-mah1:mah2)
real(sp),    intent(INOUT):: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(sp)    :: vj
!=============================================================================
do j=1,m
   vj=v(j)
   do i=j+1,min(m,j+mah1); v(i)=v(i)-a(i,j-i)*vj; enddo; v(j)=a(j,0)*vj
enddo
do j=m,2,-1
   vj=v(j)
   do i=max(1,j-mah2),j-1; v(i)=v(i)-a(i,j-i)*vj; enddo
enddo
end subroutine UDLBV

!>  Back-substitution step of linear inversion involving banded LDU factored 
!!  matrix and input vector, v. Output vector is [U**-1]*[D**-1]*[L**-1]*v
!!
!! @param[in] m the number of rows assumed for A and for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] a encodes the L*D*U factorization of the linear-system
!! matrix, as supplied by subroutine LDUB, in double precision
!! @param[inout] v input as right-hand-side vector, output as solution vector
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine dudlbv(m,mah1,mah2,a,v)!                                    [udlbv]
implicit none
integer(spi),intent(IN   ) :: m, mah1, mah2
real(dp),    intent(IN   ) :: a(m,-mah1:mah2)
real(dp),    intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(dp)    :: vj
!=============================================================================
do j=1,m
   vj=v(j)
   do i=j+1,min(m,j+mah1); v(i)=v(i)-a(i,j-i)*vj; enddo; v(j)=a(j,0)*vj
enddo
do j=m,2,-1
   vj=v(j)
   do i=max(1,j-mah2),j-1; v(i)=v(i)-a(i,j-i)*vj; enddo
enddo
end subroutine dudlbv

!> Like udlbv, except assuming a is the ldlt decomposition of a
!! SYMMETRIC banded matrix, with only the non-upper part provided (to
!! avoid redundancy)
!!
!! @param[in] m the number of rows assumed for A and for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] a encodes the L*D*L^T factorization of the linear-system
!! matrix, as supplied by subroutine LDLTB
!! @param[inout] v input as right-hand-side vector, output as solution vector
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine ltdlbv(m,mah1,a,v)!                                        [ltdlbv]
implicit none
integer(spi),intent(IN   ) :: m, mah1
real(sp),    intent(IN   ) :: a(m,-mah1:0)
real(sp),    intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(sp)    :: vj
!=============================================================================
do j=1,m
   vj=v(j)
   do i=j+1,min(m,j+mah1); v(i)=v(i)-a(i,j-i)*vj; enddo; v(j)=a(j,0)*vj
enddo
do j=m,2,-1
   vj=v(j)
   do i=max(1,j-mah1),j-1; v(i)=v(i)-a(j,i-j)*vj; enddo
enddo
end subroutine ltdlbv


!> Like udlbv, except assuming a is the ltdl decomposition of a
!! SYMMETRIC banded matrix, with only the non-upper part provided (to
!! avoid redundancy)
!!
!! @param[in] m the number of rows assumed for A and for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[inout] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!! matrix, as supplied by subroutine LDUB
!! @param[inout] v input as right-hand-side vector, output as solution vector
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine dltdlbv(m,mah1,a,v)!                                       [ltdlbv]
implicit none
integer(spi),intent(IN   ) :: m, mah1
real(dp),    intent(IN   ) :: a(m,-mah1:0)
real(dp),    intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(dp)    :: vj
!=============================================================================
do j=1,m
   vj=v(j)
   do i=j+1,min(m,j+mah1); v(i)=v(i)-a(i,j-i)*vj; enddo; v(j)=a(j,0)*vj
enddo
do j=m,2,-1
   vj=v(j)
   do i=max(1,j-mah1),j-1; v(i)=v(i)-a(j,i-j)*vj; enddo
enddo
end subroutine dltdlbv

!>  Back-substitution step of parallel linear inversion involving
!!  Banded matrix and X-Vectors.
!!
!! @param[in] mx the number of rows assumed for A and length of
!!     X-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] my number of parallel X-vectors inverted
!! @param[in] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!!     matrix, as supplied by subroutine LDUB or, if N=NA, by LDUB
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLBX(mx,mah1,mah2,my,a,v)!                                [UDLBX]
implicit none
integer(spi),intent(IN   ) :: mx, mah1, mah2, my
real(sp),    intent(IN   ) :: a(mx,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: jx, ix
!=============================================================================
do jx=1,mx
   do ix=jx+1,min(mx,jx+mah1); v(ix,:) = v(ix,:) - a(ix,jx-ix)*v(jx,:); enddo
   v(jx,:) = a(jx,0) * v(jx,:)
enddo
do jx=mx,2,-1
   do ix=max(1,jx-mah2),jx-1; v(ix,:) = v(ix,:) - a(ix,jx-ix)*v(jx,:); enddo
enddo
end subroutine UDLBX

!>  Back-substitution step of parallel linear inversion involving
!!  Banded matrix and Y-Vectors.
!!
!! @param[in] my the number of rows assumed for A and length of
!!     Y-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] mx number of parallel Y-vectors inverted
!! @param[in] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!!     matrix, as supplied by subroutine LDUB or, if N=NA, by LDUB
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLBY(my,mah1,mah2,mx,a,v)!                                [UDLBY]
implicit none
integer(spi),intent(IN   ) :: my, mah1, mah2, mx
real(sp),    intent(IN   ) :: a(my,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: iy, jy
!=============================================================================
do jy=1,my
   do iy=jy+1,min(my,jy+mah1); v(:,iy) = v(:,iy)-a(iy,jy-iy)*v(:,jy); enddo
   v(:,jy)=a(jy,0)*v(:,jy)
enddo
do jy=my,2,-1
   do iy=max(1,jy-mah2),jy-1; v(:,iy)=v(:,iy)-a(iy,jy-iy)*v(:,jy); enddo
enddo
end subroutine UDLBY

!>  Back-substitution step of linear inversion involving
!!  row-Vector and Banded matrix.
!!
!! @param[in] m the number of rows assumed for A and columns for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[inout] v input as right-hand-side row-vector, output as solution vector
!! @param[in] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!!     matrix, as supplied by subroutine LDUB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLVB(m,mah1,mah2,v,a)!                                    [UDLVB]
implicit none
integer(spi), intent(IN   ) :: m, mah1, mah2
real(sp),     intent(IN   ) :: a(m,-mah1:mah2)
real(sp),     intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(sp)    :: vi
!=============================================================================
do i=1,m
   vi=v(i)
   do j=i+1,min(m,i+mah2); v(j)=v(j)-vi*a(i,j-i); enddo
   v(i)=vi*a(i,0)
enddo
do i=m,2,-1
   vi=v(i)
   do j=max(1,i-mah1),i-1; v(j)=v(j)-vi*a(i,j-i); enddo
enddo
end subroutine UDLVB

!>  Back-substitution step of parallel linear inversion involving
!!  Banded matrix and row-X-Vectors.
!!
!! @param[in] mx the number of rows assumed for A and length of
!!     X-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] my number of parallel X-vectors inverted
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @param[in] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!!     matrix, as supplied by subroutine LDUB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLXB(mx,mah1,mah2,my,v,a)!                                [UDLXB]
implicit none
integer(spi),intent(IN   ) :: mx, mah1, mah2, my
real(sp),    intent(IN   ) :: a(mx,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: ix, jx
!=============================================================================
do ix=1,mx
   do jx=ix+1,min(mx,ix+mah2); v(jx,:)=v(jx,:)-v(ix,:)*a(ix,jx-ix); enddo
   v(ix,:)=v(ix,:)*a(ix,0)
enddo
do ix=mx,2,-1
   do jx=max(1,ix-mah1),ix-1; v(jx,:)=v(jx,:)-v(ix,:)*a(ix,jx-ix); enddo
enddo
end subroutine UDLXB

!>  BACk-substitution step of parallel linear inversion involving
!!  Banded matrix and row-Y-Vectors.
!!
!! @param[in] my the number of rows assumed for A and length of
!!     Y-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] mx number of parallel Y-vectors inverted
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @param[in] a encodes the (L)*(D**-1)*(U) factorization of the linear-system
!!     matrix, as supplied by subroutine LDUB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine UDLYB(my,mah1,mah2,mx,v,a)!                                [UDLYB]
implicit none
integer(spi),intent(IN   ) :: my, mah1, mah2, mx
real(sp),    intent(IN   ) :: a(my,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: iy, jy
!=============================================================================
do iy=1,my
   do jy=iy+1,min(my,iy+mah2); v(:,jy)=v(:,jy)-v(:,iy)*a(iy,jy-iy); enddo
   v(:,iy)=v(:,iy)*a(iy,0)
enddo
do iy=my,2,-1
   do jy=max(1,iy-mah1),iy-1; v(:,jy)=v(:,jy)-v(:,iy)*a(iy,jy-iy); enddo
enddo
end subroutine UDLYB

!>  Back-substitution step ((U**-1)*(L**-1)) of linear inversion involving
!!  special Banded matrix and right-Vector.
!!
!! @param[in] m the number of rows assumed for A and for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] a encodes the [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @param[inout] v input as right-hand-side vector, output as solution vector
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LBV(m,mah1,mah2,a,v)!                                    [U1LBV]
implicit none
integer(spi),intent(IN   ) :: m, mah1, mah2
real(sp),    intent(IN   ) :: a(m,-mah1:mah2)
real(sp),    intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(sp)    :: vj
!=============================================================================
do j=1,m
   vj=v(j)
   do i=j+1,min(m,j+mah1); v(i)=v(i)-a(i,j-i)*vj; enddo
enddo
do j=m,2,-1
   vj=v(j)
   do i=max(1,j-mah2),j-1; v(i)=v(i)-a(i,j-i)*vj; enddo
enddo
end subroutine U1LBV

!>  Special back-substitution step of parallel linear inversion involving
!!  Banded matrix and X-right-Vectors.
!!
!! @param[in] mx the number of rows assumed for A and length of
!!     X-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] my number of parallel X-vectors inverted
!! @param[in] a encodes the [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LBX(mx,mah1,mah2,my,a,v)!                                [U1LBX]
implicit none
integer(spi),intent(IN   ) :: mx, mah1, mah2, my
real(sp),    intent(IN   ) :: a(mx,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: ix, jx
!=============================================================================
do jx=1,mx
   do ix=jx+1,min(mx,jx+mah1); v(ix,:)=v(ix,:)-a(ix,jx-ix)*v(jx,:); enddo
enddo
do jx=mx,2,-1
   do ix=max(1,jx-mah2),jx-1; v(ix,:)=v(ix,:)-a(ix,jx-ix)*v(jx,:); enddo
enddo
end subroutine U1LBX

!>  Special Back-substitution step of parallel linear inversion involving
!!  Banded matrix and Y-right-Vectors.
!!
!! @param[in] my the number of rows assumed for A and length of
!!     Y-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] mx number of parallel Y-vectors inverted
!! @param[in] a encodes the [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LBY(my,mah1,mah2,mx,a,v)!                                [U1LBY]
implicit none
integer(spi),intent(IN   ) :: my, mah1, mah2, mx
real(sp),    intent(IN   ) :: a(my,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: iy, jy
!=============================================================================
do jy=1,my
   do iy=jy+1,min(my,jy+mah1); v(:,iy)=v(:,iy)-a(iy,jy-iy)*v(:,jy); enddo
enddo
do jy=my,2,-1
   do iy=max(1,jy-mah2),jy-1; v(:,iy)=v(:,iy)-a(iy,jy-iy)*v(:,jy); enddo
enddo
end subroutine U1LBY

!>  Special Back-substitution step of linear inversion involving
!!  left-Vector and Banded matrix.
!!
!! @param[in] m the number of rows assumed for A and columns for V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[inout] v input as right-hand-side row-vector, output as solution vector
!! @param[in] a encodes the special [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LVB(m,mah1,mah2,v,a)!                                    [U1LVB]
implicit none
integer(spi),intent(IN   ) :: m, mah1, mah2
real(sp),    intent(IN   ) :: a(m,-mah1:mah2)
real(sp),    intent(INOUT) :: v(m)
!-----------------------------------------------------------------------------
integer(spi):: i, j
real(sp)    :: vi
!=============================================================================
do i=1,m
   vi=v(i)
   do j=i+1,min(m,i+mah2); v(j)=v(j)-vi*a(i,j-i); enddo
enddo
do i=m,2,-1
   vi=v(i)
   do j=max(1,i-mah1),i-1; v(j)=v(j)-vi*a(i,j-i); enddo
enddo
end subroutine U1LVB

!>  Special Back-substitution step of parallel linear inversion involving
!!  Banded matrix and X-left-Vectors.
!!
!! @param[in] mx the number of rows assumed for A and length of
!!     X-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] my number of parallel X-vectors inverted
!! @param[in] v input as right-hand-side vectors, output as solution vectors
!! @param[in] a encodes the special [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LXB(mx,mah1,mah2,my,v,a)!                                [U1LXB]
implicit none
integer(spi),intent(IN   ) :: mx, mah1, mah2, my
real(sp),    intent(IN   ) :: a(mx,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: ix, jx
!=============================================================================
do ix=1,mx
   do jx=ix+1,min(mx,ix+mah2); v(jx,:)=v(jx,:)-v(ix,:)*a(ix,jx-ix); enddo
enddo
do ix=mx,2,-1
   do jx=max(1,ix-mah1),ix-1;  v(jx,:)=v(jx,:)-v(ix,:)*a(ix,jx-ix); enddo
enddo
end subroutine U1LXB

!>  Special Back-substitution step of parallel linear inversion involving
!!  special Banded matrix and Y-left-Vectors.
!!
!! @param[in] my the number of rows assumed for A and length of
!!     Y-vectors stored in V
!! @param[in] mah1 the left half-bandwidth of fortran array A
!! @param[in] mah2 the right half-bandwidth of fortran array A
!! @param[in] mx number of parallel Y-vectors inverted
!! @param[inout] v input as right-hand-side vectors, output as solution vectors
!! @param[in] a encodes the [L]*[U] factorization of the linear-system
!!     matrix, as supplied by subroutine L1UBB
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1996
subroutine U1LYB(my,mah1,mah2,mx,v,a)!                                [U1LYB]
implicit none
integer(spi),intent(IN   ) :: my, mah1, mah2, mx
real(sp),    intent(IN   ) :: a(my,-mah1:mah2)
real(sp),    intent(INOUT) :: v(mx,my)
!-----------------------------------------------------------------------------
integer(spi):: iy, jy
!=============================================================================
do iy=1,my
   do jy=iy+1,min(my,iy+mah2); v(:,jy)=v(:,jy)-v(:,iy)*a(iy,jy-iy); enddo
enddo
do iy=my,2,-1
   do jy=max(1,iy-mah1),iy-1;  v(:,jy)=v(:,jy)-v(:,iy)*a(iy,jy-iy); enddo
enddo
end subroutine U1LYB

!> Solve LINear system with square Banded-matrix and vector V.
!!
!! @param[in] m order of matrix A
!! @param[in] mah1 left half-bandwidth of A
!! @param[in] mah2 right half-bandwidth of A
!! @param[inout] a system matrix on input, its [L]*[D**-1]*[U] factorization on exit
!! @param[inout] v vector of right-hand-sides on input, solution vector on exit
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1994
subroutine LINBV(m,mah1,mah2,a,v)!                                     [LINBV]
implicit none
integer(spi),intent(IN   ) :: m, mah1, mah2
real(sp),    intent(INOUT) :: a(m,-mah1:mah2), v(m)
!=============================================================================
call LDUB(m,mah1,mah2,a)
call UDLBV(m,mah1,mah2,a,v)
end subroutine LINBV

!> Convenient routine for interactively writing out the real contents of a band 
!! matrix
!!
!! @param[in] m1 number of rows of full matrix
!! @param[in] m2 number of columns of full matrix
!! @param[in] mah1 number of sub-diagonals
!! @param[in] mah2 number of super-diagonals
!! @param[in] a contents of single precision real band matrix.
!! @author R. J. Purser, Tsukasa Fujita (JMA) @date 1999
subroutine WRTB(m1,m2,mah1,mah2,a)!                                     [WRTB]
implicit none
integer(spi),intent(IN) :: m1, m2, mah1, mah2
real(sp),    intent(IN) :: a(m1,-mah1:mah2)
!-----------------------------------------------------------------------------
integer(spi):: i1, i2, i, j1, j2, j, nj1
!=============================================================================
do i1=1,m1,20
   i2=min(i1+19,m1)
   print '(7x,6(i2,10x))',(j,j=-mah1,mah2)
   do i=i1,i2
      j1=max(-mah1,1-i)
      j2=min(mah2,m2-i)
      nj1=j1+mah1
      if(nj1==0)print '(1x,i3,6(1x,e12.5))',    i,(a(i,j),j=j1,j2)
      if(nj1==1)print '(1x,i3,12x,5(1x,e12.5))',i,(a(i,j),j=j1,j2)
      if(nj1==2)print '(1x,i3,24x,4(1x,e12.5))',i,(a(i,j),j=j1,j2)
      if(nj1==3)print '(1x,i3,36x,3(1x,e12.5))',i,(a(i,j),j=j1,j2)
      if(nj1==4)print '(1x,i3,48x,2(1x,e12.5))',i,(a(i,j),j=j1,j2)
      if(nj1==5)print '(1x,i3,60x,1(1x,e12.5))',i,(a(i,j),j=j1,j2)
   enddo
   read(*,*)
enddo
end subroutine WRTB

end module pmat2
