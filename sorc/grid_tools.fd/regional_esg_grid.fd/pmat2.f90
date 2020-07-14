!
!                                **********************************************
!                                *             MODULE pmat2                   *
!                                *  R. J. Purser, NOAA/NCEP/EMC  1994/1999    *
!                                *  jim.purser@noaa.gov                       *
!                                *  Tsukasa Fujita (JMA)             1999     *
!                                *                                            *
!                                **********************************************
!
! Routines dealing with the operations of banded matrices
! The three special routines allow the construction of compact or
! conventional interpolation and differencing stencils to a general
! order of accuracy. These are:
! AVCO:  Averaging, or interpolating;
! DFCO:  Differentiating (once);
! DFCO2: Differentiating (twice).
!
! Other routines provide the tools for applying compact schemes, and for
! the construction and application of recursive filters.
!
! Programmers:  R. J. Purser and T. Fujita
!               National Centers for Environmental Prediction.
! Last modified (Purser):                              January 6th 2005
!  added nonredundant ldltb and ltdlbv routines for symmetric matrices,
!  and remove obsolescent routines.
!                                                      January 6rd 2014
!
! DIRECT DEPENDENCIES
! Libraries[their modules]: pmat[pmat]
! Additional Modules      : pkind
!
!=============================================================================
module pmat2
!============================================================================
use    pkind, only: spi,sp,dp,dpc
implicit none
private
public:: avco,dfco,dfco2, clipb,cad1b,csb1b,cad2b,csb2b,                    &
         ldub,ldltb,udlb,l1ubb,l1ueb,ltdlbv,                                &
         udlbv,udlbx,udlby,udlvb,udlxb,udlyb,u1lbv,u1lbx,u1lby,u1lvb,u1lxb, &
         u1lyb,linbv,wrtb
real(dp),parameter:: zero=0

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

!=============================================================================
subroutine AVCO(na,nb,za,zb,z0,a,b) !                                   [AVCO]
!=============================================================================
!		    SUBROUTINE AVCO
!   R.J.Purser, National Centers for Environmental Prediction, Washington D.C.
!   jim.purser@noaa.gov					      1999
!
!  Compute one row of the coefficients for the compact mid-interval
!  interpolation scheme characterized by matrix equation of the form,
!			 A.t = B.s			       (*)
!  Where s is the vector of "source" values, t the staggered "target" values.
!
! --> NA:   number of t-points operated on by this row of the A of (*)
! --> NB:   number of s-points operated on by this row of the B of (*)
! --> ZA:   coordinates of t-points used in this row of (*)
! --> ZB:   coordinates of s-points used in this row of (*)
! --> Z0:   nominal point of application of this row of (*)
! <-- A:    the NA coefficients A for this scheme
! <-- B:    the NB coefficients B for this scheme
!=============================================================================
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
!=============================================================================
subroutine DAVCO(na,nb,za,zb,z0,a,b) !                                  [AVCO]
!=============================================================================
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
!=============================================================================
subroutine TAVCO(xa,xb,a,b)!                                            [AVCO]
!=============================================================================
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

!=============================================================================
subroutine DFCO(na,nb,za,zb,z0,a,b)!                                    [DFCO]
!=============================================================================
!   R.J.Purser, National Centers for Environmental Prediction, Washington D.C.
!   jim.purser@noaa.gov					      1999
!		    SUBROUTINE DFCO
!
!  Compute one row of the coefficients for either the compact differencing or
!  quadrature scheme characterized by matrix equation of the form,
!			 A.d = B.c			       (*)
!  In either case, d is the derivative of c.
!
! --> NA:   number of d-points operated on by this row of the A of (*)
! --> NB:   number of c-points operated on by this row of the B of (*)
! --> ZA:   coordinates of d-points used in this row of (*)
! --> ZB:   coordinates of c-points used in this row of (*)
! --> Z0:   nominal point of application of this row of (*)
! <-- A:    the A-coefficients for this scheme
! <-- B:    the B-coefficients for this scheme
!=============================================================================
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
!=============================================================================
subroutine DDFCO(na,nb,za,zb,z0,a,b) ! Real(dp) version of              [DFCO]
!=============================================================================
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
!=============================================================================
subroutine TDFCO(xa,xb,a,b)!                                            [DFCO]
!=============================================================================
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

!=============================================================================
subroutine DFCO2(na,nb,za,zb,z0,a,b)!                                  [DFCO2] 
!=============================================================================
!		    SUBROUTINE DFCO2
!   R.J.Purser, National Centers for Environmental Prediction, Washington D.C.
!   jim.purser@noaa.gov					      1999
!
!  Compute one row of the coefficients for either the compact second-
!  differencing scheme characterized by matrix equation of the form,
!			 A.d = B.c			       (*)
!  Where d is the second-derivative of c.
!
! --> NA:   number of d-points operated on by this row of the A of (*)
! --> NB:   number of c-points operated on by this row of the B of (*)
! --> ZA:   coordinates of d-points used in this row of (*)
! --> ZB:   coordinates of c-points used in this row of (*)
! --> Z0:   nominal point of application of this row of (*)
! <-- A:    the NA coefficients A for this scheme
! <-- B:    the NB coefficients B for this scheme
!=============================================================================
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
!=============================================================================
subroutine DDFCO2(na,nb,za,zb,z0,a,b) ! Real(dp) version of            [DFCO2]
!=============================================================================
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
!=============================================================================
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


!=============================================================================
pure subroutine CLIB(m1,m2,mah1,mah2,a)!                               [CLIPB]
!=============================================================================
use pietc_s, only: u0
implicit none
integer(spi), intent(IN   ) :: m1, m2, mah1, mah2
real(sp),     intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=u0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=u0; enddo
end subroutine CLIB
!=============================================================================
pure subroutine clib_d(m1,m2,mah1,mah2,a)!                             [CLIPB]
!=============================================================================
use pietc, only: u0
implicit none
integer(spi),intent(IN   ) :: m1, m2, mah1, mah2
real(dp),    intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=u0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=u0; enddo
end subroutine clib_d
!=============================================================================
pure subroutine clib_c(m1,m2,mah1,mah2,a)!                              [CLIPB]
!=============================================================================
use pietc, only: c0
implicit none
integer(spi), intent(IN   ) :: m1, m2, mah1, mah2
complex(dpc), intent(INOUT) :: a(m1,-mah1:mah2)
integer(spi):: j
do j=1,mah1; a(1:min(m1,j),-j)=c0; enddo
do j=m2-m1+1,mah2; a(max(1,m2-j+1):m1,j)=c0; enddo
end subroutine clib_c

!=============================================================================
subroutine CAD1B(m1,mah1,mah2,mirror2,a)!                              [CAD1B]
!=============================================================================
! Incorporate operand symmetry near end-1 of a band matrix operator
!
! <-> A:      Input as unclipped operator, output as symmetrized and clipped.
! m1, m2:     Sizes of implied full matrix
! mah1, mah2: Left and right semi-bandwidths of A.
! mirror2:    2*location of symmetry axis relative to end-1 operand element.
!      Note: although m2 is not used here, it IS used in companion routines
!            cad2b and csb2b; it is retained in the interests of uniformity.
!=============================================================================
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

!=============================================================================
subroutine CSB1B(m1,mah1,mah2,mirror2,a)!                              [CSB1B]
!=============================================================================
! Like cad1b, but for antisymmetric operand
!=============================================================================
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

!=============================================================================
subroutine CAD2B(m1,m2,mah1,mah2,mirror2,a)!                           [CAD2B]
!=============================================================================
! Incorporate operand symmetry near end-2 of a band matrix operator
!
! <-> A:      Input as unclipped operator, output as symmetrized and clipped.
! m1, m2:     Sizes of implied full matrix
! mah1, mah2: Left and right semi-bandwidths of A.
! mirror2:    2*location of symmetry axis relative to end-2 operand element.
!=============================================================================
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

!=============================================================================
subroutine CSB2B(m1,m2,mah1,mah2,mirror2,a)!                           [CSB2B]
!=============================================================================
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

!=============================================================================
subroutine LDUB(m,mah1,mah2,a)!                                         [LDUB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE LDUB
!  Compute [L]*[D**-1]*[U] decomposition of asymmetric band-matrix
!
! <-> A: input as the asymmetric band matrix. On output, it contains
!     the [L]*[D**-1]*[U] factorization of the input matrix, where
!     [L] is lower triangular with unit main diagonal
!     [D] is a diagonal matrix
!     [U] is upper triangular with unit main diagonal
! --> M:    The number of rows of array A
! --> MAH1: the left half-bandwidth of fortran array A
! --> MAH2: the right half-bandwidth of fortran array A
!=============================================================================
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
!=============================================================================
subroutine DLDUB(m,mah1,mah2,a) ! Real(dp) version of                   [LDUB]
!=============================================================================
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

!=============================================================================
subroutine LDLTB(m,mah1,a) ! Real(sp) version of                       [LDLTB]
!=============================================================================
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
!=============================================================================
subroutine DLDLTB(m,mah1,a) ! Real(dp) version of                   [LDLTB]
!=============================================================================
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

!=============================================================================
subroutine UDLB(m,mah1,mah2,a) ! Reversed-index version of ldub        [UDLB]
!=============================================================================
implicit none
integer(spi),                    intent(IN   ) :: m,mah1,mah2
real(sp),dimension(m,-mah1:mah2),intent(INOUT) :: a(m,-mah1:mah2)
!-----------------------------------------------------------------------------
real(sp),dimension(m,-mah2:mah1):: at
!=============================================================================
at=a(m:1:-1,mah2:-mah1:-1); call LDUB(m,mah2,mah1,at)
a=at(m:1:-1,mah1:-mah2:-1)
end subroutine UDLB 
!=============================================================================
subroutine DUDLB(m,mah1,mah2,a) ! real(dp) version of udlb              [UDLB]
!=============================================================================
implicit none
integer(spi),                    intent(IN   ) :: m,mah1,mah2
real(dp),dimension(m,-mah1:mah2),intent(INOUT) :: a(m,-mah1:mah2)
!-----------------------------------------------------------------------------
real(dp),dimension(m,-mah2:mah1):: at
!=============================================================================
at=a(m:1:-1,mah2:-mah1:-1); call DLDUB(m,mah2,mah1,at)
a=at(m:1:-1,mah1:-mah2:-1)
end subroutine DUDLB 

!=============================================================================
subroutine L1UBB(m,mah1,mah2,mbh1,mbh2,a,b)!                           [L1UBB] 
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE L1UBB
!  Form the [L]*[D]*[U] decomposition of asymmetric band-matrix  [A] replace
!  lower triangular elements of [A] by [D**-1]*[L]*[D], the upper by [U],
!  replace matrix [B] by [D**-1]*[B].
!
! <-> A input as band matrix, output as lower and upper triangulars with 1s
!     implicitly assumed to lie on the main diagonal. The product of these
!     triangular matrices is [D**-1]*[A], where [D] is a diagonal matrix.
! <-> B in as band matrix, out as same but premultiplied by diagonal [D**-1]
! --> M    Number of rows of A and B
! --> MAH1 left half-width of fortran array A
! --> MAH2 right half-width of fortran array A
! --> MBH1 left half-width of fortran array B
! --> MBH2 right half-width of fortran array B
!=============================================================================
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
!=============================================================================
subroutine DL1UBB(m,mah1,mah2,mbh1,mbh2,a,b) ! Real(dp) version of     [L1UBB]
!=============================================================================
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

!=============================================================================
subroutine L1UEB(m,mah1,mah2,mbh1,mbh2,a,b)!                           [L1UEB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1998
!		    SUBROUTINE L1UEB
!  Form the [L]*[D]*[U] decomposition of asymmetric band-matrix  [A] replace
!  all but row zero of the
!  lower triangular elements of [A] by [D**-1]*[L]*[D], the upper by [U],
!  replace matrix [B] by [D**-1]*[B].
!  This is a special adaptation of L1UBB used to process quadarature weights
!  for QEDBV etc in which the initial quadrature value is provided as input
!  instead of being implicitly assumed zero (which is the case for QZDBV etc).
!
! <-> A input as band matrix, output as lower and upper triangulars with 1s
!     implicitly assumed to lie on the main diagonal. The product of these
!     triangular matrices is [D**-1]*[A], where [D] is a diagonal matrix.
! <-> B in as band matrix, out as same but premultiplied by diagonal [D**-1]
! --> M    number of rows of B, one less than the rows of A (which has "row 0")
! --> MAH1 left half-width of fortran array A
! --> MAH2 right half-width of fortran array A
! --> MBH1 left half-width of fortran array B
! --> MBH2 right half-width of fortran array B
!=============================================================================
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
!=============================================================================
subroutine DL1UEB(m,mah1,mah2,mbh1,mbh2,a,b) ! Real(dp) version of     [L1UEB]
!=============================================================================
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

!=============================================================================
subroutine UDLBV(m,mah1,mah2,a,v)!                                    [UDLBV]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLBV
!  BACk-substitution step of linear inversion involving
!  Banded matrix and Vector.
!
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB
! <-> V input as right-hand-side vector, output as solution vector
! --> M the number of rows assumed for A and for V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
!=============================================================================
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
!=============================================================================
subroutine dudlbv(m,mah1,mah2,a,v)!                                    [udlbv]
!=============================================================================
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

!=============================================================================
subroutine ltdlbv(m,mah1,a,v)!                                        [ltdlbv]
!=============================================================================
! Like udlbv, except assuming a is the ltdl decomposition of a SYMMETRIC
! banded matrix, with only the non-upper part provided (to avoid redundancy)
!=============================================================================
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
!=============================================================================
subroutine dltdlbv(m,mah1,a,v)!                                       [ltdlbv]
!=============================================================================
! Like udlbv, except assuming a is the ltdl decomposition of a SYMMETRIC
! banded matrix, with only the non-upper part provided (to avoid redundancy)
!=============================================================================
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

!=============================================================================
subroutine UDLBX(mx,mah1,mah2,my,a,v)!                                [UDLBX]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLBX
!  BACk-substitution step of parallel linear inversion involving
!  Banded matrix and X-Vectors.
!
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB or, if N=NA, by LDUB
! <-> V input as right-hand-side vectors, output as solution vectors
! --> MX the number of rows assumed for A and length of
!     X-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MY number of parallel X-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine UDLBY(my,mah1,mah2,mx,a,v)!                                [UDLBY]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLBY
!  BACk-substitution step of parallel linear inversion involving
!  Banded matrix and Y-Vectors.
!
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB or, if N=NA, by LDUB
! <-> V input as right-hand-side vectors, output as solution vectors
! --> MY the number of rows assumed for A and length of
!     Y-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MX number of parallel Y-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine UDLVB(m,mah1,mah2,v,a)!                                    [UDLVB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLVB
!  BACk-substitution step of linear inversion involving
!  row-Vector and Banded matrix.
!
! <-> V input as right-hand-side row-vector, output as solution vector
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB
! --> M the number of rows assumed for A and columns for V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
!=============================================================================
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

!=============================================================================
subroutine UDLXB(mx,mah1,mah2,my,v,a)!                                [UDLXB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLXB
!  BACk-substitution step of parallel linear inversion involving
!  Banded matrix and row-X-Vectors.
!
! <-> V input as right-hand-side vectors, output as solution vectors
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB
! --> MX the number of rows assumed for A and length of
!     X-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MY number of parallel X-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine UDLYB(my,mah1,mah2,mx,v,a)!                                [UDLYB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE UDLYB
!  BACk-substitution step of parallel linear inversion involving
!  Banded matrix and row-Y-Vectors.
!
! <-> V input as right-hand-side vectors, output as solution vectors
! --> A encodes the (L)*(D**-1)*(U) factorization of the linear-system
!     matrix, as supplied by subroutine LDUB
! --> MY the number of rows assumed for A and length of
!     Y-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MX number of parallel Y-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine U1LBV(m,mah1,mah2,a,v)!                                    [U1LBV]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LBV
!  BACk-substitution step ((U**-1)*(L**-1)) of linear inversion involving
!  special Banded matrix and right-Vector.
!
! --> A encodes the [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! <-> V input as right-hand-side vector, output as solution vector
! --> M the number of rows assumed for A and for V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
!=============================================================================
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

!=============================================================================
subroutine U1LBX(mx,mah1,mah2,my,a,v)!                                [U1LBX]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LBX
!  Special BaCk-substitution step of parallel linear inversion involving
!  Banded matrix and X-right-Vectors.
!
! --> A encodes the [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! <-> V input as right-hand-side vectors, output as solution vectors
! --> MX the number of rows assumed for A and length of
!     X-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MY number of parallel X-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine U1LBY(my,mah1,mah2,mx,a,v)!                                [U1LBY]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LBY
!  Special BaCk-substitution step of parallel linear inversion involving
!  Banded matrix and Y-right-Vectors.
!
! --> A encodes the [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! <-> V input as right-hand-side vectors, output as solution vectors
! --> MY the number of rows assumed for A and length of
!     Y-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MX number of parallel Y-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine U1LVB(m,mah1,mah2,v,a)!                                    [U1LVB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LVB
!  Special BaCk-substitution step of linear inversion involving
!  left-Vector and Banded matrix.
!
! <-> V input as right-hand-side row-vector, output as solution vector
! --> A encodes the special [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! --> M the number of rows assumed for A and columns for V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
!=============================================================================
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

!=============================================================================
subroutine U1LXB(mx,mah1,mah2,my,v,a)!                                [U1LXB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LXB
!  Special BaCk-substitution step of parallel linear inversion involving
!  Banded matrix and X-left-Vectors.
!
! <-> V input as right-hand-side vectors, output as solution vectors
! --> A encodes the special [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! --> MX the number of rows assumed for A and length of
!     X-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MY number of parallel X-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine U1LYB(my,mah1,mah2,mx,v,a)!                                [U1LYB]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1996
!		    SUBROUTINE U1LYB
!  Special BaCk-substitution step of parallel linear inversion involving
!  special Banded matrix and Y-left-Vectors.
!
! <-> V input as right-hand-side vectors, output as solution vectors
! --> A encodes the [L]*[U] factorization of the linear-system
!     matrix, as supplied by subroutine L1UBB
! --> MY the number of rows assumed for A and length of
!     Y-vectors stored in V
! --> MAH1 the left half-bandwidth of fortran array A
! --> MAH2 the right half-bandwidth of fortran array A
! --> MX number of parallel Y-vectors inverted
!=============================================================================
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

!=============================================================================
subroutine LINBV(m,mah1,mah2,a,v)!                                     [LINBV]
!=============================================================================
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!		    SUBROUTINE LINBV
!   Solve LINear system with square Banded-matrix and vector V
!
! <-> A system matrix on input, its [L]*[D**-1]*[U] factorization on exit
! <-> V vector of right-hand-sides on input, solution vector on exit
! --> M order of matrix A
! --> MAH1 left half-bandwidth of A
! --> MAH2 right half-bandwidth of A
!=============================================================================
implicit none
integer(spi),intent(IN   ) :: m, mah1, mah2
real(sp),    intent(INOUT) :: a(m,-mah1:mah2), v(m)
!=============================================================================
call LDUB(m,mah1,mah2,a)
call UDLBV(m,mah1,mah2,a,v)
end subroutine LINBV

!=============================================================================
subroutine WRTB(m1,m2,mah1,mah2,a)!                                     [WRTB]
!=============================================================================
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
