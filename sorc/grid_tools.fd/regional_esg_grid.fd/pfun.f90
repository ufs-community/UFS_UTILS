!
!                                         *********************************
!                                         *          module pfun          *
!                                         *          R. J. Purser         *
!                                         *          NCEP/EMC 2017        *
!                                         *********************************
! Direct dependencies:
! Modules: pkind, pietc_s, pietc
!
!=============================================================================
module pfun
!=============================================================================
use pkind, only: sp,dp
implicit none
private
public:: gd,gdi,hav,havh,ahav,ahavh,sech,sechs,atanh,sinoxm,sinox,&
         sinhoxm,sinhox

interface gd;      module procedure gd_s,    gd_d;     end interface
interface gdi;     module procedure gdi_s,   gdi_d;    end interface
interface hav;     module procedure hav_s,   hav_d;    end interface
interface havh;    module procedure havh_s,  havh_d;   end interface
interface ahav;    module procedure ahav_s,  ahav_d;   end interface
interface ahavh;   module procedure ahavh_s, ahavh_d;  end interface
interface atanh;   module procedure atanh_s, atanh_d;  end interface
interface sech;    module procedure sech_s,  sech_d;   end interface
interface sechs;   module procedure sechs_s, sechs_d;  end interface
interface sinoxm;  module procedure          sinoxm_d; end interface
interface sinox;   module procedure          sinox_d;  end interface
interface sinhoxm; module procedure          sinhoxm_d;end interface
interface sinhox;  module procedure          sinhox_d; end interface

contains

!=============================================================================
function gd_s(x) result(y)!                                               [gd]
!=============================================================================
! Gudermannian function
implicit none
real(sp),intent(in ):: x
real(sp)            :: y
y=atan(sinh(x))
end function gd_s
!=============================================================================
function gd_d(x) result(y)!                                               [gd]
!=============================================================================
implicit none
real(dp),intent(in ):: x
real(dp)            :: y
y=atan(sinh(x))
end function gd_d

!=============================================================================
function gdi_s(y) result(x)!                                             [gdi]
!=============================================================================
! Inverse Gudermannian function
implicit none
real(sp),intent(in ):: y
real(sp)            :: x
x=atanh(sin(y))
end function gdi_s
!=============================================================================
function gdi_d(y) result(x)!                                             [gdi]
!=============================================================================
implicit none
real(dp),intent(in ):: y
real(dp)            :: x
x=atanh(sin(y))
end function gdi_d

!=============================================================================
function hav_s(t) result(a)!                                             [hav]
!=============================================================================
! Haversine function
use pietc_s, only: o2
implicit none
real(sp),intent(in ):: t
real(sp)            :: a
a=(sin(t*o2))**2
end function hav_s
!=============================================================================
function hav_d(t) result(a)!                                             [hav]
!=============================================================================
use pietc, only: o2
implicit none
real(dp),intent(in ):: t
real(dp)            :: a
a=(sin(t*o2))**2
end function hav_d

!=============================================================================
function havh_s(t) result(a)!                                           [havh]
!=============================================================================
! Note the minus sign in the hyperbolic-haversine definition
use pietc_s, only: o2
implicit none
real(sp),intent(in ):: t
real(sp)            :: a
a=-(sinh(t*o2))**2
end function havh_s
!=============================================================================
function havh_d(t) result(a)!                                           [havh]
!=============================================================================
use pietc, only: o2
implicit none
real(dp),intent(in ):: t
real(dp)            :: a
a=-(sinh(t*o2))**2
end function havh_d

!=============================================================================
function ahav_s(a) result(t)!                                           [ahav]
!=============================================================================
use pietc_s, only: u2
! Arc-haversine function
implicit none
real(sp),intent(in ):: a
real(sp)            :: t
t=u2*asin(sqrt(a))
end function ahav_s
!=============================================================================
function ahav_d(a) result(t)!                                           [ahav]
!=============================================================================
use pietc, only: u2
implicit none
real(dp),intent(in ):: a
real(dp)            :: t
t=u2*asin(sqrt(a))
end function ahav_d

!=============================================================================
function ahavh_s(a) result(t)!                                         [ahavh]
!=============================================================================
use pietc_s, only: u2
! Note the minus sign in the hyperbolic arc-haversine definition
implicit none
real(sp),intent(in ):: a
real(sp)            :: t
t=u2*asinh(sqrt(-a))
end function ahavh_s
!=============================================================================
function ahavh_d(a) result(t)!                                         [ahavh]
!=============================================================================
use pietc, only: u2
implicit none
real(dp),intent(in ):: a
real(dp)            :: t
t=u2*asinh(sqrt(-a))
end function ahavh_d

!=============================================================================
function atanh_s(t) result(a)!                                         [atanh]
!=============================================================================
use pietc_s, only: u1,o2,o3,o5
implicit none
real(sp),intent(IN ):: t
real(sp)            :: a,tt
real(sp),parameter  :: o7=u1/7_sp,o9=u1/9_sp
!=============================================================================
if(abs(t)>=u1)stop 'In atanh; no solution'
if(abs(t)>1.e-3_sp)then; a=log((u1+t)/(u1-t))*o2
else; tt=t*t;            a=t*(u1+tt*(o3+tt*(o5+tt*(o7+tt*o9))))
endif
end function atanh_s
!=============================================================================
function atanh_d(t) result(a)!                                         [atanh]
!=============================================================================
use pietc, only: u1,o2,o3,o5
implicit none
real(dp),intent(IN ):: t
real(dp)            :: a,tt
real(dp),parameter  :: o7=u1/7_dp,o9=u1/9_dp
!=============================================================================
if(abs(t)>=u1)stop 'In atanh; no solution'
if(abs(t)>1.e-3_dp)then; a=log((u1+t)/(u1-t))*o2
else; tt=t*t;            a=t*(u1+tt*(o3+tt*(o5+tt*(o7+tt*o9))))
endif
end function atanh_d

!=============================================================================
function sech_s(x)result(r)!                                            [sech]
!=============================================================================
! This indirect way of computing 1/cosh(x) avoids overflows at large x
use pietc_s, only: u1,u2
implicit none
real(sp),intent(in ):: x
real(sp)            :: r
real(sp)            :: e,ax
ax=abs(x)
e=exp(-ax)
r=e*u2/(u1+e*e)
end function sech_s
!=============================================================================
function sech_d(x)result(r)!                                            [sech]
!=============================================================================
use pietc, only: u1,u2
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
real(dp)            :: e,ax
ax=abs(x)
e=exp(-ax)
r=e*u2/(u1+e*e)
end function sech_d

!=============================================================================
function sechs_s(x)result(r)!                                          [sechs]
!=============================================================================
implicit none
real(sp),intent(in ):: x
real(sp)            :: r
r=sech(x)**2
end function sechs_s
!=============================================================================
function sechs_d(x)result(r)!                                          [sechs]
!=============================================================================
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
r=sech(x)**2
end function sechs_d

!=============================================================================
function sinoxm_d(x) result(r)!                                       [sinoxm]
!=============================================================================
! Evaluate the symmetric real function sin(x)/x-1
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!-----------------------------------------------------------------------------
real(dp):: xx
!=============================================================================
xx=x*x
if(xx > .05_dp)then; r=sin(x)/x-u1
else             ; r=-xx*(u1-xx*(u1-xx*(u1-xx*(u1-xx*(u1-xx/&
                     156._dp)/110._dp)/72._dp)/42._dp)/20._dp)/6._dp
endif
end function sinoxm_d

!=============================================================================
function sinox_d(x) result(r)!                                         [sinox]
!=============================================================================
! Evaluate the symmetric real function sin(x)/x
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!=============================================================================
r=sinoxm(x)+u1
end function sinox_d

!=============================================================================
function sinhoxm_d(x) result(r)!                                     [sinhoxm]
!=============================================================================
! Evaluate the symmetric real function sinh(x)/x-1
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!-----------------------------------------------------------------------------
real(dp):: xx
!=============================================================================
xx=x*x
if(xx > .05_dp)then; r=sinh(x)/x-u1
else;              r=xx*(u1+xx*(u1+xx*(u1+xx*(u1+xx*(u1+xx/&
                     156._dp)/110._dp)/72._dp)/42._dp)/20._dp)/6._dp
endif
end function sinhoxm_d

!=============================================================================
function sinhox_d(x) result(r)!                                       [sinhox]
!=============================================================================
! Evaluate the symmetric real function sinh(x)/x
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!=============================================================================
r=sinhoxm(x)+u1
end function sinhox_d

end module pfun
