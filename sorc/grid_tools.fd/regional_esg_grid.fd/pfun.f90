!> @file
!! @brief Routine evaluating useful functions not always available in Fortran.
!! @author R. J. Purser

!> This module is for evaluating several useful real-valued functions
!! that are not always available in Fortran compilers. These have applications
!! in map transformations, spherical and pseudo-spherical surface geometry,
!! probability distributions, and splines.
!!
!! @author R. J. Purser
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

!> Gudermannian function (related to Mercator map projections)
!!
!! @param[in] x single precision real argument of function
!! @return y Gudermannian function of x
!! @author R. J. Purser  
function gd_s(x) result(y)!                                               [gd]
! Gudermannian function
implicit none
real(sp),intent(in ):: x
real(sp)            :: y
y=atan(sinh(x))
end function gd_s

!> Gudermannian function (related to Mercator map projections)
!!
!! @param[in] x double precision real argument of function
!! @return y Gudermannian function of x
!! @author R. J. Purser  
function gd_d(x) result(y)!                                               [gd]
implicit none
real(dp),intent(in ):: x
real(dp)            :: y
y=atan(sinh(x))
end function gd_d

!> Inverse Gudermannian function for single precision real.
!!
!! @param[in] y single precision real argument
!! @return x inverse Gudermannian function of y
!! @author R. J. Purser  
function gdi_s(y) result(x)!                                             [gdi]
implicit none
real(sp),intent(in ):: y
real(sp)            :: x
x=atanh(sin(y))
end function gdi_s

!> Inverse Gudermannian function for double precision real.
!!
!! @param[in] y double precision real argument
!! @return x inverse Gudermannian function of y
!! @author R. J. Purser  
function gdi_d(y) result(x)!                                             [gdi]
implicit none
real(dp),intent(in ):: y
real(dp)            :: x
x=atanh(sin(y))
end function gdi_d

!> Haversine function for single precision real (for geometry on the sphere).
!!
!! @param[in] t single precision real argument
!! @return a haversine function of t
!! @author R. J. Purser  
function hav_s(t) result(a)!                                             [hav]
use pietc_s, only: o2
implicit none
real(sp),intent(in ):: t
real(sp)            :: a
a=(sin(t*o2))**2
end function hav_s

!>  Haversine function for double precision real (for geometry on the sphere).
!!
!! @param[in] t double precision real argument
!! @return a haversine function of t
!! @author R. J. Purser  
function hav_d(t) result(a)!                                             [hav]
use pietc, only: o2
implicit none
real(dp),intent(in ):: t
real(dp)            :: a
a=(sin(t*o2))**2
end function hav_d

!> Hyperbolic-haversine for single precision real (for pseudosphere geometry).
!!
!! @note The minus sign in the hyperbolic-haversine definition.
!!
!! @param[in] t single precision real argument
!! @return a hyperbolic-haversine function of t
!! @author R. J. Purser  
function havh_s(t) result(a)!                                           [havh]
use pietc_s, only: o2
implicit none
real(sp),intent(in ):: t
real(sp)            :: a
a=-(sinh(t*o2))**2
end function havh_s

!> Hyperbolic-haversine for double precision real (for pseudosphere geometry).
!!
!! @param[in] t double precision real argument
!! @return a hyperbolic-haversine function of t
!! @author R. J. Purser  
function havh_d(t) result(a)!                                           [havh]
use pietc, only: o2
implicit none
real(dp),intent(in ):: t
real(dp)            :: a
a=-(sinh(t*o2))**2
end function havh_d

!> Arc-haversine function for single precision real.
!!
!! @param[in] a single precision real argument
!! @return t arc-haversine function of a
!! @author R. J. Purser  
function ahav_s(a) result(t)!                                           [ahav]
use pietc_s, only: u2
implicit none
real(sp),intent(in ):: a
real(sp)            :: t
t=u2*asin(sqrt(a))
end function ahav_s

!> Arc-haversine function for double precision real.
!!
!! @param[in] a double precision real argument
!! @return t arc-haversine function of a
!! @author R. J. Purser  
function ahav_d(a) result(t)!                                           [ahav]
use pietc, only: u2
implicit none
real(dp),intent(in ):: a
real(dp)            :: t
t=u2*asin(sqrt(a))
end function ahav_d

!> Hyperbolic arc-haversine for single precision real.
!!
!! @note The minus sign in the hyperbolic arc-haversine definition.
!!
!! @param[in] a single precision real argument
!! @return t hyperbolic arc-haversine of a
!! @author R. J. Purser  
function ahavh_s(a) result(t)!                                         [ahavh]
use pietc_s, only: u2
implicit none
real(sp),intent(in ):: a
real(sp)            :: t
t=u2*asinh(sqrt(-a))
end function ahavh_s

!> Hyperbolic arc-haversine for double precision real.
!!
!! @param[in] a double precision real argument
!! @return t hyperbolic arc-haversine of a
!! @author R. J. Purser  
function ahavh_d(a) result(t)!                                         [ahavh]
use pietc, only: u2
implicit none
real(dp),intent(in ):: a
real(dp)            :: t
t=u2*asinh(sqrt(-a))
end function ahavh_d

!> Hyperbolic arc-tangent for single precision real. (compilers now have this)
!!
!! @param[in] t single precision real argument
!! @return a hyperbolic arc-tangent of t
!! @author R. J. Purser  
function atanh_s(t) result(a)!                                         [atanh]
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

!> Hyperbolic arc-tangent for double precision real.
!!
!! @param[in] t double precision real argument
!! @return a hyperbolic arc-tangent of t
!! @author R. J. Purser  
function atanh_d(t) result(a)!                                         [atanh]
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

!> Hyperbolic secant for single precision real.
!!
!! @param[in] x single precision real argument
!! @return r hyperbolic secant of x
!! @author R. J. Purser  
function sech_s(x)result(r)!                                            [sech]
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

!> Hyperbolic secant for double precision real.
!!
!! @param[in] x double precision real argument
!! @return r hyperbolic secant of x
!! @author R. J. Purser  
function sech_d(x)result(r)!                                            [sech]
use pietc, only: u1,u2
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
real(dp)            :: e,ax
ax=abs(x)
e=exp(-ax)
r=e*u2/(u1+e*e)
end function sech_d

!> Hyperbolic secant-squared function (logistic distribution).
!!
!! @param[in] x single precision real argument
!! @return r sech-squared of x
!! @author R. J. Purser  
function sechs_s(x)result(r)!                                          [sechs]
implicit none
real(sp),intent(in ):: x
real(sp)            :: r
r=sech(x)**2
end function sechs_s

!> Hyperbolic secant-squared function (logistic distribution).
!!
!! @param[in] x double precision real argument
!! @return r sech-squared of x
!! @author R. J. Purser  
function sechs_d(x)result(r)!                                          [sechs]
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
r=sech(x)**2
end function sechs_d

!> Evaluate the symmetric real function sin(x)/x-1, still accurate near x=0.
!!
!! @param[in] x double precision real argument
!! @return r sin(x)/x-1
!! @author R. J. Purser  
function sinoxm_d(x) result(r)!                                       [sinoxm]
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

!> Evaluate the symmetric real function sin(x)/x.
!!
!! @param[in] x double precision real argument
!! @return r sin(x)/x
!! @author R. J. Purser  
function sinox_d(x) result(r)!                                         [sinox]
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!=============================================================================
r=sinoxm(x)+u1
end function sinox_d

!> Evaluate the symmetric real function sinh(x)/x-1. still accurate near x=0.
!!
!! @param[in] x double precision real argument
!! @return r sinh(x)-1
!! @author R. J. Purser  
function sinhoxm_d(x) result(r)!                                     [sinhoxm]
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

!> Evaluate the symmetric real function sinh(x)/x.
!!
!! @param[in] x double precision real argument
!! @return r sinh(x)/x
!! @author R. J. Purser  
function sinhox_d(x) result(r)!                                       [sinhox]
use pietc, only: u1
implicit none
real(dp),intent(in ):: x
real(dp)            :: r
!=============================================================================
r=sinhoxm(x)+u1
end function sinhox_d

end module pfun
