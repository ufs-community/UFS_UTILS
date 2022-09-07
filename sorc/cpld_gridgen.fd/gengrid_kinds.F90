!> @file
!! @brief Define type-kind variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module defines the type-kind variables
!! @author Denise.Worthen@noaa.gov

module gengrid_kinds

  implicit none

  integer,parameter :: real_kind = selected_real_kind( 6) !< 4 byte real
  integer,parameter ::  dbl_kind = selected_real_kind(12) !< 8 byte real

  integer,parameter ::  int_kind = selected_int_kind ( 6) !< 4 byte integer
  integer,parameter :: int8_kind = selected_int_kind (13) !< 8 byte integer

  integer, parameter :: CL = 256                          !< a long length character string
  integer, parameter :: CM =  64                          !< a medium length character string
  integer, parameter :: CS =  24                          !< a short length character string

end module gengrid_kinds
