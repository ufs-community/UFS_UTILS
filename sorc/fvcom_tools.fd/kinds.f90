!> @file
!! @brief Hold specification kinds for variable declaration.
!! @author David Wright, University of Michigan

!> Module to hold specification kinds for variable declaration.
!!
!! This module is based on (copied from) Paul vanDelst's type_kinds
!! module found in the community radiative transfer model
!!
!! @author David Wright, University of Michigan
module kinds
  implicit none
  private
!
! for name string
  integer, parameter, public  :: len_sta_name  = 8 !< Name length.

! Integer type definitions below

! Integer types
  integer, parameter, public  :: i_kind = 4 !< generic specification kind for default integer.
  integer, parameter, public  :: i_short = 2 !< generic specification kind for default short.
  integer, parameter, public  :: i_byte  = 1 !< generic specification kind for default byte.
! Real types
  integer, parameter, public  :: r_single = 4  !< specification kind for single precision (4-byte) real variable.
  integer, parameter, public  :: r_kind = 8 !< generic specification kind for default floating point

!
  real(r_single),parameter,public :: rmissing=-99999.0 !< Fill value for single real missing data.
  real(i_kind),parameter,public :: imissing=-99999 !< Fill value for integer missing data.
  real(r_kind),parameter,public :: drmissing=-99999.0 !< Fill value for double real missing data.

end module kinds
