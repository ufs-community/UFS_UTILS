module gengrid_kinds

  implicit none

  integer,parameter :: real_kind = selected_real_kind( 6) ! 4 byte real
  integer,parameter ::  dbl_kind = selected_real_kind(12) ! 8 byte real

  integer,parameter ::  int_kind = selected_int_kind ( 6) ! 4 byte integer
  integer,parameter :: int8_kind = selected_int_kind (13) ! 8 byte integer

  integer, parameter :: CL = 256
  integer, parameter :: CM =  64
  integer, parameter :: CS =  24

end module gengrid_kinds
