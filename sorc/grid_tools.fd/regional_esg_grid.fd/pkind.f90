!> @file
!! @brief Standard integer, real, and complex single and double precision kinds.
!! @author R. J. Purser

!> Standard integer, real, and complex single and double precision kinds.
!! @author R. J. Purser
module pkind
integer,parameter:: spi=selected_int_kind(6) !< Single precision integer kind.
integer,parameter:: dpi=selected_int_kind(12) !< Double precision integer kind.
integer,parameter:: sp =selected_real_kind(6,30) !< Single precision real kind.
integer,parameter:: dp =selected_real_kind(15,300) !< Double precision real kind.
integer,parameter:: spc=sp  !< Single precision real kind.
integer,parameter:: dpc=dp  !< Double precision real kind.
!private:: one_dpi; integer(8),parameter:: one_dpi=1
!integer,parameter:: dpi=kind(one_dpi)
!integer,parameter:: sp=kind(1.0)
!integer,parameter:: dp=kind(1.0d0)
!integer,parameter:: spc=kind((1.0,1.0))
!integer,parameter:: dpc=kind((1.0d0,1.0d0))
end module pkind
