!> @file
!! @brief Holds machine dependent constants for global_cycle.
!! @author M. Iredell, xuli, Hang Lei, George Gayno

!> Holds machine dependent constants for global_cycle.
!! @author M. Iredell, xuli, Hang Lei, George Gayno
MODULE MACHINE

 IMPLICIT NONE
 SAVE
!  Machine dependant constants
 integer, parameter :: kind_io4  = 4
 integer, parameter :: kind_io8  = 8
 integer, parameter :: kind_ior = 8
 integer, parameter :: kind_evod = 8
 integer, parameter :: kind_dbl_prec = 8
 integer, parameter :: kind_rad  = selected_real_kind(13,60) !< the '60' maps to 64-bit real
 integer, parameter :: kind_phys = selected_real_kind(13,60) !< the '60' maps to 64-bit real
 integer, parameter :: kind_REAL = 8 !< used in cmp_comm
 integer, parameter :: kind_INTEGER = 4 !< -,,-
!
 real(kind=kind_evod), parameter :: mprec = 1.e-12 !< machine precision to restrict dep

 END MODULE MACHINE
