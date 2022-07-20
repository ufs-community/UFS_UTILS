!> @file
!! @brief Holds machine dependent constants for global_cycle.
!! @author Mark Iredell NOAA/EMC

!> Holds machine dependent constants for global_cycle.
!! @author Mark Iredell NOAA/EMC
MODULE MACHINE

 IMPLICIT NONE
 SAVE
 integer, parameter :: kind_sngl_prec  = 4 !< Kind type for 4-byte float point
                                           !! variables.
 integer, parameter :: kind_io4  = kind_sngl_prec !< Kind type for 4-byte float point
                                                  !! variables.
 integer, parameter :: kind_dbl_prec  = 8 !< Kind type for 8-byte float point
                                          !! variables.
 integer, parameter :: kind_io8  = kind_dbl_prec !< Kind type for 8-byte float point
                                          !! variables.
 END MODULE MACHINE
