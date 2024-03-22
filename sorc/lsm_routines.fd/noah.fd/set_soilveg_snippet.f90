!> @file
!> @brief Routine to set Noah LSM soil and veg params needed for sflx_snippet
!> @author Clara Draper

!> Below was extracted from namelist_soilveg.f and set_soilveg.f 
!! (couldn't get above to compile for doxygen)

!> Add Noah-MP LSM soil and veg params needed for global_cycle
!> Noah-MP related parameters were extracted from noahmp_table.f
!> isot (soil type) = 1: STATSGO must be selected if NoahMP is used
!> ivet (vegetation type) = 1: IBGP is used by UFS offline Land DA for Noah-MP
!> as of 07/13/2023
!> @author Yuan Xue

module set_soilveg_snippet_mod

 implicit none

 private

 public set_soilveg_noah
 public set_soilveg_noahmp

contains

!> This subroutine initializes soil and vegetation
!! parameters needed in global_cycle/land_increment.f90 
!! @param[in] isot Soil type
!! @param[in] ivet Vegetation type
!! @param[out] maxsmc Maximum soil moisture for each soil type
!! @param[out] bb B exponent for each soil type
!! @param[out] satpsi Saturated matric potential for each soil type
!! @param[out] iret Return integer
subroutine set_soilveg_noah(isot,ivet, maxsmc, bb, satpsi, iret) 
  implicit none

  integer, intent(in) :: isot,ivet
  real, dimension(30), intent(out)  :: maxsmc, bb, satpsi
  integer, intent(out) :: iret

! set vegetation-dependent params (May 2021, UFS uses ivet=1) 
! Draper, not needed for now, but might need SNUPX 
! for SWE-> SCF calculation later
! 
!      if(ivet.eq.1)then

  !defined_veg=20
! might want this later
 ! SNUPX  =(/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,
 !*             0.020, 0.060, 0.040, 0.020, 0.010, 0.020,
 !*             0.020, 0.020, 0.013, 0.013, 0.010, 0.020,
 !&             0.020, 0.020, 0.000, 0.000, 0.000, 0.000,
 !&             0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

!      endif

! set soil-dependent params (May 2021, UFS uses isot=1) 

  if (isot .eq. 1) then

! using stasgo table
  BB         =(/4.05,  4.26, 4.74, 5.33, 5.33,  5.25, &
             6.77,  8.72,  8.17, 10.73, 10.39,  11.55,&
             5.25,  4.26,  4.05, 4.26,  11.55,  4.05, & 
             4.05,  0.00,  0.00, 0.00,  0.00,  0.00,  & 
             0.00,  0.00,  0.00, 0.00,  0.00,  0.00/)
! Draper, these are provided for reference only, and 
! may be useful for later SMC updates
!      DRYSMC=(/0.010, 0.025, 0.010, 0.010, 0.010, 0.010,
!     &            0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
!     &            0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
!     &            0.010, 0.000, 0.000, 0.000, 0.000, 0.000,
!     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

  MAXSMC=(/0.395, 0.421, 0.434, 0.476, 0.476, 0.439,   & 
             0.404, 0.464, 0.465, 0.406, 0.468, 0.457, &
             0.464, 0.421, 0.200, 0.421, 0.457, 0.200, & 
             0.395, 0.000, 0.000, 0.000, 0.000, 0.000, & 
             0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

  SATPSI=(/0.035, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548,   & 
             0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677,&
             0.3548, 0.0363, 0.0350, 0.0363, 0.4677, 0.0350,&
             0.0350, 0.00, 0.00, 0.00, 0.00, 0.00,          &
             0.00, 0.00, 0.00, 0.00, 0.00, 0.00/)

 !defined_soil=19
  else 
        print *, 'set_soilveg_snippet not coded for soil type ', isot
        iret = -1
        return
  endif 
  
  iret = 0

end subroutine set_soilveg_noah

!> This subroutine initializes soil and vegetation
!! parameters needed in global_cycle/land_increment.f90 for noah-mp
!! @param[in] isot Soil type
!! @param[in] ivet Vegetation type
!! @param[out] maxsmc Maximum soil moisture for each soil type
!! @param[out] bb B exponent for each soil type
!! @param[out] satpsi Saturated matric potential for each soil type
!! @param[out] iret Return integer
subroutine set_soilveg_noahmp(isot,ivet, maxsmc, bb, satpsi,iret)

  implicit none

  integer, intent(in) :: isot,ivet !ivet is *not* used for now
  real, dimension(30), intent(out)  :: maxsmc, bb, satpsi
  integer, intent(out) :: iret

 if (isot .eq. 1) then

! set soil-dependent params (STATSGO is the only option for UFS, 07/13/2023)
  maxsmc= (/0.339, 0.421, 0.434, 0.476, 0.484,&
     &   0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468,                    &
     &   0.439, 1.000, 0.200, 0.421, 0.468, 0.200,                           &
     &   0.339, 0.339, 0.000, 0.000, 0.000, 0.000,                           &
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
  bb= (/2.79,  4.26, 4.74, 5.33, 3.86,  5.25,&
     &    6.77,  8.72,  8.17, 10.73,  10.39, 11.55,                          &
     &    5.25,  0.0,  2.79, 4.26,  11.55,  2.79,                            &
     &    2.79,  0.00,  0.00, 0.00,  0.00,  0.00,                            &
     &    0.00,  0.00,  0.00, 0.00,  0.00,  0.00/)
  satpsi= (/0.069, 0.036, 0.141, 0.759, 0.955, &
     &   0.355, 0.135, 0.617, 0.263, 0.098, 0.324, 0.468,                    &
     &   0.355, 0.00, 0.069, 0.036, 0.468, 0.069,                            &
     &   0.069, 0.00, 0.00, 0.00, 0.00, 0.00,                                &
     &   0.00, 0.00, 0.00, 0.00, 0.00, 0.00/)

 else
    print*, 'For Noah-MP, set_soilveg is not supported for soil type ', isot
    iret = -1
    return

 endif

 iret = 0
end subroutine set_soilveg_noahmp

end module set_soilveg_snippet_mod
