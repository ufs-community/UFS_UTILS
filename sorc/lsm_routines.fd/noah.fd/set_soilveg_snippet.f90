!> @file
!> @brief Routine to set Noah LSM soil and veg params needed for sflx_snippet
!> @author Clara Draper

!> Below was extracted from namelist_soilveg.f and set_soilveg.f 
!! (couldn't get above to compile for doxygen)

module set_soilveg_snippet_mod

 implicit none

 private

 public set_soilveg

contains

!> This subroutine initializes soil and vegetation
!! parameters needed in global_cycle/land_increment.f90 
!! @param[in] isot Soil type
!! @param[in] ivet Vegetation type
!! @param[out] maxsmc Maximum soil moisture for each soil type
!! @param[out] bb B exponent for each soil type
!! @param[out] satpsi Saturated matric potential for each soil type
!! @param[out] iret Return integer
subroutine set_soilveg(isot,ivet, maxsmc, bb, satpsi, iret) 
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

end subroutine set_soilveg

end module set_soilveg_snippet_mod
