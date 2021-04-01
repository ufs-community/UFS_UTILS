!> @file
!! @brief Utilities for use when reading grib2 data.
!! @author George Gayno NCEP/EMC

!> Utilities for use when reading grib2 data.
!!
!! This module contains routines to:
!! - convert from RH to specific humidity
!! - convert from omega to dzdt.
!!
!! George Gayno NCEP/EMC
module grib2_util

use esmf

use model_grid, only      : i_input, j_input

implicit none

contains 

!> Convert relative humidity to specific humidity.
!!
!! @param[inout] rh_sphum rel humidity on input. spec hum on output.
!! @param[in] p pressure in Pa
!! @param[in] t temperature
!! @author Larissa Reames
!! @author Jeff Beck
 subroutine rh2spfh(rh_sphum,p,t)
    
  implicit none
  real,parameter      :: alpha=-9.477E-4 , & !K^-1,
                         Tnot=273.15, &  !K
                         Lnot=2.5008E6, & !JKg^-1
                         Rv=461.51, & !JKg^-1K^-1
                         esnot=611.21 !Pa
  
  real(esmf_kind_r4), intent(inout), dimension(i_input,j_input) ::rh_sphum
  real(esmf_kind_r8), intent(in)                  :: p, t(i_input,j_input)

  real, dimension(i_input,j_input)  :: es, e, rh

  print*,"- CONVERT RH TO SPFH AT LEVEL ", p

  rh = rh_sphum
  !print *, 'T = ', T, ' RH = ', RH, ' P = ', P
  es = esnot * exp( Lnot/Rv * ((t-Tnot)/(t*tnot) + alpha * LOG(t/Tnot) - alpha * (t-Tnot)/ t))
  !print *, 'es = ', es
  e = rh * es / 100.0
  !print *, 'e = ', e
  rh_sphum = 0.622 * e / p
  !print *, 'q = ', sphum
  
  !if (P .eq. 100000.0) THEN
  ! print *, 'T = ', T, ' RH = ', RH, ' P = ', P, ' es = ', es, ' e = ', e, ' q = ', sphum
  !end if

end subroutine RH2SPFH

!> Convert omega to vertical velocity.
!!
!! @param[inout] omega on input, vertical velocity on output
!! @param[in] p pressure
!! @param[in] t temperature
!! @param[in] q specific humidity
!! @param[in] clb lower bounds of indices processed by this mpi task
!! @param[in] cub upper bounds of indices processed by this mpi task
!! @author Larissa Reames
!! @author Jeff Beck
subroutine convert_omega(omega,p,t,q,clb,cub)

  implicit none
  real(esmf_kind_r8), pointer     :: omega(:,:,:), p(:,:,:), t(:,:,:), q(:,:,:),omtmp,ptmp
  
  integer                         :: clb(3), cub(3), i ,j, k
  
  real, parameter                 :: Rd = 287.15_esmf_kind_r8, &  !JKg^-1K^-1
                                     Rv=461.51_esmf_kind_r8, & !JKg^-1K^-1
                                     g = 9.81_esmf_kind_r8 ! ms^-2
                                     
  real(esmf_kind_r8)              :: tv, w
  
  do k = clb(3),cub(3)
    do j = clb(2),cub(2)
      do i = clb(1),cub(1)
        tv = t(i,j,k)*(1+Rd/Rv*q(i,j,k))
        omtmp=>omega(i,j,k)
        ptmp=>p(i,j,k)

        w = -1 * omtmp * Rd * tv / (ptmp * g)
        omega(i,j,k)=w
      enddo
    enddo
  enddo

end subroutine convert_omega

 end module grib2_util
