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

public :: rh2spfh
public :: convert_omega
public :: rh2spfh_gfs
public :: fpvsnew 

contains 

!> Convert relative humidity to specific humidity.
!> Calculation of saturation water vapor pressure is based on
!> Brock and Richardson 2001 (Meterological Measurement 
!> Systems, p. 86, equation 5.1) 
!!
!! @param[inout] rh_sphum rel humidity (%) on input. spec hum (kg/kg) on output.
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
  rh_sphum = real((0.622 * e / p),kind=esmf_kind_r4)
  !print *, 'q = ', sphum
  
  !if (P .eq. 100000.0) THEN
  !print *, 'T = ', t, ' RH = ', rh, ' P = ', p, ' es = ', es, ' e = ', e, ' q = ', rh_sphum
  !end if

end subroutine RH2SPFH

!> Convert relative humidity to specific humidity (GFS formula)
!> Calculation of saturation water vapor pressure is based on
!> GFS function fvpsnew (Phillips 1982). The model does account for the variation of the
!> latent heat of condensation with temperature. A linear interpolation is done
!> between values in a calculated lookup table. Ice and water are considered
!> separately. This option provides a consistent conversion for GFS grib2 data
!> to be ingested.
!!
!! @param[inout] rh_sphum rel humidity on input. spec hum on output.
!! @param[in] p pressure in Pa
!! @param[in] t temperature
!! @author Jili Dong NCEP/EMC 
 subroutine rh2spfh_gfs(rh_sphum,p,t)

  implicit none

 integer kind_phys

 parameter (kind_phys = selected_real_kind(13,60)) ! the '60' maps to 64-bit real


 real(kind=kind_phys),parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)                         
 real(kind=kind_phys),parameter:: con_rv      =4.6150e+2 ! gas constant H2O    (J/kg/K)

 real(kind=kind_phys),parameter:: con_eps     =con_rd/con_rv
 real(kind=kind_phys),parameter:: con_epsm1   =con_rd/con_rv-1.




  real(esmf_kind_r4), intent(inout), dimension(i_input,j_input) ::rh_sphum
  real(esmf_kind_r8), intent(in)                  :: p, t(i_input,j_input)

  real, dimension(i_input,j_input)  :: QC, rh
  real :: ES
  integer :: i,j

  print*,"- CONVERT RH TO SPFH AT LEVEL ", p

  rh = rh_sphum

do j=1,j_input
  do i=1,i_input
    ES = MIN(FPVSNEW(T(I,J)),P)
    QC(i,j) = CON_EPS*ES/(P+CON_EPSM1*ES)
    rh_sphum(i,j) = real((rh(i,j)*QC(i,j)/100.0),kind=esmf_kind_r4)
  end do
end do


  !print *, 'T = ', T, ' RH = ', RH, ' P = ', P
  !print *, 'q = ', sphum


end subroutine RH2SPFH_GFS


!> Compute saturation vapor pressure 
!!
!! @param[in] t temperature in Kelvin
!! @return fpvsnew Saturation vapor pressure
!! @author N Phillips w/NMC2X2  
!! @date 30 dec 82

!
!-------------------------------------------------------------------------------------
!
      elemental function fpvsnew(t)
!
      implicit none
      integer,parameter:: nxpvs=7501
      real,parameter:: con_ttp     =2.7316e+2 ! temp at H2O 3pt
      real,parameter:: con_psat    =6.1078e+2 ! pres at H2O 3pt
      real,parameter:: con_cvap    =1.8460e+3 ! spec heat H2O gas   (J/kg/K)
      real,parameter:: con_cliq    =4.1855e+3 ! spec heat H2O liq
      real,parameter:: con_hvap    =2.5000e+6 ! lat heat H2O cond
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter:: con_csol    =2.1060e+3 ! spec heat H2O ice
      real,parameter:: con_hfus    =3.3358e+5 ! lat heat H2O fusion
      real,parameter:: tliq=con_ttp
      real,parameter:: tice=con_ttp-20.0
      real,parameter:: dldtl=con_cvap-con_cliq
      real,parameter:: heatl=con_hvap
      real,parameter:: xponal=-dldtl/con_rv
      real,parameter:: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
      real,parameter:: dldti=con_cvap-con_csol
      real,parameter:: heati=con_hvap+con_hfus
      real,parameter:: xponai=-dldti/con_rv
      real,parameter:: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
      real tr,w,pvl,pvi
      real fpvsnew
      real(esmf_kind_r8),intent(in):: t
      integer jx
      real  xj,x,tbpvs(nxpvs),xp1
      real xmin,xmax,xinc,c2xpvs,c1xpvs
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xmin=180.0
      xmax=330.0
      xinc=(xmax-xmin)/(nxpvs-1)
!   c1xpvs=1.-xmin/xinc
      c2xpvs=1./xinc
      c1xpvs=1.-xmin*c2xpvs
!    xj=min(max(c1xpvs+c2xpvs*t,1.0),real(nxpvs,krealfp))
      xj=min(max(c1xpvs+c2xpvs*t,1.0),float(nxpvs))
      jx=nint(min(xj,float(nxpvs)-1.0))
      x=xmin+(jx-1)*xinc

      tr=con_ttp/x
      if(x>=tliq) then
        tbpvs(jx)=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
      elseif(x<tice) then
        tbpvs(jx)=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
      else
        w=(t-tice)/(tliq-tice)
        pvl=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
        pvi=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
        tbpvs(jx)=w*pvl+(1.-w)*pvi
      endif

      xp1=xmin+(jx-1+1)*xinc

      tr=con_ttp/xp1
      if(xp1>=tliq) then
        tbpvs(jx+1)=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
      elseif(xp1<tice) then
        tbpvs(jx+1)=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
      else
        w=(t-tice)/(tliq-tice)
        pvl=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
        pvi=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
        tbpvs(jx+1)=w*pvl+(1.-w)*pvi
      endif

      fpvsnew=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end function fpvsnew



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
