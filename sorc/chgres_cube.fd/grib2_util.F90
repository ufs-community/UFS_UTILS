module grib2_util

!--------------------------------------------------------------------------
! Module: grib2_util
!
! Abstract: Utilities for use when reading grib2 data.
!
! Main Subroutines:
! -------------------
! read_vcoord           Reads vertical coordinate data
!
!--------------------------------------------------------------------------

use esmf

use program_setup, only   : external_model, base_install_dir

use model_grid, only      : i_input, j_input

implicit none

contains 

subroutine read_vcoord(isnative,vcoordi,vcoordo,lev_input,levp1_input,pt,metadata,iret)

  implicit none
  integer, intent(in)                       :: lev_input, levp1_input
  logical, intent(in)                       :: isnative
  character (len=500), intent(in)           :: metadata
  real(esmf_kind_r8), intent(in)            :: vcoordi(lev_input)
  real(esmf_kind_r8), intent(inout), allocatable  :: vcoordo(:,:)
  real(esmf_kind_r8), intent(out)						:: pt
  integer, intent(out)                      :: iret
  
  integer :: k, idate(3)
  character (len=1000)											:: fname_coord
  character (len=20)												:: lev_type
  real(esmf_kind_r8)                        :: sigma, sigflat

  !desc: D=YYYYMMDDHHmmss:RH:xxx mb:etc
  read(metadata(3:6),'(I4)') idate(1)
  read(metadata(7:8),'(I2)') idate(2)
  read(metadata(9:10),'(I2)') idate(3)

  vcoordo(:,:) = 0.0
  sigflat = 0.1
  if (isnative) then 
    if ((trim(external_model) .eq. 'HRRR' .or. trim(external_model) .eq. 'RAP') & 
        .and. lev_input == 50) then 
      if (idate(1) .le. 2018 .and. idate(2) .le. 7 .and. idate(3) .lt. 12) then !old sigma coordinates
        lev_type = "sigma"
      else  !new hybrid levels
        lev_type = "hybrid"
      endif
    elseif (trim(external_model) .eq. 'NAM' .and. lev_input == 60) then
      lev_type = "hybrid"
    else  
      iret = 1
      call error_handler("This code only supports rap/hrrr data w/ 50 sigma/hybrid coordinate levels &
      or NAM data with 60 hybrid coordinate levels", iret)
    endif ! end check for mname and num levels
    
    fname_coord = trim(base_install_dir)//"/fix/fix_chgres/vertical_coordinate_"// &
    						trim(external_model)//"-"//trim(lev_type)//".txt"
    print*, fname_coord
		open(14, file=trim(fname_coord), form='formatted', iostat=iret)
		
		if (iret /= 0) call error_handler("OPENING VERTICAL COORDINATE FILE", iret)


		read(14, *, iostat=iret) pt
		if (iret /= 0) call error_handler("READING VERTICAL COORDINATE FILE", iret)

		do k = 1, levp1_input
			read(14, *, iostat=iret) vcoordo(k,1), vcoordo(k,2)
		enddo    
  else ! create sigma from isobaric levels
		!vcoordo(2:levp1_input,2) = vcoordi / 100000.0
		
		!vcoordo(1,2) = 0.0
		
		! convert isobaric to hybrid levels
		do k=2,levp1_input

			sigma = vcoordi(k-1)/100000.0

			if (sigma .le. sigflat) then
				vcoordo(k,1)=sigma
				vcoordo(k,2)=0
			else
				vcoordo(k,1)=sigflat*(1-sigma)/(1-sigflat)
				vcoordo(k,2)=(sigma-sigflat)/(1-sigflat)
			end if

		end do
  endif ! end native vs. non-native check
  iret=0

end subroutine read_vcoord

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
