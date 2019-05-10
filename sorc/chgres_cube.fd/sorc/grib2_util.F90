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
! iso2sig               Main code for conversion of isobaric to sigma coordinates
!
! p2hyo, p2hyb          Routines for pressure to sigma conversion; used by iso2sig 
!--------------------------------------------------------------------------

use esmf
use netcdf
use nemsio_module

use program_setup, only   : tracers_input,num_tracers, external_model

use model_grid, only      : i_input,j_input, ip1_input, jp1_input

use atmdata_type

implicit none


contains 


subroutine read_vcoord(isnative,vcoordi,vcoordo,lev_input,levp1_input,metadata,iret)

  implicit none
  integer, intent(in)                       :: lev_input, levp1_input
  logical, intent(in)                       :: isnative
  character (len=500), intent(in)           :: metadata
  real(esmf_kind_r8), intent(in)                          :: vcoordi(lev_input)
  real(esmf_kind_r8), intent(inout), allocatable  :: vcoordo(:,:)
  integer, intent(out)                      :: iret
  
  integer :: k, idate(3)
  real :: b1, b2, b3, b4, b5, eta(51), etac
  
  eta = (/1.0000, 0.9980, 0.9940, 0.9870, 0.9750, 0.9590, &
              0.9390, 0.9160, 0.8920, 0.8650, 0.8350, 0.8020, 0.7660, &
              0.7270, 0.6850, 0.6400, 0.5920, 0.5420, 0.4970, 0.4565, &
              0.4205, 0.3877, 0.3582, 0.3317, 0.3078, 0.2863, 0.2670, &
              0.2496, 0.2329, 0.2188, 0.2047, 0.1906, 0.1765, 0.1624, &
              0.1483, 0.1342, 0.1201, 0.1060, 0.0919, 0.0778, 0.0657, &
              0.0568, 0.0486, 0.0409, 0.0337, 0.0271, 0.0209, 0.0151, &
              0.0097, 0.0047, 0.0000/)

  !desc: D=YYYYMMDDHHmmss:RH:xxx mb:etc
  read(metadata(3:6),'(I4)') idate(1)
  read(metadata(7:8),'(I2)') idate(2)
  read(metadata(9:10),'(I2)') idate(3)


  if (isnative) then 
    if ((trim(external_model) .eq. 'HRRR' .or. trim(external_model) .eq. 'RAP') & 
        .and. lev_input == 50) then 
      if (idate(1) .le. 2018 .and. idate(2) .le. 7 .and. idate(3) .lt. 12) then !old sigma coordinates
        vcoordo(:,2) = eta
      else  !new hybrid levels
        etac = 0.2
        b1 = 2. * etac**2 * ( 1. - etac )
        b2 = -etac * ( 4. - 3. * etac - etac**3 )
        b3 = 2. * ( 1. - etac**3 )
        b4 = - ( 1. - etac**2 )
        b5 = (1.-etac)**4
        do k=1,lev_input
          vcoordo(k,2) = ( b1 + b2*eta(k) + b3*eta(k)**2 + b4*eta(k)**3 ) / b5
          if ( eta(k) .lt. etac ) then
            vcoordo(k,2) = 0.
          end if
          if ( k .eq. 0 ) then
            vcoordo(k,2) = 1.
          else if ( k .eq. lev_input ) then
            vcoordo(k,2) = 0.
          end if
        enddo
      endif !end date check for sigma vs. hybrid
    elseif(trim(external_model) == "NAM" .and. lev_input == 60) then
      vcoordo(:,1) = (/0.0,1.1310602*10**-2,2.2672068*10**-2,3.4080435*10**-2,4.570657*10**-2,5.7557303*10**-2, &
      6.9626026*10**-2,8.1914276*10**-2,9.4423383*10**-2,0.1073834,0.1207936,0.1346522,0.1491913, & 
      0.1644059,0.1802934,0.1970789,0.215691,0.2363433,0.261099,0.288769,0.320698,0.3568189, &
      0.3972864,0.44291,0.493114,0.5477551,0.6066551,0.6687292,0.733094,0.7977924,0.8627368, &
      0.9268976,0.9893017,1.048898,1.104175,1.152405,1.189998,1.211094,1.20741,1.169137, &
      1.094331,0.0,4.0614128*10**-2.,8.1287362*10**-2,0.1221084,0.1630771,0.2042676,0.2460121, &
      0.2884955,0.3317178,0.3764177,0.4222259,0.4695117,0.5190141,0.5703636,0.6254073, &
      0.6852535,0.7521186,0.8263722,0.9083837,1.0/)
      vcoordo(:,2) = (/1.0,0.9921913,0.9843418,0.9764536,0.9684085,0.9602016,0.9518366, &
      0.9433123,0.9346272,0.925621,0.9162929,0.9066435,0.8965098,0.8858933,0.8747947, &
      0.8630539,0.8500177,0.8355303,0.8181331,0.798646,0.776103,0.7505253,0.7217709, &
      0.6892209,0.653231,0.61384,0.5710961,0.525693,0.4781699,0.4298564,0.3806763, &
      0.3312084,0.281916,0.2331837,0.1855005,0.1399449,9.756881*10**-2,6.0258172*10**-2, &
      3.0155413*10**-2,1.0215664*10**-2,1.4652965E-03,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
    else  
      write(*,*) "this code only supports rap/hrrr data w/ 50 sigma/hybrid coordinate levels &
      or NAM data with 60 levels"
      iret = 1
      return
    endif ! end check for mname and num levels
    
  else ! create sigma from isobaric levels
      vcoordo(2:levp1_input,2) = vcoordi / 100000.0
      vcoordo(1,2) = 0.0
  endif ! end native vs. non-native check

  vcoordo(:,1) = 0.0
  iret=0

end subroutine read_vcoord

subroutine iso2sig(pi,sigma,lev_input,levp1_input,psptr,atm,clb,cub,nvars,iret)

  implicit none
  
  real(esmf_kind_r8), intent(inout)                       :: pi(lev_input), & 
                                                           sigma(levp1_input,2)
  integer, intent(inout)                                  :: lev_input, levp1_input, &
                                                              clb(3),cub(3)
  type(atmdata), intent(inout)   :: atm(:)
  real(esmf_kind_r8), pointer, intent(inout)               :: psptr(:,:)
  integer, intent(out)                                  :: iret
  
  
  real                                      :: msg
  integer                                   :: nvars

 
  call p2hyo(pi,atm,psptr,100000.0,sigma(2:levp1_input,1),sigma(2:levp1_input,2), &
              lev_input,clb,cub,nvars,msg,4,iret)
  if (iret /= 0) call error_handler(" CONVERTING TO SIGMA COORDINATES. ONE OR BOTH PRESSURE &
                                    &ARRAYS ARE NOT IN TOP TO BOTTOM ORDER.", iret)
  iret = 0
   
end subroutine iso2sig

 subroutine p2hyo(pi,xi,psfc,p0,hyao,hybo,lev_input,clb,cub,nvars,xmsg,kflag,iret) 
      implicit none
      
!      this routine interploates constant pres levels to hybrid
!     the formula for the pressure of a hybrid surface is;
!          phy(k) = pout(k) = hya(k)*p0 + hyb(k)*psfc
!
!     input  ["i" input ... "o" output]
!          pi     - pressure level                     [input]
!          psfc   - is the surface pressure pa         [input]
!          mlon   - longitude dimension
!          nlat   - latitude  dimension
!          hyao   - is the "a" or pressure hybrid coef
!          hybo   - is the "b" or sigma coeficient
!          klevo  - number of output levels
!          kflag  - specify how values outside the "pi" will be handled
!                   by "outside" i mean [pout < pi(1)] or [pout > pi(klevi)]
!                   extrapolation is via log-linear extrapolation.
!                   =0   no extrapolation. values set to _fillvalue
!                   =1   values set to nearest valid value
!                   =2   values at pout less    than pi(1) set to nearest value
!                        values at pout greater than pi(klevi) are extrapolated 
!                   =3   values at pout less    than pi(1) are extrapolated
!                        values at pout greater than pi(klevi) set to nearest value 
!                   =4   values at pout less    than pi(1) are extrapolated
!                        values at pout greater than pi(klevi) are extrapolated
!               
!          iret    - error code  [=0 no error detected]
!                               [.ne.0 error detected: one or both
!                                      pressure arrays not top->bot order]
!     output
!          xo     - pressure at hybrid levels [pa]
      
      integer, intent(inout)       :: lev_input
      integer                      :: kflag,clb(3),cub(3)
      real, intent(in)             :: p0,xmsg
      real, intent(inout)          :: pi(lev_input), hyao(lev_input), hybo(lev_input)
      type(atmdata), intent(inout):: xi(:)
      real(esmf_kind_r8), pointer, intent(inout)       :: psfc(:,:)
      
      integer, intent(out)      :: iret
      integer                   :: iflag, klevo, nvars
      
      real          :: po(lev_input)
      
      iflag = 0
!                                                 ! ? input asending order     
      iret   = 0
      klevo = lev_input
      
      !print *, hyao(klevo), ' ', hybo(klevo), ' ', psfc(1,1), ' ', p0
      ! print *, pi
      
      if (pi(1).gt.pi(lev_input)) then
          iret = 1
          return
      end if
    

      po(1)     = hyao(1)*p0     + hybo(1)*psfc(clb(1),clb(2))
      
      po(klevo) = hyao(klevo)*p0 + hybo(klevo)*psfc(clb(1),clb(2))
!                                                 ! ? output ascending order     
      if (po(1).gt.po(klevo)) then
          iret = 20 + iret
      end if

      if (iret /= 0) return

      call p2hyb(pi,xi,psfc,p0,hyao,hybo,po, lev_input,clb,cub,nvars,iflag, kflag, xmsg)
      
      return
end subroutine p2hyo

subroutine p2hyb(pi,xi,psfc,p0,hyao,hybo,po,lev_input,clb,cub,nvar,iflag, kflag, xmsg)
  implicit none
  
  integer, intent(in)         :: lev_input,clb(3),cub(3)
  integer, intent(inout)      :: kflag
  real, intent(in)          :: p0,pi(lev_input), hyao(lev_input), hybo(lev_input), xmsg
  real, intent(inout)       :: po(lev_input)
  type(atmdata), intent(inout):: xi(:)
  real(esmf_kind_r8), pointer, intent(inout)       :: psfc(:,:)
  
  integer, intent(out)      :: iflag
  
  integer     :: nl,nlo,ml,mlo,ki,ko,nv,  dims(3),nvar
  real        :: pimin, pimax, pomin, pomax, dxdp
  real, allocatable :: xo(:,:,:,:)


  dims = shape(xi(1)%var)
  allocate(xo(dims(1),dims(2),dims(3),nvar+1))

  pimin = pi(1)   
  pimax = pi(lev_input)

  do nl = clb(2),cub(2)
    nlo = nl-clb(2)+1
    
    do ml = clb(1),cub(1)
      mlo = ml-clb(1)+1
      
      do ko = clb(3),cub(3)
        po(ko) = hyao(ko)*p0 + hybo(ko)*psfc(ml,nl)
        xo(mlo,nlo,ko,1) = po(ko)     
      end do

      pomin = po(1)
      pomax = po(lev_input)
    

      do ko = 1,lev_input
        xo(mlo,nlo,ko,2:nvar) = xmsg
    
        do ki = 1,lev_input-1
          if (po(ko) >= pimin .and. po(ko) <= pimax ) then     
            if (po(ko).ge.pi(ki) .and. po(ko).lt.pi(ki+1)) then
              do nv = 2,nvar
                xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,ki)            &
                               +(xi(nv)%var(ml,nl,ki+1) -xi(nv)%var(ml,nl,ki) )*    &
                                (log(po(ko))  -log(pi(ki)))/    &
                                (log(pi(ki+1))-log(pi(ki)))
              enddo 
            end if
           else
            if (kflag == 0) then
              iflag = 1
            elseif (kflag == 1) then
              if (po(ko) < pimin) then
                do nv = 2,nvar
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,1)
                enddo
              elseif (po(ko) > pimax) then
                do nv = 2,nvar
                  xo(mlo,nlo,ko,nv-1) = xi(nv)%var(ml,nl,lev_input)
                enddo
              end if
            elseif (kflag == 2) then
              if (po(ko) < pimin) then
                do nv = 2,nvar
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,1)
                enddo
              elseif (po(ko) > pimax) then
                do nv = 2,nvar
                  dxdp = (xi(nv)%var(ml,nl,lev_input) -xi(nv)%var(ml,nl,lev_input-1))*   &
                     (log(pi(lev_input))-log(pi(lev_input-1)))
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,lev_input)       &
                     + (log(po(ko))-log(pi(lev_input)))*dxdp
                enddo

              end if
            elseif (kflag == 3) then
              if (po(ko) < pimin) then
                do nv = 2,nvar
                  dxdp = (xi(nv)%var(ml,nl,2) -xi(nv)%var(ml,nl,1))*      &
                       (log(pi(2))-log(pi(1)))
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,1)          &
                       + (log(po(ko))-log(pi(1)))*dxdp
                enddo
  
              elseif (po(ko) > pimax) then
                do nv = 2,nvar       
                  xo(mlo,nlo,ko,nv-1) = xi(nv)%var(ml,nl,lev_input)
                enddo   
              end if
            elseif (kflag == 4) then
              if (po(ko) < pimin) then
                do nv = 2,nvar

                  dxdp = (xi(nv)%var(ml,nl,2) -xi(nv)%var(ml,nl,1))*       &
                         (log(pi(2))-log(pi(1)))
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,1)         &
                         + (log(po(ko))-log(pi(1)))*dxdp
                enddo
              elseif (po(ko) > pimax) then
                do nv = 2,nvar
                  dxdp = (xi(nv)%var(ml,nl,lev_input) -xi(nv)%var(ml,nl,lev_input-1))*   &
                         (log(pi(lev_input))-log(pi(lev_input-1)))
                  xo(mlo,nlo,ko,nv) = xi(nv)%var(ml,nl,lev_input)       &
                         + (log(po(ko))-log(pi(lev_input)))*dxdp   
                enddo
              end if
            end if !kflag 
          endif ! outside pi bounds
         end do ! loop over klevs input
      end do  !loop over klevs output
    end do !loop over lon
  end do  !loop over lat
  
  do nv = 1,nvar
    xi(nv)%var(clb(1):cub(1),clb(2):cub(2),clb(3):cub(3)) = xo(:,:,lev_input:1:-1,nv)
  enddo
  

  deallocate(xo)
  
 end subroutine p2hyb
 
 subroutine rh2spfh(rh_sphum,p,tptr,lev_cur)
    
  implicit none
  real,parameter      :: alpha=-9.477E-4 , & !K^-1,
                         Tnot=273.15, &  !K
                         Lnot=2.5008E6, & !JKg^-1
                         Rv=461.51, & !JKg^-1K^-1
                         esnot=611.21 !Pa
  
  real(esmf_kind_r4), intent(inout)    :: rh_sphum(i_input,j_input)
  real(esmf_kind_r8), intent(in)                  :: p
  real, pointer, intent(inout)      :: tptr(:,:,:)
  integer                           :: lev_cur
  real, dimension(i_input,j_input)  :: es, e
  real, pointer                     :: t(:,:)  

    
  t => tptr(:,:,lev_cur)
  !print *, 'T = ', T, ' RH = ', RH, ' P = ', P
  es = esnot * exp( Lnot/Rv * ((T-Tnot)/(T*tnot) + alpha * LOG(T/Tnot) - alpha * (T-Tnot)/ T))
  !print *, 'es = ', es
  e = rh_sphum * es / 100.0
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