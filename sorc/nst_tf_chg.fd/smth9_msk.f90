!> @file
!! @brief Contains 9-point smoother subroutine.

!> 9-point smoother subroutine.
!!
!! smth9_msk --- 9-point smoother with surface mask
!! by spatial average from a to b with ancillary surface mask (e.g.,
!! analysis grid => surface grid) for global arrays
!! to gurantee the averaged value (b) is determined by the
!! candidates (a) with the identical surface type from (a)
!!
!! @param[in] a real: 2-d array such as analysis increment at analysis grids
!! @param[out] b real: 2-d array 
!! @param[in] isli integer: 2-d array: surface mask (0 = water, 1 = land, 2 = sea ice) for a grids
!! @param[in] nlon integer: number of longitude of a & b
!! @param[in] nlat integer: number of latitude of a & b
!! @param[in] istyp integer: target surface type value
!!
!! @author  Xu Li @date 2017-03-13
 subroutine smth9_msk(a,b,isli,nlon,nlat,istyp)
! USES:

 implicit none

! INPUT:
 real   , dimension(nlon,nlat), intent(in   ) :: a
 integer, dimension(nlon,nlat), intent(in   ) :: isli

 integer, intent(in   ) :: nlat,nlon,istyp

!OUTPUT:
 real   , dimension(nlon,nlat), intent(inout) :: b

!Declare local variables
 integer :: i,j,ii,jj,ix,iy
 integer :: num

 real    :: bout,dlat,dlon

 b=a
!Loop over all grids of array b to get interpolated value
 do j = 1, nlat
   do i = 1, nlon
     if ( isli(i,j) == istyp ) then
!
!      get the averraged value with the nearby grids (in a) which has
!      the identical surface mask (in b) only

       num = 0
       bout  = 0.0
       do ii = i-1, i+1
          ix = ii
         if ( ix == 0 ) ix = nlon
         if ( ix ==  nlon + 1 ) ix = 1
         do jj = j-1, j+1
           iy = jj
           iy = max(min(nlat,iy),1)
           if ( isli(ix,iy) == istyp ) then
             bout  = bout  + a(ix,iy)
             num = num + 1
           endif
         enddo
       enddo

       if ( num > 0 ) then
          bout = bout/real(num)
       else
          write(*,'(a,I4,I3,I5)') 'smth9_msk: warning, no data to averageistyp, keep it as is (orginal value) num ; ',istyp,num
       endif              

       b(i,j)=bout

     endif                ! if ( isli(i,j) == istyp ) then
   enddo                    ! do i = 1, nlon
 enddo                      ! do j = 1, nlat



 end subroutine smth9_msk

