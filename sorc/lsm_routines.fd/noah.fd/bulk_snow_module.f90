!> @file
!> @brief Routines to make DA updates to a bulk (single layer) snow model
!! such as that in Noah.
!> @author Clara Draper

module bulk_snow_module

 implicit none

 private

 public calc_density

contains

!> This subroutine calculates snow density from forecast fields.
!! density = SWE/SND where snow present.
!!         = average from snow forecasts over land, where snow not present
!! @param[in] lensfc Number of sfc grid cells
!! @param[in] rank Processor rank
!! @param[in] landmask Mask for land increments
!! @param[in] swe Snow Water Equivalent
!! @param[in] snd Snow Depth
!! @param[out] density Snow density [-]

 subroutine calc_density(lensfc, landmask, swe, snd, rank, density)

       implicit none

       integer, intent(in) :: lensfc, rank
       integer, intent(in) :: landmask(lensfc)
       real, intent(in)    :: swe(lensfc), snd(lensfc)
       real, intent(out)   :: density(lensfc)

       real :: dens_mean
       integer :: n

        ! density = swe/snd
        do n =1,lensfc
                if (snd(n) > 0.001 ) then
                        density(n) = swe(n)/snd(n)
                else 
                        density(n)=0.1 
                endif 
        enddo

        where (density < 0.0001) density = 0.1

        ! calculate mean density of snow over land
        if (count (landmask==2) > 0) then
                ! mean density over snow-covered land
                dens_mean = sum(density, mask = (landmask==2 )) &
                         / count (landmask==2)
                print *, 'mean density on rank ', rank,': ', dens_mean
        else
                dens_mean = 0.1  ! default value if have no snow in tile
                print *, 'no snow on rank ', rank, ' using default density ', dens_mean
        endif

        ! for grid cells with no valid density, fill in the average snodens
        where( swe <= 0.001 ) density = dens_mean

 end subroutine calc_density

end module bulk_snow_module
