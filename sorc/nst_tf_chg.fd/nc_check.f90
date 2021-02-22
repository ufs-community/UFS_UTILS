!> @file
!! @brief Contains subroutine to Check netCDF return code.

!> Check netCDF return code, printing error message for errors.
!!
!! @param status the error code to check
!! @author Xu Li @date Mar, 2017
subroutine nc_check(status)

  use netcdf

  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine nc_check

