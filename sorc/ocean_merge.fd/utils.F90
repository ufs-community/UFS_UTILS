!> Check NetCDF return code. If an error is indicated,
!! stop program.
!!
!! @param[in] ret NetCDF return code.
!! @author Shan Sun
subroutine handle_err (ret)
  use netcdf
  implicit none
  integer, intent(in) :: ret

  if (ret /= NF90_NOERR) then
    write(6,*) '- FATAL ERROR.'
    write(6,*) nf90_strerror (ret)
    stop 999
  end if
end subroutine handle_err
