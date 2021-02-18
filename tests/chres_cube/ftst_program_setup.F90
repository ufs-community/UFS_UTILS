! Unit test for chres_cube utility, input_data.F90.
!
! Ed Hartnett 2/16/21

program ftst_program_setup
  use esmf
  use netcdf
  use program_setup
  implicit none

  
  print*, "Starting test of program_setup."

  print*, "testing read_setup_namelist..."

  
  print*, "OK"
  call read_setup_namelist()
  print*, "SUCCESS!"
end program ftst_program_setup
