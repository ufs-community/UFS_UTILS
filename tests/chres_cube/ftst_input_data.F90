! Unit test for chres_cube utility, input_data.F90.
!
! Ed Hartnett 2/16/21

program ftst_input_data
  use esmf
  use netcdf
  use model_grid, only : input_grid, i_input, j_input,  &
       ip1_input, jp1_input, num_tiles_input_grid, &
       latitude_input_grid, longitude_input_grid, inv_file
  use input_data, only : read_fv3_grid_data_netcdf
  implicit none

  integer :: tile
  integer :: lsoil_input=4
  real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
  
  print*, "Starting test of input_data."

  print*, "testing read_fv3_grid_data_netcdf..."
  allocate(data_one_tile(i_input,j_input))
  call read_fv3_grid_data_netcdf('slc', tile, i_input, j_input, &
       lsoil_input, sfcdata=data_one_tile)
  deallocate(data_one_tile)

  print*, "OK"

  print*, "SUCCESS!"
end program ftst_input_data
