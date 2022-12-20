! Unit test for sfc_climo_gen utility, program_setup.f90.

program ftst_pgm_setup
  use program_setup
  use mpi
  implicit none

  integer :: my_rank, ierr
  
  call mpi_init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  if (my_rank .eq. 0) print*,"Starting test of program_setup."
  if (my_rank .eq. 0) print*,"Call routine read_setup_namelist."

  call read_setup_namelist(my_rank)

  if (trim(input_leaf_area_index_file) /= "leaf") stop 2
  if (trim(input_facsf_file) /= "facsf") stop 4
  if (trim(input_substrate_temperature_file) /= "substrate_temp") stop 6
  if (trim(input_maximum_snow_albedo_file) /= "maxsnow_alb") stop 8
  if (trim(input_snowfree_albedo_file) /= "snowfree_alb") stop 10
  if (trim(input_slope_type_file) /= "slope") stop 12
  if (trim(input_soil_type_file) /= "soil_type") stop 14
  if (trim(input_soil_color_file) /= "soil_color") stop 16
  if (trim(input_vegetation_type_file) /= "veg_type") stop 18
  if (trim(input_vegetation_greenness_file) /= "greenness") stop 20
  if (trim(mosaic_file_mdl) /= "mosaic.nc") stop 22
  if (trim(orog_dir_mdl) /= "./orog") stop 24
  if (trim(orog_files_mdl(1)) /= "oro.tile1.nc") stop 26
  if (trim(orog_files_mdl(2)) /= "oro.tile2.nc") stop 28
  if (trim(orog_files_mdl(3)) /= "oro.tile3.nc") stop 30
  if (trim(orog_files_mdl(4)) /= "oro.tile4.nc") stop 32
  if (trim(orog_files_mdl(5)) /= "oro.tile5.nc") stop 34
  if (trim(orog_files_mdl(6)) /= "oro.tile6.nc") stop 36
  if (trim(leaf_area_index_method) /= "bilinear") stop 38
  if (trim(maximum_snow_albedo_method) /= "conserve") stop 40
  if (trim(snowfree_albedo_method) /= "bilinear") stop 42
  if (trim(vegetation_greenness_method) /= "conserve") stop 44
  if (halo /= 4) stop 46

  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "SUCCESS!"
  
  call mpi_finalize(ierr)

end program ftst_pgm_setup
