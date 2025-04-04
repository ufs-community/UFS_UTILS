! Unit test for the readin_setup routine.
!
! Read a sample namelist and check the data in each
! variable against expected values.
!
 program read_namelist

 use  grids_IO, only     : grid_setup_type

 implicit none

 integer                  :: ierr

 type(grid_setup_type)    :: grid_setup_in, grid_setup_out

 print*,'Starting test of readin_setup.'

 open (12, file="./regrid.nml", iostat=ierr)
 if (ierr /= 0) stop 66

 print*,'Read input namelist.'
 call readin_setup(12, "input", grid_setup_in)

 if (trim(grid_setup_in%descriptor) /= "gau_inc") stop 2
 if (trim(grid_setup_in%dir) /= "./") stop 4
 if (trim(grid_setup_in%fname) /= "sfcincr_gsi") stop 6
 if (trim(grid_setup_in%dir_mask) /= "./") stop 8
 if (trim(grid_setup_in%fname_mask) /= "sfcincr_gsi") stop 10
 if (trim(grid_setup_in%mask_variable(1)) /= "soilsnow_mask") stop 12
 if (trim(grid_setup_in%dir_coord) /= "./") stop 14
 if (trim(grid_setup_in%fname_coord) /= "gaussian_scrip.nc") stop 16
 if (grid_setup_in%ires /= 768) stop 18
 if (grid_setup_in%jres /= 384) stop 20

 print*,'Read output namelist.'
 call readin_setup(12, "output", grid_setup_out)

 if (trim(grid_setup_out%descriptor) /= "fv3_rst") stop 32
 if (trim(grid_setup_out%dir) /= "./") stop 34
 if (trim(grid_setup_out%fname) /= "sfci") stop 36
 if (trim(grid_setup_out%dir_mask) /= "./") stop 38
 if (trim(grid_setup_out%fname_mask) /= "vegetation_type") stop 40
 if (trim(grid_setup_out%mask_variable(1)) /= "vegetation_type") stop 42
 if (trim(grid_setup_out%dir_coord) /= "./") stop 44
 if (trim(grid_setup_out%fname_coord) /= "NULL") stop 46
 if (grid_setup_out%ires /= 192) stop 48
 if (grid_setup_out%jres /= 192) stop 50

 close (12)

 print*, "OK"

 print*, "SUCCESS!"

 end program read_namelist
