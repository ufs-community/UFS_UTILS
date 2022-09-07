! Unit test for filter_topo routine "read_namelist".
!
! Reads a sample namelist from file input.nml.
! If any namelist variable does not match expected values,
! the test fails.
!
! Author George Gayno 7/23/2021

 program readnml

 use utils

 implicit none

 print*, "Starting test of filter_topo routine read_namelist."
 print*, "Testing with file input.nml..."

 call read_namelist()

 if (trim(topo_file) /= "/dir1/dir2/orography.nc") stop 2
 if (trim(topo_field) /= "orog_filter") stop 4
 if (trim(mask_field) /= "landmask") stop 6
 if (trim(grid_file) /= "/dir1/dir2/dir3/mosaic.nc") stop 8
 if (zero_ocean) stop 10
 if (stretch_fac /= 2.0) stop 14
 if (res /= 96.0) stop 14
 if (.not. nested) stop 16
 if (grid_type /= 1) stop 18
 if (.not. regional) stop 20

 print*, "OK"
 print*, "SUCCESS!"

 end program readnml
