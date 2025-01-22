 program read_model_dims

! Test routine "read_mdl_dims", which gets the
! i/j dimensions of the model tile from the
! 'grid' file. This test uses a global uniform
! C12 'grid' file. The dimensions returned from
! the routine should be i=12, j=12.
!
! Author George Gayno NCEP/EMC

 use io_utils, only     : read_mdl_dims

 implicit none

 character(len=19)    :: mdl_grid_file

 integer             :: im, jm

 print*,"- Begin test of routine read_mdl_dims."

 mdl_grid_file="./C12_grid.tile1.nc"

 call read_mdl_dims(mdl_grid_file,im,jm)

 if (im /= 12) stop 3
 if (jm /= 12) stop 6

 print*,"OK"

 print*,"SUCCESS"

 end program read_model_dims
