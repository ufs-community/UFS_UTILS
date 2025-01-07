 program read_model_dims

 use io_utils, only     : read_mdl_dims

 implicit none

 character(len=19)    :: mdl_grid_file

 integer             :: im, jm

 print*,"- Begin test of routine read_mdl_dims."

 mdl_grid_file="./C12_grid.tile1.nc"

 call read_mdl_dims(mdl_grid_file,im,jm)

 print*,'im/jm ',im,jm

 print*,"OK"

 print*,"SUCCESS"

 end program read_model_dims
