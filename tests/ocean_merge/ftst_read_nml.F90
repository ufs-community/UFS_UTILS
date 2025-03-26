! Unit test for the read_nml routine.
!
! Read a sample namelist and check the data in each
! variable against expected values.
!
 program read_namelist

 implicit none

 character(len=200) :: ocean_mask_dir
 character(len=200) :: lake_mask_dir
 character(len=200) :: out_dir
 character(len=10)  :: atmres,ocnres
 
 integer            :: binary_lake

 print*,"Call routine read_nml."

 call read_nml(ocean_mask_dir, lake_mask_dir, atmres,ocnres,out_dir,binary_lake)

 if (trim(ocean_mask_dir) /= "/ocean/mask/dir") stop 2
 if (trim(lake_mask_dir) /= "/lake/mask/dir") stop 4
 if (trim(atmres) /= "C96") stop 6
 if (trim(ocnres) /= "mx500") stop 8
 if (trim(out_dir) /= "/out/dir") stop 10
 if (binary_lake /= 1) stop 12

 print*, "OK"

 print*, "SUCCESS!"

 end program read_namelist
