 program read_mask_test

! Test routine 'read_mask' by reading a sample C48 global
! uniform tile file.
!
! Author George Gayno NCEP/EMC

 use io_utils, only      : read_mask

 implicit none

 character(len=20)       :: merge_file

 integer, parameter      :: im=48
 integer, parameter      :: jm=48
 integer, parameter      :: chk_pts = 4
 integer                 :: i, j, ij

 real, parameter         :: EPSILON = 0.001
 real                    :: land_frac(im,jm)
 real                    :: lake_frac(im,jm)
 real                    :: slm(im,jm)
 real                    :: land_frac_chk(chk_pts)
 real                    :: lake_frac_chk(chk_pts)
 real                    :: slm_chk(chk_pts)

 type :: indices
  integer :: i
  integer :: j
 end type indices

 type(indices)           :: idx(chk_pts)

! The i/j indicies of the check points.

 data idx%i /1,29,47,48/
 data idx%j /1,25,23,48/

! The expected values of the check points.

 data land_frac_chk /0.58463, 0.8377, 0.0, 0.415374/
 data lake_frac_chk /0.0, 0.0, 1.0, 0.0/
 data slm_chk /1.0, 1.0, 0.0, 0.0/

 print*,"- Starting test of routine read_mask."

 merge_file = "./C48.mx500.tile1.nc"

! Initialize output fields to flag value.

 land_frac = -999.9
 lake_frac = -999.9
 slm       = -999.9

 call read_mask(merge_file,slm,land_frac,lake_frac,im,jm)

! All flag values should have been replaced. All fields
! should have values greater than zero.

 do j = 1, jm
 do i = 1, im
   if (land_frac(i,j) < 0.0) stop 3
   if (lake_frac(i,j) < 0.0) stop 5
   if (slm(i,j) < 0.0) stop 7
 enddo
 enddo
 
! Check some points against expected values.

 do ij = 1, chk_pts
   i = idx(ij)%i
   j = idx(ij)%j
   if (abs(land_frac(i,j)-land_frac_chk(ij)) > EPSILON) stop 12
   if (abs(lake_frac(i,j)-lake_frac_chk(ij)) > EPSILON) stop 15
   if (abs(slm(i,j)-slm_chk(ij)) > EPSILON) stop 18
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program read_mask_test
