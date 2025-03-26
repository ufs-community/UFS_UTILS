 program read_gbl_mask

! Test routine "read_global_mask" using a reduced-size
! version (6 x 3 vs 43200 x 21600) of the umd land mask file. 
!
! Author George Gayno NCEP/EMC

 use io_utils, only : read_global_mask

 implicit none

 integer, parameter :: im=6
 integer, parameter :: jm=3

 integer            :: i, j

 integer(kind=1)    :: mask(im,jm)
 integer(kind=1)    :: mask_expected(im,jm)

! Note the routine flips the i and j directions.
 data mask_expected /1,1,1,0,0,1, &
                     1,1,1,0,0,0, &
                     1,1,1,0,0,0/

 print*,"- Begin test of read_global_mask"

 call read_global_mask(im, jm, mask)

 do j = 1, jm
 do i = 1, im
   if (mask(i,j) /= mask_expected(i,j)) stop 3
 enddo
 enddo

 print*,"OK"
 print*,"SUCCESS"

 end program read_gbl_mask
