 program read_gbl_orog

! Test routine "read_global_orog" using a reduced-size
! version (6 x 3 vs 43200 x 21600) of the gmted 2010 orog file.
!
! Author George Gayno NCEP/EMC

 use io_utils, only : read_global_orog

 implicit none

 integer, parameter :: im=6
 integer, parameter :: jm=3

 integer            :: i, j

 integer(kind=2)    :: glob(im,jm)
 integer(kind=2)    :: glob_expected(im,jm)

! Note the routine flips the i and j directions.
 data glob_expected /161,162,163,167,162,162, &
                     157,153,148,166,165,162, &
                     169,163,155,169,170,171/

 print*,"- Begin test of read_global_orog"

 call read_global_orog(im, jm, glob)

 do j = 1, jm
 do i = 1, im
   if (glob(i,j) /= glob_expected(i,j)) stop 3
 enddo
 enddo

 print*,"OK"
 print*,"SUCCESS"

 end program read_gbl_orog
