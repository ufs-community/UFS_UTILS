 program minmax_test

! Unit test for routine minmax, which finds the
! minimum and maximum value of an array and
! the indices of the maximum.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : minmax

 implicit none

 character(len=8)   :: title

 integer, parameter :: im = 3
 integer, parameter :: jm = 2
 integer            :: imax, jmax

 real               :: a(im,jm)

 print*,"Starting test of minmax."

! Test array.

 a(1,1) = 3.
 a(2,1) = 4.
 a(3,1) = 2.
 a(1,2) = 1.
 a(2,2) = 4.
 a(3,2) = -1.

 title = 'test    '

! Call the routine to unit test.

 call minmax(im,jm,a,title,imax,jmax)

 if (imax /= 2 .or. jmax /= 2) stop 3

 print*,"OK"

 print*,"SUCCESS"

 end program minmax_test
