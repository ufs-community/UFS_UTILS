 program transpose

! Unit test for routines transpose_mask and transpose_orog.
! They are essentially the same routine - one takes 
! one byte integer data, and one takes two byte integer
! data.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : transpose_mask, transpose_orog

 implicit none

 integer, parameter :: imn = 360
 integer, parameter :: jmn = 181

 integer(1)         :: mask(imn,jmn)
 integer(2)         :: mask2(imn,jmn)
 integer            :: i, j, jj

 print*,"Starting test of transpose routines."

! Transpose is from S to N to the NCEP standard N to S,
! and, in the E/W direction, from the dateline to the 
! NCEP standard Greenwich.

! Test N/S transpose. Set up a one-degree global mask.
! Although mask is a yes/no flag, for this test set each
! row to the latitude to simplify checking the answer.

 jj=0
 do j = -90, 90  ! row 1 is South Pole.
   jj = jj + 1
   mask(:,jj) = j
   mask2(:,jj) = j
 enddo

 print*,"Call transpose routines to test N/S flip."

 call transpose_mask(imn, jmn, mask)
 call transpose_orog(imn, jmn, mask2)

 jj=0
 do j = 90, -90, -1  ! row 1 is North Pole.
   jj = jj + 1
   do i = 1, imn
     if (mask(i,jj) /= j) stop 2
     if (mask2(i,jj) /= j) stop 3
   enddo
 enddo
 
! Now test the transpose in the E/W direction.
! Here, the East half of the domain is a flag value
! of minus 1 and the West half is plus 1.

 do i = 1, 180
   mask(i,:) = -1
   mask2(i,:) = -1
 enddo
 do i = 181, 360
   mask(i,:) = +1
   mask2(i,:) = +1
 enddo
 
 print*,"Call transpose routines to test E/W flip."

 call transpose_mask(imn, jmn, mask)
 call transpose_orog(imn, jmn, mask2)

! After the transpose, the East half should be plus 1 
! and the West half should be minus 1.

 do i = 1, 180
   do j = 1, jmn
    if (mask(i,j) /= 1) stop 4
    if (mask2(i,j) /= 1) stop 5
   enddo
 enddo
 do i = 181, 360
   do j = 1, jmn
     if (mask(i,j) /= -1) stop 6
     if (mask2(i,j) /= -1) stop 7
   enddo
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program transpose
