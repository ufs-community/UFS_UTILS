 program transpose

 use orog_utils, only : transpose_mask

 implicit none

 integer, parameter :: imn = 360
 integer, parameter :: jmn = 181

 integer(1)         :: mask(imn,jmn)
 integer            :: i, ii, j, jj

 print*,"Starting test of transpose routines."

! Transpose is from S to N to the NCEP standard N to S,
! and from the dateline to the NCEP standard Greenwich.

! Set up a one-degree global mask. Although mask is a
! yes/no flag, for this test set each row to the 
! latitude to simplify checking the answer.

 jj=0
 do j = -90, 90  ! row 1 is South Pole.
   jj = jj + 1
   mask(:,jj) = j
 enddo

 call transpose_mask(imn, jmn, mask)

 jj=0
 do j = 90, -90, -1  ! row 1 is North Pole.
   jj = jj + 1
   do i = 1, imn
     if (mask(i,jj) /= j) stop 2
   enddo
 enddo
 
! Now test the transpose in the E/W direction.
! Here, the East half of the domain is a flag value
! of minus 1 and the West half is plus 1.

 do i = 1, 180
   mask(i,:) = -1
 enddo
 do i = 181, 360
   mask(i,:) = +1
 enddo
 
 call transpose_mask(imn, jmn, mask)

! After the transpose, the East half should be plus 1 
! and the West half should be minus 1.

 do i = 1, 180
   do j = 1, jmn
    if (mask(i,j) /= 1) stop 4
   enddo
 enddo
 do i = 181, 360
   do j = 1, jmn
     if (mask(i,j) /= -1) stop 6
   enddo
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program transpose
