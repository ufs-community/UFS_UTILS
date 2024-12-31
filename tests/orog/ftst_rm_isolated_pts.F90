 program rm_isolated_pts

! Unit test for subroutine remove_isolated_pts.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only    : remove_isolated_pts

 implicit none

 integer                :: im, jm, i, j, k

 real, allocatable      :: slm(:,:), oro(:,:), var(:,:), &
                           var4(:,:), oa(:,:,:), ol(:,:,:)      

 print*,"Starting test of remove_isolated_pts."

 im = 3
 jm = 3
 
 allocate (slm(im,jm))
 allocate (oro(im,jm))
 allocate (var(im,jm))
 allocate (var4(im,jm))
 allocate (oa(im,jm,4))
 allocate (ol(im,jm,4))

! Initialize grid to all ocean.

 slm = 0.0
 oro = 0.0
 var = 0.0
 var4 = 0.0
 oa = 0.0
 ol = 0.0

! This is an isolated island. The island should be
! removed (slm set to 0.0) and all other fields 
! should be the average of the surrounding points (which
! for this test case is zero.

 slm(2,2) = 1.0
 oro(2,2) = 50.0
 var(2,2) = 10.0
 var4(2,2) = 5.0

 oa(2,2,1) = -1.0
 oa(2,2,2) = -0.5
 oa(2,2,3) = 0.5
 oa(2,2,4) = 1.0

 ol(2,2,1) = 0.1
 ol(2,2,2) = 0.25
 ol(2,2,3) = 0.5
 ol(2,2,4) = 1.0

 call remove_isolated_pts(im,jm,slm,oro,var,var4,oa,ol)

 do j = 1, jm
 do i = 1, im
   if (slm(i,j) /= 0.0) stop 2
   if (oro(i,j) /= 0.0) stop 4
   if (var(i,j) /= 0.0) stop 6
   if (var4(i,j) /= 0.0) stop 8
   do k = 1, 4
     if (oa(i,j,k) /= 0.0) stop 10
     if (ol(i,j,k) /= 0.0) stop 12
   enddo
 enddo
 enddo

 deallocate (slm, oro, var, var4, oa, ol)

 print*,"OK"
 print*,"SUCCSSS"

 end program rm_isolated_pts
