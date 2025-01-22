 program rm_isolated_pts

! Unit test for subroutine remove_isolated_pts.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only    : remove_isolated_pts

 implicit none

 integer                :: im, jm, i, j, k

 real, parameter        :: EPSILON=0.001

 real, allocatable      :: slm(:,:), oro(:,:), var(:,:), &
                           var4(:,:), oa(:,:,:), ol(:,:,:)      

 real, allocatable      :: slm_expected(:,:), oro_expected(:,:), var_expected(:,:), &
                           var4_expected(:,:), oa_expected(:,:,:), ol_expected(:,:,:)      

 print*,"- Starting test of remove_isolated_pts."

! A 5x4 grid.

 im = 5
 jm = 4
 
 allocate (slm(im,jm))
 allocate (oro(im,jm))
 allocate (var(im,jm))
 allocate (var4(im,jm))
 allocate (oa(im,jm,4))
 allocate (ol(im,jm,4))

! Test Point 1.

! This is an isolated island. The island should be
! removed (slm set to 0.0) and all other fields 
! should be the average of the surrounding points (which
! for this test case is zero). All surrounding points
! should be unchanged.

 print*,'- Test point 1.'

! Initialize grid to all ocean.

 slm = 0.0    ! land/sea mask
 oro = 0.0    ! orography
 var = 0.0    ! standard deviation of orography
 var4 = 0.0   ! convexity
 oa = 0.0     ! orographic asymetry
 ol = 0.0     ! orographic length scale

! Initialize an island in the middle of the grid.

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

 allocate (slm_expected(im,jm))
 allocate (oro_expected(im,jm))
 allocate (var_expected(im,jm))
 allocate (var4_expected(im,jm))
 allocate (oa_expected(im,jm,4))
 allocate (ol_expected(im,jm,4))

! All points should be ocean (all fields zero).

 slm_expected = 0.0
 oro_expected = 0.0
 var_expected = 0.0
 var4_expected = 0.0
 oa_expected = 0.0
 ol_expected = 0.0

 call remove_isolated_pts(im,jm,slm,oro,var,var4,oa,ol)

 do j = 1, jm
 do i = 1, im
   if (abs(slm(i,j)-slm_expected(i,j)) > EPSILON) stop 2
   if (abs(oro(i,j)-oro_expected(i,j)) > EPSILON) stop 4
   if (abs(var(i,j)-var_expected(i,j)) > EPSILON) stop 6
   if (abs(var4(i,j)-var4_expected(i,j)) > EPSILON) stop 8
   do k = 1, 4
     if (abs(oa(i,j,k)-oa_expected(i,j,k)) > EPSILON) stop 10
     if (abs(ol(i,j,k)-ol_expected(i,j,k)) > EPSILON) stop 12
   enddo
 enddo
 enddo

 deallocate (slm, oro, var, var4, oa, ol)
 deallocate (slm_expected, oro_expected, var_expected, var4_expected, oa_expected, ol_expected)

! Test Point 2

! Now remove an isolated water point. The point should be changed
! to land (slm=1.0) and all other fields should be the average
! of the eight surrounding points.

 print*,'- Test point 2.'

! A 3x3 grid.

 im = 3
 jm = 3
 
 allocate (slm(im,jm))
 allocate (oro(im,jm))
 allocate (var(im,jm))
 allocate (var4(im,jm))
 allocate (oa(im,jm,4))
 allocate (ol(im,jm,4))

 slm = 1.0      ! Initialize grid to all land.
 slm(2,2) = 0.0 ! Set interior point to water.

 oro(1,1) = 5.0
 oro(2,1) = 10.0
 oro(3,1) = 17.0
 oro(1,2) = 22.0
 oro(2,2) = 0.0  ! water point.
 oro(3,2) = 100.0
 oro(1,3) = 34.0
 oro(2,3) = 55.0
 oro(3,3) = 68.0

 var = 0.5 * oro
 var4 = 0.25 * oro

 ol(1,1,1) = 0.1
 ol(2,1,1) = 0.2
 ol(3,1,1) = 0.3
 ol(1,2,1) = 0.4
 ol(2,2,1) = 0.0 ! water point
 ol(3,2,1) = 0.6
 ol(1,3,1) = 0.7
 ol(2,3,1) = 0.8
 ol(3,3,1) = 0.9

 do j = 1, jm
 do i = 1, im
   ol(i,j,2) = ol(i,j,1) * .75
   ol(i,j,3) = ol(i,j,1) * .5
   ol(i,j,4) = ol(i,j,1) * .25
 enddo
 enddo

 oa = -(ol)

! All points should be land. The former isolated ocean
! point should have values that are the average of the
! surrounding land points. All other points should remain
! unchanged.

 allocate (slm_expected(im,jm))
 allocate (oro_expected(im,jm))
 allocate (var_expected(im,jm))
 allocate (var4_expected(im,jm))
 allocate (oa_expected(im,jm,4))
 allocate (ol_expected(im,jm,4))

 slm_expected = slm
 oro_expected = oro
 var_expected = var
 var4_expected = var4
 ol_expected = ol
 oa_expected = oa

! Former ocean point should be average of the surrounding
! ocean points.

 slm_expected(2,2) = 1.0 ! land point
 oro_expected(2,2) = 38.875
 var_expected(2,2) = 0.5 * oro_expected(2,2)
 var4_expected(2,2) = 0.25 * oro_expected(2,2)
 ol_expected(2,2,1) = 0.5
 ol_expected(2,2,2) = 0.375
 ol_expected(2,2,3) = 0.25
 ol_expected(2,2,4) = 0.125
 oa_expected(2,2,1) = -0.5
 oa_expected(2,2,2) = -0.375
 oa_expected(2,2,3) = -0.25
 oa_expected(2,2,4) = -0.125

 call remove_isolated_pts(im,jm,slm,oro,var,var4,oa,ol)

 do j = 1, jm
 do i = 1, im
   if (abs(slm(i,j)-slm_expected(i,j)) > EPSILON) stop 22
   if (abs(oro(i,j)-oro_expected(i,j)) > EPSILON) stop 24
   if (abs(var(i,j)-var_expected(i,j)) > EPSILON) stop 26
   if (abs(var4(i,j)-var4_expected(i,j)) > EPSILON) stop 28
   do k = 1, 4
     if (abs(oa(i,j,k)-oa_expected(i,j,k)) > EPSILON) stop 30
     if (abs(ol(i,j,k)-ol_expected(i,j,k)) > EPSILON) stop 32
   enddo
 enddo
 enddo

 deallocate (slm, oro, var, var4, oa, ol)
 deallocate (slm_expected, oro_expected, var_expected, var4_expected, oa_expected, ol_expected)

 print*,"OK"
 print*,"SUCCESS"

 end program rm_isolated_pts
