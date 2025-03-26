 program find_nearest_pole_pts

! Unit test for routine find_nearest_pole_points.
!
! The routine is passed the supergrid indices of
! the pole (from the 'grid' files) and returns the 
! nearest pole points (using a logical array) on 
! the model tile.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only     : find_nearest_pole_points

 implicit none

 integer, parameter      :: im=48
 integer, parameter      :: jm=48

 integer                 :: i, j
 integer                 :: i_north_pole, j_north_pole
 integer                 :: i_south_pole, j_south_pole

 logical                 :: is_north_pole(im,jm)
 logical                 :: is_south_pole(im,jm)

 print*,"Starting test of find_nearest_pole_points routine."

! Test 1 - C48 uniform tile containing north pole.

 print*,"Test number 1."

 i_north_pole = 49 ! Supergrid index
 j_north_pole = 49
 i_south_pole = 0
 j_south_pole = 0

 call find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 do j = 1, im
 do i = 1, jm
   if((i == 24 .and. j == 24) .or. &
      (i == 24 .and. j == 25) .or. &
      (i == 25 .and. j == 24) .or. &
      (i == 25 .and. j == 25)) then
     if (.not.is_north_pole(i,j)) stop 2
   else
     if (is_north_pole(i,j)) stop 4
   endif
   if (is_south_pole(i,j)) stop 8
 enddo
 enddo

! Test 2 - C48 uniform tile containing south pole.

 print*,"Test number 2."

 i_north_pole = 0
 j_north_pole = 0
 i_south_pole = 49 ! Supergrid index.
 j_south_pole = 49

 call find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 do j = 1, im
 do i = 1, jm
   if((i == 24 .and. j == 24) .or. &
      (i == 24 .and. j == 25) .or. &
      (i == 25 .and. j == 24) .or. &
      (i == 25 .and. j == 25)) then
     if (.not.is_south_pole(i,j)) stop 12
   else
     if (is_south_pole(i,j)) stop 14
   endif
   if (is_north_pole(i,j)) stop 18
 enddo
 enddo

! Test 3 - C48 uniform tile containing no pole.

 print*,"Test number 3."

 i_north_pole = 0  ! Zero indicates no pole.
 j_north_pole = 0
 i_south_pole = 0
 j_south_pole = 0

 call find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 do j = 1, im
 do i = 1, jm
   if (is_south_pole(i,j)) stop 24
   if (is_north_pole(i,j)) stop 26
 enddo
 enddo

! Test 4 - C48 stretched grid tile containing south pole.

 print*,"Test number 4."

 i_north_pole = 0
 j_north_pole = 0
 i_south_pole = 10 ! Supergrid index.
 j_south_pole = 49

 call find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 do j = 1, im
 do i = 1, jm
   if((i == 5 .and. j == 24) .or. &
      (i == 5 .and. j == 25)) then
     if (.not.is_south_pole(i,j)) stop 32
   else
     if (is_south_pole(i,j)) stop 34
   endif
   if (is_north_pole(i,j)) stop 38
 enddo
 enddo

! Test 5 - C48 stretched grid tile containing north pole.

 print*,"Test number 5."

 i_north_pole = 62  ! Supergrid index.
 j_north_pole = 49
 i_south_pole = 0
 j_south_pole = 0

 call find_nearest_pole_points(i_north_pole, j_north_pole, &
      i_south_pole, j_south_pole, im, jm, is_north_pole, &
      is_south_pole)

 do j = 1, im
 do i = 1, jm
   if((i == 31 .and. j == 24) .or. &
      (i == 31 .and. j == 25)) then
     if (.not.is_north_pole(i,j)) stop 32
   else
     if (is_north_pole(i,j)) stop 34
   endif
   if (is_south_pole(i,j)) stop 38
 enddo
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program find_nearest_pole_pts
