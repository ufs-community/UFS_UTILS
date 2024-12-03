 program find_nearest_pole_pts

 use orog_utils, only     : find_nearest_pole_points

 implicit none

 integer, parameter      :: im=48
 integer, parameter      :: jm=48

 integer                 :: i, j
 integer                 :: i_north_pole, j_north_pole
 integer                 :: i_south_pole, j_south_pole

 logical                 :: is_north_pole(im,jm)
 logical                 :: is_south_pole(im,jm)

! Test 1 - C48 uniform tile containing north pole.

 i_north_pole = 49 ! supergrid index
 j_north_pole = 49 ! uniform grid
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

 i_north_pole = 0
 j_north_pole = 0
 i_south_pole = 49
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

 i_north_pole = 0
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

 i_north_pole = 0
 j_north_pole = 0
 i_south_pole = 10
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

 end program find_nearest_pole_pts
