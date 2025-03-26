 program read_model_grid_file

! Test routine 'read_mdl_grid_file' by reading a 
! C12 'grid' file from a global uniform grid.

! Author George Gayno NCEP/EMC

 use io_utils, only     : read_mdl_grid_file

 implicit none

 character(len=19)    :: mdl_grid_file

 integer, parameter   :: im = 12
 integer, parameter   :: jm = 12
 integer              :: i, j
 
 logical              :: is_north_pole(im,jm)
 logical              :: is_south_pole(im,jm)

 real, parameter      :: EPSILON=0.01
 real                 :: geolat(im,jm)
 real                 :: geolat_c(im+1,jm+1)
 real                 :: geolon(im,jm)
 real                 :: geolon_c(im+1,jm+1)
 real                 :: dx(im,jm), dy(im,jm)

 print*,"- Begin test of routine read_mdl_grid_file."

 mdl_grid_file="./C12_grid.tile1.nc"

! Initialize all output variables to flag values.

 is_north_pole=.true.
 is_south_pole=.true.

 geolat   = -999.9
 geolat_c = -999.9
 geolon   = -999.9
 geolon_c = -999.9
 dx       = -999.9
 dy       = -999.9

 call read_mdl_grid_file(mdl_grid_file, im, jm, &
             geolon, geolon_c, geolat, geolat_c, dx, dy, &
             is_north_pole, is_south_pole)
 
! Check values at selected points.

 if (abs(geolon_c(1,1)-305.0) > EPSILON) stop 2
 if (abs(geolon_c(13,13)-35.0) > EPSILON) stop 4
 if (abs(geolat_c(13,1)-(-35.2644)) > EPSILON) stop 6
 if (abs(geolat_c(1,13)-35.2644) > EPSILON) stop 8
 if (abs(geolon(5,5)-337.656) > EPSILON) stop 10
 if (abs(geolon(10,10)-17.9245) > EPSILON) stop 12
 if (abs(geolat(7,3)-(-27.84)) > EPSILON) stop 14
 if (abs(geolat(2,9)-16.8678) > EPSILON) stop 16
 if (abs(dx(1,1)-631596.076) > EPSILON) stop 18
 if (abs(dx(12,12)-631596.076) > EPSILON) stop 20

! There is no pole on this tile, so the pole variables
! should be .false. The routine sets dx equal to dy,
! so they should match.

 do j = 1, jm
 do i = 1, im
   if (abs(dx(i,j)-dy(i,j)) > EPSILON) stop 22 
   if (is_north_pole(i,j)) stop 24
   if (is_south_pole(i,j)) stop 26
 enddo
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program read_model_grid_file
