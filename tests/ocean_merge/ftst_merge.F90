! Unit test for the merge routine.
! 
! Test several combinations of lake and ocean attributes
! and check the 'merged' value for correctness.
 program ftst_merge

 implicit none

 integer, parameter   :: lon = 1
 integer, parameter   :: lat = 1

 integer              :: binary_lake

 real, parameter      :: epsilon=0.00001
 real                 :: lat2d(lon,lat)
 real                 :: ocn_frac(lon,lat)
 real                 :: lake_frac(lon,lat)
 real                 :: lake_depth(lon,lat)
 real                 :: land_frac(lon,lat)
 real                 :: slmsk(lon,lat)

 print*,"Begin test of merge routine"

! Test point 1
! Some ocean. No lake. Ocean info retained.

 print*,'Test point 1.'

 binary_lake = 0         ! Keep fractional lake.
 lat2d(1,1) = 30.0       ! Point at 30N.
 ocn_frac(1,1) = .75     ! .75 ocean/.25 land
 lake_frac(1,1) = 0.0    ! No lake.
 lake_depth(1,1) = 0.0   ! No lake.
 land_frac(1,1) = -99.   ! Based on ocn_frac, should be .25.
 slmsk(1,1) = -99.       ! Should be zero (nint of land_frac)

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 2
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 2
 if ( abs(land_frac(1,1) - .25) > epsilon) stop 2
 if ( abs(slmsk(1,1) - 0.) > epsilon) stop 2

! Test point 2
! No ocean. Some lake. Lake info retained.

 print*,'Test point 2.'

 binary_lake = 0          ! Keep fractional lake.
 lat2d(1,1) = 30.0        ! Point at 30N.
 ocn_frac(1,1) = 0.0      ! 0 ocean/1 land.
 lake_frac(1,1) = .75     ! Some lake.
 lake_depth(1,1) = 100.0  ! Lake depth.
 land_frac(1,1) = -99.    ! Should be .25
 slmsk(1,1) = -99.        ! Should be zero (nint of land_frac).

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - .75) > epsilon) stop 4
 if ( abs(lake_depth(1,1) - 100.) > epsilon) stop 4
 if ( abs(land_frac(1,1) - .25) > epsilon) stop 4
 if ( abs(slmsk(1,1) - 0.) > epsilon) stop 4

! Test point 3
! Some lake and some ocean. Ocean should dominate. Lake
! should be removed.

 print*,'Test point 3.'

 binary_lake = 0           ! Keep fractional lake.
 lat2d(1,1) = 30.0         ! Point at 30N.
 ocn_frac(1,1) = .45       ! Some ocean.
 lake_frac(1,1) = .75      ! Some lake.
 lake_depth(1,1) = 100.0   ! Lake depth
 land_frac(1,1) = -99.     ! Should be 1 minus ocn_frac
 slmsk(1,1) = -99.         ! Should be 1 (nint of land_frac).

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 6
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 6
 if ( abs(land_frac(1,1) - .55) > epsilon) stop 6
 if ( abs(slmsk(1,1) - 1.) > epsilon) stop 6

! Test point 4
! No ocean. Some lake. Lake will be retainted. Lake has 
! a missing depth that should be given a default value.

 print*,'Test point 4.'

 binary_lake = 0           ! Keep fractional lake.
 lat2d(1,1) = 30.0         ! Point at 30N.
 ocn_frac(1,1) = 0.0       ! No ocean.
 lake_frac(1,1) = .75      ! Some lake.
 lake_depth(1,1) = -9.     ! Lake depth is missing.
 land_frac(1,1) = -99.     ! Should be 1 minus lake_frac
 slmsk(1,1) = -99.         ! Should be 0 (nint of land_frac).

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - .75) > epsilon) stop 8
 if ( abs(lake_depth(1,1) - 10.) > epsilon) stop 8
 if ( abs(land_frac(1,1) - .25) > epsilon) stop 8
 if ( abs(slmsk(1,1) - 0.) > epsilon) stop 8

! Test point 5
! Some ocean (but very small percentage). No lake. 
! The ocean % is below the min threshold, so it will be
! removed and the point becomes all land.

 print*,'Test point 5.'

 binary_lake = 0          ! Keep fractional lake.
 lat2d(1,1) = 30.0        ! Point at 30N.
 ocn_frac(1,1) = 1.e-6    ! Very small % of ocean.
 lake_frac(1,1) = 0.0     ! No lake.
 lake_depth(1,1) = 0.0    ! No lake.
 land_frac(1,1) = -99.    ! Should be 1.
 slmsk(1,1) = -99.        ! Should be 1.

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 10
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 10
 if ( abs(land_frac(1,1) - 1.) > epsilon) stop 10
 if ( abs(slmsk(1,1) - 1.) > epsilon) stop 10

! Test point 6
! Almost all ocean. Some lake. Since the ocean %
! is very close to 1, it will be rounded up to
! exactly 1. Since ocean dominates, the lake
! will be removed.

 print*,'Test point 6.'

 binary_lake = 0           ! Keep fractional lake.
 lat2d(1,1) = 30.0         ! Point at 30N
 ocn_frac(1,1) = 0.99999   ! Very high % of ocean.
 lake_frac(1,1) = 0.24     ! Some lake.
 lake_depth(1,1) = 50.0    ! Lake depth.
 land_frac(1,1) = -99.     ! Should be 0.
 slmsk(1,1) = -99.         ! Should be 0.

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 12
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 12
 if ( abs(land_frac(1,1) - 0.) > epsilon) stop 12
 if ( abs(slmsk(1,1) - 0.) > epsilon) stop 12

! Test point 7
! No ocean. Some lake, but it is near Antarctica.
! Lakes near Antarctica are to be removed. Point
! becomes all land.

 print*,'Test point 7.'

 binary_lake = 0           ! Keep fractional lake.
 lat2d(1,1) = -70.0        ! Point at 70S.
 ocn_frac(1,1) = 0.0       ! No ocean.
 lake_frac(1,1) = .75      ! Some lake.
 lake_depth(1,1) = 100.0   ! Lake depth.
 land_frac(1,1) = -99.     ! Should be 1.
 slmsk(1,1) = -99.         ! Should be 1.

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 14
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 14
 if ( abs(land_frac(1,1) - 1.) > epsilon) stop 14
 if ( abs(slmsk(1,1) - 1.) > epsilon) stop 14

! Test point 8
! No ocean. Some lake. Assume binary yes/no lake.
! In this case, the lake fraction will become zero and
! the land fraction will become 1.

 print*,'Test point 8.'

 binary_lake = 1           ! Keep fractional lake.
 lat2d(1,1) = 30.0         ! Point at 30N.
 ocn_frac(1,1) = 0.0       ! No ocean.
 lake_frac(1,1) = .15      ! Less than 50% lake.
 lake_depth(1,1) = 100.0   ! Lake depth.
 land_frac(1,1) = -99.     ! Should be 1.
 slmsk(1,1) = -99.         ! Should be 1.

 call merge(lon, lat, binary_lake, lat2d, ocn_frac, &
            lake_frac, lake_depth, land_frac, slmsk)

 if ( abs(lake_frac(1,1) - 0.) > epsilon) stop 16
 if ( abs(lake_depth(1,1) - 0.) > epsilon) stop 16
 if ( abs(land_frac(1,1) - 1.) > epsilon) stop 16
 if ( abs(slmsk(1,1) - 1.) > epsilon) stop 16

 print*, "OK"

 print*, "SUCCESS!"

 end program ftst_merge
