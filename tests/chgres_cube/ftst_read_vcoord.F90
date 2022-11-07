 program vcoord

! Test the 'read_vcoord_info' routine, which reads
! the 'global_hyblev' files. These files contain
! the definition of the model vertical hybrid levels.
! The levels are defined by two coefficients 'ak'
! and 'bk' as follows:
!
! pressure(level) = ak + (bk * surface pressure)
!
! The 'read_vcoord_info' routine returns the
! number of vertical levels, the number of level
! interfaces (# of levels plus 1), the number of 
! level coordinates (two for 'ak' and 'bk') and
! the 'ak' and 'bk' values. All are compared to
! expected values. The sample 'global_hyblev'
! file is set by variable "vcoord_file_target_grid'.
!
! @author George Gayno

 use atmosphere_target_data, only : vcoord_target, &
                             nvcoord_target, &
                             lev_target, &
                             levp1_target

 use atmosphere, only : read_vcoord_info

 use program_setup, only : vcoord_file_target_grid

 implicit none

 integer :: j

 integer, parameter :: LEV_TARGET_EXPECTED=28    ! number of levels.
 integer, parameter :: LEVP1_TARGET_EXPECTED=29  ! number of level interfaces.
 integer, parameter :: NVCOORD_TARGET_EXPECTED=2 ! number of level coordinates.

 real :: VCOORD_AK_EXPECTED(LEVP1_TARGET_EXPECTED) ! 'ak' values for each interface.
 real :: VCOORD_BK_EXPECTED(LEVP1_TARGET_EXPECTED) ! 'bk' values for each interface.

 data VCOORD_AK_EXPECTED /.000,    .009, 11.635,  &
                           86.457, 292.577,  701.453, &
                           1381.526, 2389.205, 3757.142, &
                           5481.144, 7508.532,  9731.807, & 
                           11991.428, 14089.535, 15812.926, &
                           16959.994, 17364.658, 16912.130, & 
                           15613.564, 13671.427, 11343.654, &
                           8913.767, 6678.608, 4844.036, &
                           3376.516, 2210.979,  1290.533, &
                           566.898,  0.000/

 data VCOORD_BK_EXPECTED /1.00000000, .98872586, .97440184, &
                          .95587239,  .93174961, .90058087, &
                          .86097487,  .81178485, .75234710, &
                          .68274682,  .60405491, .51845667, &
                          .42919536,  .34029321, .25608430, &
                          .18066704,  .11741761, .06867500, &
                          .03495010,  .01432627, .00393276, &
                          .00037868,   0.0,       0.0, &
                          0.0,         0.0,       0.0, &
                          0.0,         0.0 /

 print*,'Starting test of read_vcoord_info routine'

 vcoord_file_target_grid="./data/global_hyblev.l28.txt"
 
 call read_vcoord_info

 if (lev_target /= LEV_TARGET_EXPECTED) stop 2
 if (levp1_target /= LEVP1_TARGET_EXPECTED) stop 4
 if (nvcoord_target /= NVCOORD_TARGET_EXPECTED) stop 6

 do j = 1, levp1_target
   if (vcoord_target(j,1) /= VCOORD_AK_EXPECTED(j)) stop 8
   if (vcoord_target(j,2) /= VCOORD_BK_EXPECTED(j)) stop 10
 enddo

 print*,"OK"

 print*,"SUCCESS!"

 end program vcoord
