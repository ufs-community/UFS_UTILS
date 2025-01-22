 program qc_orog_ramp

! Test routine "qc_orog_by_ramp", which adjusts
! the global terrain in the vicinity of Antarctica
! using 'RAMP' data.
!
! In OPS, the global data is 30-sec with dimensions
! 43200 x 21600. The RAMP data is 30-sec with dimension
! 43201 x 3601 and only covers Antarctica. For this test, 
! reduced versions of both grids are used:
! global - 9x7; RAMP 10x5.
!
! Author George Gayno NCEP/EMC

 use io_utils, only        : qc_orog_by_ramp

 implicit none

! Dimensions of the global grid.

 integer, parameter       :: imn = 9
 integer, parameter       :: jmn = 7

 integer                  :: i, j

! The terrain height (zavg) and land mask (zslm)
! of the global grid.

 integer                  :: zavg(imn,jmn)
 integer                  :: zslm(imn,jmn)
 
! The expected values for a successful test.

 integer                  :: zavg_expected(imn,jmn)
 integer                  :: zslm_expected(imn,jmn)

 print*,'- Starting test of qc_orog_by_ramp.'

! Initialize the global grid to all ocean.

 zslm = 0    ! water mask
 zavg = 0    ! sea level

! These global grid points are outside the 'ramp' grid,
! so they should not change.

 zslm_expected(:,6:7) = 0
 zavg_expected(:,6:7) = 0

! These global grid points are located within the 'ramp'
! grid. For this test, the first two rows of the 'ramp'
! data have non-zero terrain. So, rows 1 and 2 of the global
! grid will become land, located above sea level. Rows
! 3,4 and 5 of the RAMP data are ocean.  So, rows 3,4 and
! 5 of the global grid will remain ocean.

 zslm_expected(:,3:5) = 0 ! ocean mask
 zavg_expected(:,3:5) = 0 ! terrain height

 zslm_expected(:,1:2) = 1 ! becomes land

 zavg_expected(1,1)   = 5 ! acquires non-zero terrain.
 zavg_expected(2,1)   = 5
 zavg_expected(3,1)   = 5
 zavg_expected(4,1)   = 5
 zavg_expected(5,1)   = 5
 zavg_expected(6,1)   = 4
 zavg_expected(7,1)   = 4
 zavg_expected(8,1)   = 3
 zavg_expected(9,1)   = 3

 zavg_expected(1,2)   = 2
 zavg_expected(2,2)   = 2
 zavg_expected(3,2)   = 3
 zavg_expected(4,2)   = 2
 zavg_expected(5,2)   = 2
 zavg_expected(6,2)   = 2
 zavg_expected(7,2)   = 1
 zavg_expected(8,2)   = 1
 zavg_expected(9,2)   = 0 ! Note: this 'ramp' point has non-zero terrain of
                          ! 0.14, which rounds down to zero.

! Note: The path/name of the RAMP data is set in the routine.

 call qc_orog_by_ramp(imn, jmn, zavg, zslm)

 do i = 1, imn
 do j = 1, jmn
   if (zavg(i,j) /= zavg_expected(i,j)) stop 4
   if (zslm(i,j) /= zslm_expected(i,j)) stop 8
 enddo
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program qc_orog_ramp
