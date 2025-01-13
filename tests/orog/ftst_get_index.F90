 program test_get_index

! Unit test for routine "get_index", which finds the
! i/j location of the model point on the high-resolution
! mask and orography grids.
!
! Author George Gayno NCEP/EMC

 use orog_utils, only : get_index

 implicit none

 integer, parameter  :: imn=360*120
 integer, parameter  :: jmn=180*120
 integer, parameter  :: npts=4

 integer             :: i, ii, jst, jen, ilist(imn), numx

 real                :: lono(npts), lato(npts)
 real, parameter     :: delxn=360.0/float(imn)

 print*,"Starting test of get_index."

! The get_index routine contains a main 'if/elseif/else' statement. The test
! points are designed to check each branch.

! Point 1 - At equator. Western edge of grid cell at
! Greenwich. This tests the 'else' branch.

 print*,"Test point 1."

 lato(1) = 0.0; lono(1) = 0.0
 lato(2) = 1.0; lono(2) = 0.5
 lato(3) = 0.0; lono(3) = 1.0
 lato(4) = -1.0; lono(4) = 0.5

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 2  ! Starting 'j' index on hi-res grid.
 if (jen /= 10925) stop 4  ! Ending 'j' index on hi-res grid.
 if (numx /= 121)  stop 6  ! Number of 'i' indices on the hi-res grid.
 do i = 1, numx
   if (ilist(i) /= i) stop 8  ! List of 'i' indices on the hi-res grid.
 enddo

! Point 2 - At equator. Grid cell crosses Greenwich.
! This tests the 'elseif' branch.

 print*,"Test point 2."

 lato(1) = 0.0; lono(1) = 359.5
 lato(2) = 1.0; lono(2) = 0.0
 lato(3) = 0.0; lono(3) = 0.5
 lato(4) = -1.0; lono(4) = 0.0

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 12  ! Starting 'j' index on hi-res grid.
 if (jen /= 10925) stop 14  ! Ending 'j' index on hi-res grid.
 if (numx /= 121)  stop 16  ! Number of 'i' indices on the hi-res grid.
 ii = 1
 do i = 43141, 43200  ! East of greenwich
   if (ilist(ii) /= i) stop 18  ! List of 'i' indices on the hi-res grid.
   ii = ii + 1
 enddo
 do i = 1, 61  ! West of greenwich
   if (ilist(ii) /= i) stop 18
   ii = ii + 1
 enddo

! Point 3 - At equator. Grid cell centered at
! the dateline. This tests the 'else' branch.

 print*,"Test point 3."

 lato(1) = -1.0; lono(1) = 179.0
 lato(2) = 1.0; lono(2) = 179.0
 lato(3) = 1.0; lono(3) = 181.0
 lato(4) = -1.0; lono(4) = 181.0

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 22  ! Starting 'j' index on hi-res grid.
 if (jen /= 10925) stop 24  ! Ending 'j' index on hi-res grid.
 if (numx /= 241)  stop 26  ! Number of 'i' indices on the hi-res grid.
 ii = 1
 do i = 21481, 21721
   if (ilist(ii) /= i) stop 28 ! List of 'i' indices on the hi-res grid.
   ii = ii + 1
 enddo

! Point 4 - At North Pole. Grid cell centered at
! Greenwich. This tests the 'if' branch.

 print*,"Test point 4."

 lato(1) = 89.0; lono(1) = 359.5
 lato(2) = 90.0; lono(2) = 359.5
 lato(3) = 90.0; lono(3) = 0.5
 lato(4) = 89.0; lono(4) = 0.5

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 21476) stop 32  ! Starting 'j' index on hi-res grid.
 if (jen /= 21600) stop 34  ! Ending 'j' index on hi-res grid.
 if (numx /= 43200)  stop 36  ! Number of 'i' indices on the hi-res grid.
 do i = 1, numx
   if (ilist(i) /= i) stop 38 ! List of 'i' indices on the hi-res grid.
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program test_get_index
