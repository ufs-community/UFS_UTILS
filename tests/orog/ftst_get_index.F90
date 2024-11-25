 program test_get_index

! Unit test for routine get_index, which finds the
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

! Point 1 - At equator. Western edge of grid cell at
! Greenwich.

 lato(1) = 0.0; lono(1) = 0.0
 lato(2) = 1.0; lono(2) = 0.5
 lato(3) = 0.0; lono(3) = 1.0
 lato(4) = -1.0; lono(4) = 0.5

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 2
 if (jen /= 10925) stop 4
 if (numx /= 121)  stop 6
 do i = 1, numx
   if (ilist(i) /= i) stop 8
 enddo

! Point 2 - At equator. Grid cell centered at
! Greenwich.

 lato(1) = 0.0; lono(1) = -0.5
 lato(2) = 1.0; lono(2) = 0.0
 lato(3) = 0.0; lono(3) = 0.5
 lato(4) = -1.0; lono(4) = 0.0

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 12
 if (jen /= 10925) stop 14
 if (numx /= 121)  stop 16
 ii = 1
 do i = -59, 61, 1
   if (ilist(ii) /= i) stop 18
   ii = ii + 1
 enddo

! Point 3 - At equator. Grid cell centered at
! the dateline.

 lato(1) = -1.0; lono(1) = 179.0
 lato(2) = 1.0; lono(2) = 179.0
 lato(3) = 1.0; lono(3) = 181.0
 lato(4) = -1.0; lono(4) = 181.0

 ilist = -999

 call get_index(imn,jmn,npts,lono,lato,delxn,jst,jen,ilist,numx)

 if (jst /= 10676) stop 22
 if (jen /= 10925) stop 24
 if (numx /= 241)  stop 26
 ii = 1
 do i = 21481, 21721
   if (ilist(ii) /= i) stop 28
   ii = ii + 1
 enddo

 print*,"OK"

 print*,"SUCCESS"

 end program test_get_index
