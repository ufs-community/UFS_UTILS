! Unit test for filter_topo routine "fill_regional_halo".
!
! Author Ratko Vasic

 program fill_halo

 use utils

 implicit none

 integer :: halo,ix,jx,i,j

 real, allocatable :: testdata(:,:,:)

 print*, "Starting test of filter_topo routine fill_regional_halo"

 halo = 3
 ix=6
 jx=9

 allocate(testdata(1-halo:ix+halo,1-halo:jx+halo,1))

! Initialize whole domain (including halo points) to zero
 testdata(:,:,:)=0.

! Initialize inner domain
 do j=1,jx
 do i=1,ix
   testdata(i,j,1)=10*i+j
 enddo
 enddo

!!! do j=1-halo,jx+halo
!!!  write(0,101) testdata(1-halo:ix+halo,j,1)
!!! enddo
!!!101 format(12(f9.5,x))

! Fill halo points
 call fill_regional_halo(testdata, halo)

!!! print*
!!! do j=1-halo,jx+halo
!!!  write(0,101) testdata(1-halo:ix+halo,j,1)
!!! enddo

 if(testdata(-2,-2,1)==  5.5 .and. testdata(-2,12,1)== 14.5 .and. &
    testdata( 9,-2,1)== 65.5 .and. testdata( 9,12,1)== 74.5) then
  print*, "OK"
  print*, "SUCCESS!"
 else
  stop 22
 endif

 deallocate(testdata)

 end program fill_halo
