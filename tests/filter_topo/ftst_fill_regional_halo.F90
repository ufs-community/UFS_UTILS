! Unit test for filter_topo routine "fill_regional_halo".
!
! Author Ratko Vasic

 program fill_halo

 use utils

 implicit none

 integer :: halo

 real, allocatable :: testdata(:,:,:)

 print*, "Starting test of filter_topo routine fill_regional_halo"

 call fill_regional_halo(testdata, halo)


 print*, "OK"
 print*, "SUCCESS!"

 end program fill_halo
