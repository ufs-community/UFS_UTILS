 program test_get_index

 use orog_utils, only : get_index

 implicit none

 integer, parameter  :: imn=360*120
 integer, parameter  :: jmn=180*120
 integer, parameter  :: npts=4

 integer             :: jst, jen, ilist(imn), numx

 real                :: lono(npts), lato(npts)
 real, parameter     :: delxn=360.0/imn

 print*,'hello world'

 lato(1) = 0.0; lono(1) = 0.0
 lato(2) = 1.0; lono(2) = 0.5
 lato(3) = 0.0; lono(3) = 1.0
 lato(4) = -1.0; lono(4) = 0.5

 ilist = -999

 call get_index(imn,jmn,npts,lonO,latO,delxn,jst,jen,ilist,numx)

 print*,jst,jen,numx

 print*,ilist(1:numx)


 end program test_get_index
