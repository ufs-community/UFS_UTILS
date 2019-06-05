module atmdata_type

 use esmf

 implicit none

 type atmdata
  real(esmf_kind_r8), pointer         :: var(:,:,:)
 end type atmdata
 
 end module atmdata_type