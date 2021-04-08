 program test_search_util

 use mpi
 use esmf
 use search_util

 implicit none

 integer, parameter :: idim = 3
 integer, parameter :: jdim = 3
 integer, parameter :: tile = 1

 integer :: field_num, ierr
 integer(esmf_kind_i8) :: mask(idim,jdim)

 real(esmf_kind_r8) :: field(idim,jdim)
 real(esmf_kind_r8) :: latitude(idim,jdim)
 real(esmf_kind_r8) :: terrain_land(idim,jdim)
 real(esmf_kind_r8) :: soilt_climo(idim,jdim)

 data mask /0, 0, 0,  1, 1, 1, 0, 0, 0/

 data latitude /-30.0, -30.0, -30.0, 0., 0., 0., 25.0, 25.0, 25.0/

 data terrain_land /0., 0., 0., 50.0, 75.0, 80.0, 0., 0., 0./

 data soilt_climo /285.0, 285.0, 285.0, 0., 0., 0., 278.0, 278.0, 278.0/

 data field/0., 0., 0.,  0., 0., 0.,  0., 0., 0./

 call mpi_init(ierr)
 
 print*,'hello world'

 field_num = 0
 
 call search (field, mask, idim, jdim, tile, field_num, latitude, terrain_land, &
       soilt_climo)

 call mpi_finalize(ierr)

 end program test_search_util
