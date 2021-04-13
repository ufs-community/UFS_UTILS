 program test_search_util

 use mpi
 use esmf
 use search_util

 implicit none

! Test uses a simple 3x3 fv3 tile.

 integer, parameter :: idim = 3
 integer, parameter :: jdim = 3
 integer, parameter :: tile = 1

 integer :: field_num, ierr, i
 integer(esmf_kind_i8) :: mask(idim,jdim)

 real(esmf_kind_r8) :: field(idim,jdim)
 real(esmf_kind_r8) :: field_updated(idim,jdim)
 real(esmf_kind_r8) :: latitude(idim,jdim)
 real(esmf_kind_r8) :: terrain_land(idim,jdim)
 real(esmf_kind_r8) :: soilt_climo(idim,jdim)

! These variables are used to test the 'default'
! search routine option (i.e., when the search
! fails, a default must be used)

 integer, parameter :: num_default_tests = 21
 integer :: default_field_num(num_default_tests)
 real(esmf_kind_r8) :: default_field_val(num_default_tests)

! Definition of the mask. The '1' indicates an isolated
! island or lake depending on the field type.

 data mask /0, 0, 0,  0, 1, 0, 0, 0, 0/

 data latitude /-30.0, -30.0, -30.0, 0., 0., 0., 25.0, 25.0, 25.0/

 data terrain_land /0., 0., 0., 0., 75.0, 0., 0., 0., 0./

 data soilt_climo /0., 0., 0., 0., 2., 0., 0., 0., 0./

! The field values input to the search routine. The
! flag value indicates an unmapped point that must
! be replaced.

 data field/0., 0., 0.,  0., -9999.9, 0.,  0., 0., 0./
 
! The complete list of field numbers the search
! routine works for.

 data default_field_num /0, 1, 7, 11, 21, &
                         30, 65, 66, 83, 85, &
                         86, 91, 92, 223, 224, &
                         225, 226, 227, 228, 229, 230/

! The field value that should be returned by the 
! search routine for each field value. If the returned
! value does not match, this test fails.

 data default_field_val /0.0, 1.0, 75.0, 300.0, 265.0, &
                         30.0, 0.0, 0.0, 0.01, 280.0, &
                         0.18, 0.0, 1.0, 0.0, 2.0, &
                         -99999.9, 0.5, 0.5, 0.5, 1.0, 11.0/

 call mpi_init(ierr)
 
 print*,'RUN TESTS TO CHECK DEFAULT LOGIC'

 do i = 1, num_default_tests

   field_num = default_field_num(i)
   field_updated = field
 
   print*,'CHECK DEFAULT LOGIC FOR FIELD NUMBER ',field_num

   call search (field_updated, mask, idim, jdim, tile, field_num, latitude, terrain_land, &
       soilt_climo)

   if (field_updated(2,2) /= default_field_val(i)) then
     print*,'TEST FAILED ',field_updated(2,2), default_field_val(i)
     stop 2
   endif

 enddo

 call mpi_finalize(ierr)

 print*,'DONE'

 end program test_search_util
