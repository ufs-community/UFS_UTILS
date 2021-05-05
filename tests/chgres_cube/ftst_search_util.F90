 program test_search_util

! Test search_util using a simple 3x3 fv3 tile. Two types
! of tests are performed. First, a missing value is 
! replaced with a valid neighboring value. Second, a missing
! value is replaced by a default value. This can happen
! for an isolated island located far from a valid neighbor.
!
! author: George Gayno (george.gayno@noaa.gov)

 use mpi
 use esmf
 use search_util

 implicit none

 integer, parameter :: IDIM = 3
 integer, parameter :: JDIM = 3
 integer, parameter :: TILE = 1
 integer, parameter :: NUM_DEFAULT_TESTS = 20
 integer, parameter :: NUM_DEFAULT_SST_TESTS = 4

 real, parameter :: EPSILON=0.00001

 integer :: field_num, ierr, i, j, ii, jj
 integer(esmf_kind_i8) :: mask_search1(IDIM,JDIM)
 integer(esmf_kind_i8) :: mask_search2(IDIM,JDIM)
 integer(esmf_kind_i8) :: mask_default(IDIM,JDIM)
 integer :: default_field_num(NUM_DEFAULT_TESTS)

 real(esmf_kind_r8) :: default_field_val(NUM_DEFAULT_TESTS)
 real(esmf_kind_r8) :: default_sst_val(NUM_DEFAULT_SST_TESTS)
 real(esmf_kind_r8) :: field_default(IDIM,JDIM)
 real(esmf_kind_r8) :: field_updated(IDIM,JDIM)
 real(esmf_kind_r8) :: field_search1(IDIM,JDIM)
 real(esmf_kind_r8) :: field_search2(IDIM,JDIM)
 real(esmf_kind_r8) :: latitude(IDIM,JDIM)
 real(esmf_kind_r8) :: latitude_sst(NUM_DEFAULT_SST_TESTS)
 real(esmf_kind_r8) :: terrain_land(IDIM,JDIM)
 real(esmf_kind_r8) :: soilt_climo(IDIM,JDIM)

!--------------------------------------------------------
! These variables are used to test the 'default'
! search routine option (i.e., when the search
! fails, a default must be used)
!--------------------------------------------------------

! Definition of the mask. The '1' indicates an isolated
! island or lake depending on the field type.

 data mask_default /0, 0, 0,  0, 1, 0, 0, 0, 0/

 data latitude /-30.0, -30.0, -30.0, 0., 0., 0., 25.0, 25.0, 25.0/

 data terrain_land /0., 0., 0., 0., 75.0, 0., 0., 0., 0./

 data soilt_climo /0., 0., 0., 0., 2., 0., 0., 0., 0./ ! soil type

! The field values input to the search routine. The
! flag value (-9999.9) indicates an unmapped point that must
! be replaced.

 data field_default/0., 0., 0.,  0., -9999.9, 0.,  0., 0., 0./

! A list of field numbers the default search
! works for. SST is handled with separate data statements
! below.

 data default_field_num /0, 1, 7, 21, &
                         30, 65, 66, 83, 85, &
                         86, 91, 92, 223, 224, &
                         225, 226, 227, 228, 229, 230/

! The field value that should be returned by the 
! search routine for each field value. If the returned
! value does not match this value, the test fails.

 data default_field_val /0.0, 1.0, 75.0, 265.0, &
                         30.0, 0.0, 0.0, 0.01, 280.0, &
                         0.18, 0.0, 1.0, 0.0, 2.0, &
                         -99999.9, 0.5, 0.5, 0.5, 1.0, 11.0/

! For SST, test the default for four latitudes to ensure
! all 'if' branches of routine 'sst_guess' are invoked.
! If the returned value does not match "default_sst_val",
! the test failes.

 data default_sst_val /273.16, 286.5785, 300.0, 273.16/

 data latitude_sst /75.0, 45.0, 0.0, -65.0/

!--------------------------------------------------------
! These variables are used for the two search option 
! tests. Both tests use vegetation greenness. For this
! test, the logic is is independent of the field type.
!--------------------------------------------------------

! Test 1 - The missing value at (2,2) should be replaced
! with the valid value at (1,1).

 data mask_search1 /1, 0, 0,  0, 1, 0, 0, 0, 0/
 data field_search1 /.88, 0., 0.,  0., -9999.9, 0.,  0., 0., 0./

! Test 2 - The missing value at (3,3) should be replaced
! with the valid value at (2,2).

 data mask_search2 /0, 0, 0,  0, 1, 0, 0, 0, 1/
 data field_search2 /0., 0., 0.,  0., .88, 0.,  0., 0., -9999.9/
 
 print*,"Starting test of search util."

 call mpi_init(ierr)
 
 print*,'Run test 1 to check search logic.'

 field_num = 226 ! veg greenness
 field_updated = field_search1

 call search (field_updated, mask_search1, IDIM, JDIM, TILE, field_num)

 do j = 1, JDIM
 do i = 1, IDIM
   if (i==2 .and.j==2) then ! Check the replaced value.
     if (abs(field_updated(i,j) - field_updated(1,1)) > EPSILON) stop 2
   else ! Ensure other values are unchanged.
     if (abs(field_updated(i,j) - field_search1(i,j)) > EPSILON) stop 3
   endif
 enddo
 enddo

 print*,'Run test 2 to check search logic.'

 field_num = 226
 field_updated = field_search2

 call search (field_updated, mask_search2, IDIM, JDIM, TILE, field_num)

 do j = 1, JDIM
 do i = 1, IDIM
   if (i==3 .and.j==3) then ! Check the replaced value.
     if (abs(field_updated(i,j) - field_updated(3,3)) > EPSILON) stop 4
   else ! Ensure other values are unchanged.
     if (abs(field_updated(i,j) - field_search2(i,j)) > EPSILON) stop 5
   endif
 enddo
 enddo

 print*,'Run tests to check default logic.'

 do i = 1, NUM_DEFAULT_TESTS

   field_num = default_field_num(i)
   field_updated = field_default
 
   print*,'CHECK DEFAULT LOGIC FOR FIELD NUMBER ',field_num

   call search (field_updated, mask_default, IDIM, JDIM, TILE, field_num, &
                latitude, terrain_land, soilt_climo)

   do jj = 1, JDIM
   do ii = 1, IDIM
     if (ii==2 .and. jj==2) then ! Check the replaced value.
       if (abs(field_updated(ii,jj)-default_field_val(i)) > EPSILON) stop 6
     else ! Ensure other values are unchanged.
       if (abs(field_updated(ii,jj) - field_default(ii,jj)) > EPSILON) stop 7
     endif
   enddo
   enddo

 enddo

 print*,'Run tests to check default logic for SST.'

 do i = 1, NUM_DEFAULT_SST_TESTS

   field_num = 11
   field_updated = field_default

   latitude(2,2) = latitude_sst(i)

   print*,'CHECK DEFAULT LOGIC FOR FIELD NUMBER ',field_num
   print*,'AT LATITUDE ',latitude_sst(i)

   call search (field_updated, mask_default, IDIM, JDIM, TILE, field_num, &
                latitude, terrain_land, soilt_climo)

   do jj = 1, JDIM
   do ii = 1, IDIM
     if (ii==2 .and. jj==2) then ! Check the replaced value.
       if (abs(field_updated(ii,jj) - default_sst_val(i)) > EPSILON) stop 8
     else ! Ensure other values are unchanged.
       if (abs(field_updated(ii,jj) - field_default(ii,jj)) > EPSILON) stop 9
     endif
   enddo
   enddo

 enddo

 call mpi_finalize(ierr)

 print*,"SUCCESS!"

 end program test_search_util
