 program winds

! Unit test for routine convert_winds, which
! converts u/v component winds to x/y/z component
! winds on a sphere. This test converts unit
! vectors at several latitude and longitudes.
!
! Author: George Gayno NCEP/EMC

 use esmf

 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        latitude_input_grid, &
                        longitude_input_grid

 use atm_input_data, only : lev_input, convert_winds, &
                        wind_input_grid, &
                        u_input_grid, &
                        v_input_grid

 implicit none

 integer, parameter           :: IPTS=4
 integer, parameter           :: JPTS=3

 real, parameter              :: EPSILON=0.0001

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)
 real(esmf_kind_r8), allocatable  :: u_wind(:,:,:)
 real(esmf_kind_r8), allocatable  :: v_wind(:,:,:)
 real(esmf_kind_r8), pointer      :: windptr(:,:,:,:)

 real :: expected_x_component(IPTS,JPTS)
 real :: expected_y_component(IPTS,JPTS)
 real :: expected_z_component(IPTS,JPTS)

 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)

! These are the expected x/y/z components return from
! convert_winds for our test points.

 data expected_x_component/1.0, 0.0, -1.0, 0.0, &
                           0.0, 0.0,  0.0, 0.0, &
                           0.0, 0.0,  1.0, 0.0 / 

 data expected_y_component/0.0, 1.0,  0.0, -1.0, &
                           0.0, 0.0,  0.0, 0.0, &
                           1.0, -1.0, 0.0, -0.707106/ 

 data expected_z_component/0.0, 0.0, 0.0, 0.0, &
                           1.0, 1.0, 1.0, 1.0, &
                           0.0, 0.0, 0.0, 0.707106/ 

 print*,"Starting test of convert_winds."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 i_input = IPTS
 j_input = JPTS
 lev_input = 1 ! number of vertical levels

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=rc)

 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)

 wind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_input,3/), rc=rc)

 u_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)

 v_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))

 allocate(u_wind(i_input,j_input,lev_input))
 allocate(v_wind(i_input,j_input,lev_input))

! Test point 1 - A west wind at the Equator/Greenwich.

 latitude(1,1) = 0.
 longitude(1,1) = 0.
 u_wind(1,1,1) = 1.0
 v_wind(1,1,1) = 0.0

! Test point 2 - A west wind at the Equator/90 degrees East longitude.

 latitude(2,1) = 0.
 longitude(2,1) = 90.
 u_wind(2,1,1) = 1.0
 v_wind(2,1,1) = 0.0

! Test point 3 - A west wind at the Equator/Dateline.

 latitude(3,1) = 0.
 longitude(3,1) = 180.
 u_wind(3,1,1) = 1.0
 v_wind(3,1,1) = 0.0

! Test point 4 - A west wind at the Equator/90 degrees West longitude.

 latitude(4,1) = 0.
 longitude(4,1) = 270.
 u_wind(4,1,1) = 1.0
 v_wind(4,1,1) = 0.0

! Test point 5 - A south wind at the Equator/Greenwich.

 latitude(1,2) = 0.
 longitude(1,2) = 0.
 u_wind(1,2,1) = 0.0
 v_wind(1,2,1) = 1.0

! Test point 6 - A south wind at the Equator/90 degrees East longitude.

 latitude(2,2) = 0.
 longitude(2,2) = 90.
 u_wind(2,2,1) = 0.0
 v_wind(2,2,1) = 1.0

! Test point 7 - A south wind at the Equator/Dateline.

 latitude(3,2) = 0.
 longitude(3,2) = 180.
 u_wind(3,2,1) = 0.0
 v_wind(3,2,1) = 1.0

! Test point 8 - A south wind at the Equator/90 degrees West longitude.

 latitude(4,2) = 0.
 longitude(4,2) = 270.
 u_wind(4,2,1) = 0.0
 v_wind(4,2,1) = 1.0

! Test point 9 - A south wind at the North Pole/Greenwich.

 latitude(1,3) = 90.
 longitude(1,3) = 0.
 u_wind(1,3,1) = 0.0
 v_wind(1,3,1) = 1.0

! Test point 10 - A south wind at the South Pole/Greenwich.

 latitude(2,3) = -(90.0)
 longitude(2,3) = 0.
 u_wind(2,3,1) = 0.0
 v_wind(2,3,1) = 1.0

! Test point 11 - A west wind at the 45 degrees North latitude/Greenwich.

 latitude(3,3) = 45.
 longitude(3,3) = 0.
 u_wind(3,3,1) = 1.0
 v_wind(3,3,1) = 0.0

! Test point 12 - A south wind at 45 degrees North latitude/Dateline.

 latitude(4,3) = 45.
 longitude(4,3) = 180.
 u_wind(4,3,1) = 0.0
 v_wind(4,3,1) = 1.0

 call ESMF_FieldScatter(u_input_grid, u_wind, rootpet=0, rc=rc)
 call ESMF_FieldScatter(v_input_grid, v_wind, rootpet=0, rc=rc)
 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)

! Call the routine to unit test.

 call convert_winds

 call ESMF_FieldGet(wind_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=windptr, rc=rc)

 print*,"Check results."

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
 do k = clb(3), cub(3)
   if (abs(windptr(i,j,k,1) - expected_x_component(i,j)) > EPSILON) stop 2
   if (abs(windptr(i,j,k,2) - expected_y_component(i,j)) > EPSILON) stop 3
   if (abs(windptr(i,j,k,3) - expected_z_component(i,j)) > EPSILON) stop 4
 enddo
 enddo
 enddo

 print*,"OK"

 deallocate(u_wind, v_wind, latitude, longitude)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program winds
