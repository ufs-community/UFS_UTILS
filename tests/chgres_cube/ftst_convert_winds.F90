 program winds

 use esmf

 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        latitude_input_grid, &
                        longitude_input_grid

 use input_data, only : lev_input, convert_winds, &
                        wind_input_grid, &
                        u_input_grid, &
                        v_input_grid

 implicit none

 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)
 real(esmf_kind_r8), allocatable  :: u_wind(:,:,:)
 real(esmf_kind_r8), allocatable  :: v_wind(:,:,:)
 real(esmf_kind_r8), pointer      :: windptr(:,:,:,:)

 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 i_input = 4
 j_input = 3
 lev_input = 1 ! number of vertical levels
 
 print*,'wind check ',i_input, j_input

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)
 print*,'rc ',rc

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=rc)
 print*,'rc ',rc

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)
 print*,'rc ',rc


 print*,"- CALL FieldCreate FOR INPUT wind."
 wind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_input,3/), rc=rc)
 print*,'rc ',rc

 u_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)
 print*,'rc ',rc

 v_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)
 print*,'rc ',rc

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))

 allocate(u_wind(i_input,j_input,lev_input))
 allocate(v_wind(i_input,j_input,lev_input))

! Test point 1 - A pure west wind at the Equator/Greenwich.

 latitude(1,1) = 0.
 longitude(1,1) = 0.
 u_wind(1,1,1) = 1.0
 v_wind(1,1,1) = 0.0

! Test point 2 - A pure west wind at the Equator/90 degrees E longitude.

 latitude(2,1) = 0.
 longitude(2,1) = 90.
 u_wind(2,1,1) = 1.0
 v_wind(2,1,1) = 0.0

! Test point 3 - A pure west wind at the Equator/dateline.

 latitude(3,1) = 0.
 longitude(3,1) = 180.
 u_wind(3,1,1) = 1.0
 v_wind(3,1,1) = 0.0

! Test point 4 - A pure west wind at the Equator/90 degrees W longitude.

 latitude(4,1) = 0.
 longitude(4,1) = 270.
 u_wind(4,1,1) = 1.0
 v_wind(4,1,1) = 0.0

! Test point 5 - A pure south wind at the Equator/Greenwich.

 latitude(1,2) = 0.
 longitude(1,2) = 0.
 u_wind(1,2,1) = 0.0
 v_wind(1,2,1) = 1.0

! Test point 6 - A pure south wind at the Equator/90 degrees E longitude.

 latitude(2,2) = 0.
 longitude(2,2) = 90.
 u_wind(2,2,1) = 0.0
 v_wind(2,2,1) = 1.0

! Test point 7 - A pure south wind at the Equator/dateline.

 latitude(3,2) = 0.
 longitude(3,2) = 180.
 u_wind(3,2,1) = 0.0
 v_wind(3,2,1) = 1.0

! Test point 8 - A pure south wind at the Equator/90 degrees W longitude.

 latitude(4,2) = 0.
 longitude(4,2) = 270.
 u_wind(4,2,1) = 0.0
 v_wind(4,2,1) = 1.0

! Test point 9 - A pure south wind at the North pole/Greenwich.

 latitude(1,3) = 90.
 longitude(1,3) = 0.
 u_wind(1,3,1) = 0.0
 v_wind(1,3,1) = 1.0

! Test point 10 - A pure south wind at the South pole/Greenwich.

 latitude(2,3) = -(90.0)
 longitude(2,3) = 0.
 u_wind(2,3,1) = 0.0
 v_wind(2,3,1) = 1.0

! Test point 11 - A west wind at the 45 degrees N/Greenwich.

 latitude(3,3) = 45.
 longitude(3,3) = 0.
 u_wind(3,3,1) = 1.0
 v_wind(3,3,1) = 0.0

! Test point 12 - A south wind at 45 degrees N/Dateline.

 latitude(4,3) = 45.
 longitude(4,3) = 180.
 u_wind(4,3,1) = 0.0
 v_wind(4,3,1) = 1.0

 call ESMF_FieldScatter(u_input_grid, u_wind, rootpet=0, rc=rc)
 print*,'rc ',rc
 call ESMF_FieldScatter(v_input_grid, v_wind, rootpet=0, rc=rc)
 print*,'rc ',rc
 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 print*,'rc ',rc
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)
 print*,'rc ',rc

 call convert_winds

 print*,"- CALL FieldGet FOR 3-D WIND."
 call ESMF_FieldGet(wind_input_grid, &
                    farrayPtr=windptr, rc=rc)
 print*,'rc ',rc

 do j = 1, j_input
 do i = 1, i_input
 do k = 1, lev_input
   print*,'wind ', i,j,k,windptr(i,j,k,:)
 enddo
 enddo
 enddo

 deallocate(u_wind, v_wind, latitude, longitude)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program winds
