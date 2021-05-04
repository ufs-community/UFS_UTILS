 program winds

 use esmf

 use model_grid, only : i_input, j_input, &
                        input_grid

 implicit none

 integer                      :: ierr, localpet, npets, rc

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)

 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_field) :: latitude_input_grid
 type(esmf_field) :: longitude_input_grid
 type(esmf_field) :: wind_input_grid

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 i_input = 3
 j_input = 3
 
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

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=rc)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))

 latitude = 55.
 longitude = 180.

 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)

 print*,"- CALL FieldCreate FOR INPUT wind."
 wind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/1,3/), rc=rc)

 print*,"SUCCESS!"

 end program winds
