 program surface_interp

! Unit test for surface routine interp that regrids surface 
! variables from input to target grid. 
!
! Author: Larissa Reames, OU CIMMS/NOAA NSSL

 use esmf

 use model_grid, only : i_target, j_target, &
                        target_grid, num_tiles_target_grid, &
                        latitude_target_grid, &
                        longitude_target_grid, &
                        lsoil_target

 use program_setup, only : external_model, input_type

 use surface, only : search_many

 use utilities, only : error_handler

 implicit none

 integer, parameter           :: IPTS_TARGET=3
 integer, parameter           :: JPTS_TARGET=3

 real, parameter              :: EPSILON=0.0001
 real(esmf_kind_r8)           :: deltalon

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k, num_fields
 integer                      :: isrctermprocessing

 integer(esmf_kind_i8),allocatable :: mask_target_search(:,:), &
                                      mask_default(:,:)
 integer, allocatable             :: field_nums(:)
 real(esmf_kind_r8), allocatable  :: latitude(:,:), longitude(:,:)
 real(esmf_kind_r8), allocatable  :: field1_search(:,:), & 
                                     field2_search(:,:), &
                                     field1_default(:,:), &
                                     latitude_default(:,:), &
                                     terrain_land(:,:), &
                                     soilt_climo(:,:), &
                                     soil_temp_search(:,:,:)
 real(esmf_kind_r8), allocatable  :: field1_search_correct(:,:), & 
                                     field2_search_correct(:,:), &
                                     field_default_correct(:,:), &
                                     soil_temp_correct(:,:)
 real(esmf_kind_r8), allocatable  :: dummy_2d(:,:), &
                                     dummy_3d(:,:,:)
 real(esmf_kind_r8), pointer      :: lon_ptr(:,:), &
                                     lat_ptr(:,:)

 character(len=50)                :: fname

 type(esmf_vm)                    :: vm
 type(esmf_field)                 :: field1_target_grid, &
                                     field2_target_grid, &
                                     field3_target_grid, &
                                     field4_target_grid, &
                                     field_3d_target_grid, &
                                     temp_field
 type(esmf_fieldbundle)           :: bundle_search1, &
                                     bundle_search2, &
                                     bundle_default1, &
                                     bundle_default2, &
                                     bundle_3d_search

 print*,"Starting test of surface regrid_many."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 !--------------------------------------------------------------------!
 !---------------- Setup Target Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 i_target = IPTS_TARGET
 j_target = JPTS_TARGET
 lsoil_target = 2

 num_tiles_target_grid = 1
 target_grid = ESMF_GridCreate1PeriDim(maxIndex=(/i_target,j_target/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))  

 ! Regional grid
 deltalon = 0.5
 do i = 1, i_target
   longitude(i,:) = 91.1_esmf_kind_r8 + real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 do i = 1, j_target
   latitude(:,i) = 34.1_esmf_kind_r8 - real((i-1),kind=esmf_kind_r8) * deltalon
 enddo            

  call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_ptr(i,j) = longitude(i,j)
     if (lon_ptr(i,j) > 360.0_esmf_kind_r8) lon_ptr(i,j) = lon_ptr(i,j) -360.0_esmf_kind_r8
     lat_ptr(i,j) = latitude(i,j)
   enddo
 enddo
 nullify(lat_ptr,lon_ptr)


 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=rc)

 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)

 call ESMF_FieldScatter(longitude_target_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_target_grid, latitude, rootpet=0, rc=rc)
 deallocate(latitude, longitude)

 ! Create target fields
 field1_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="field1_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 field2_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="field2_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 field3_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="soil_type_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 field4_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="field4_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 field_3d_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="field_3d_target_grid", &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_target/), rc=rc)

 ! Allocate space for arrays 
 allocate(field1_search_correct(i_target,j_target))
 allocate(field2_search_correct(i_target,j_target))
 allocate(field1_search(i_target,j_target))
 allocate(field2_search(i_target,j_target))
 allocate(mask_target_search(i_target,j_target))

 allocate(field1_default(i_target,j_target))
 allocate(field_default_correct(i_target,j_target))
 allocate(mask_default(i_target,j_target))
 allocate(latitude_default(i_target,j_target))
 allocate(dummy_2d(i_target,j_target))

 allocate(terrain_land(i_target,j_target))
 allocate(soilt_climo(i_target,j_target))

 allocate(soil_temp_search(i_target,j_target,lsoil_target))
 allocate(soil_temp_correct(i_target,j_target))
 allocate(dummy_3d(i_target,j_target,lsoil_target))
 
 ! Field values for default replacement tests
 field1_default = reshape((/0., 0., 0.,  0., -9999.9, 0.,  0., 0.,0./),(/i_target,j_target/))
 mask_default = reshape((/0, 0, 0,  0, 1, 0, 0, 0, 0/),(/i_target,j_target/))
 latitude_default = reshape((/-30.0, -30.0, -30.0, 0., 75., 0., 25.0, 25.0,25.0/),(/i_target,j_target/))
 

 ! Field values to check basic search option tests
 field1_search=reshape((/-9999.9, 0., 0.,  0., .88, 0.,  0., 0.,.1/),(/i_target,j_target/))
 field1_search_correct=reshape((/.88, 0., 0.,  0., .88, 0.,  0.,0.,.1/),(/i_target,j_target/))   
 field2_search=reshape((/3., 0., 0.,  0., 2., 0.,  0., 0., -9999.9/),(/i_target,j_target/))
 field2_search_correct=reshape((/3., 0., 0.,  0., 2., 0.,  0., 0.,2./),(/i_target,j_target/))
 mask_target_search=reshape((/1, 0, 0,  0, 1, 0, 0, 0, 1/),(/i_target,j_target/))
 soil_temp_search(:,:,1) = reshape((/-9999.9, 0., 0.,  0., 280., 0.,  0.,0.,290./),(/i_target,j_target/))
 soil_temp_search(:,:,2) = reshape((/-9999.9, 0., 0.,  0., 280., 0.,0.,0.,290./),(/i_target,j_target/))
 soil_temp_correct(:,:) = reshape((/280., 0., 0.,  0., 280.,0.,0.,0.,290./),(/i_target,j_target/))
 ! Default terrain values to check default terrain replacement
 terrain_land = reshape((/0., 0., 0., 0., 75.0, 0., 0., 0.,0./),(/i_target,j_target/))

 ! Climatology soil type values to check soil type replacement
 soilt_climo = reshape((/0., 0., 0., 0., 2., 0., 0., 0.,0./),(/i_target,j_target/))

 ! Create field bundles and assign fields to them
 bundle_default1 = ESMF_FieldBundleCreate(name="fields_default1", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)

 ! will search sst, terrain height, soil_type_target_grid
 call ESMF_FieldBundleAdd(bundle_default1, (/field1_target_grid,field2_target_grid, & 
         field3_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)                        

 bundle_default2 = ESMF_FieldBundleCreate(name="fields_default2", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 ! will search GFS grib2 soil type
 call ESMF_FieldBundleAdd(bundle_default2,(/field1_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)

  bundle_3d_search = ESMF_FieldBundleCreate(name="fields_search_3d", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 ! will search soil temperature
 call ESMF_FieldBundleAdd(bundle_3d_search,(/field_3d_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)

 bundle_search1 = ESMF_FieldBundleCreate(name="fields_search1", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 ! will search veg greeness and restart soil type
 call ESMF_FieldBundleAdd(bundle_search1,(/field1_target_grid,field2_target_grid/),rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)

 bundle_search2 = ESMF_FieldBundleCreate(name="fields_search2", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 ! will search hrrr grib2 non-target-grid soil type
 call ESMF_FieldBundleAdd(bundle_search2,(/field1_target_grid/),rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)

!-------------------------------------------------------------------------------------
! SEARCH TEST CHECKS REPLACEMENT OF VEG FRACTION AND RESTART FILE SOIL TYPE
!-------------------------------------------------------------------------------------

 ! Fill esmf fields for search test
  call ESMF_FieldScatter(field1_target_grid, field1_search, rootPet=0,tile=1, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

 call ESMF_FieldScatter(field2_target_grid, field2_search, rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc) 
                           
 call ESMF_FieldBundleGet(bundle_search1,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)  
 
 allocate(field_nums(num_fields))
 field_nums = (/226,224/)
 input_type="restart"

 !Call the search many routine to test search and replace
 call search_many(num_fields,bundle_search1,1,field_nums,localpet, &
                  soilt_climo=soilt_climo,mask=mask_target_search)

 call ESMF_FieldBundleDestroy(bundle_search1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)  
 deallocate(field_nums)

 call ESMF_FieldGather(field1_target_grid, dummy_2d, rootPet=0, rc=rc)

  print*,"Check results for field1_search."

 if (any((abs(dummy_2d - field1_search_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'field1_search SHOULD BE:', field1_search_correct
   print*,'field1_search FROM TEST:', dummy_2d
   stop 2
 endif
 call ESMF_FieldGather(field2_target_grid, dummy_2d, rootPet=0, rc=rc)

 print*,"Check results for field2_search."
 if (any((abs(dummy_2d - field2_search_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'field2_search SHOULD BE:', field2_search_correct
   print*,'field2_search FROM TEST:', dummy_2d
   stop 2
 endif

!-------------------------------------------------------------------------------------
! SEARCH TEST CHECKS REPLACEMENT OF HRRR GRIB2 SOIL NO TYPE TARGET GRID
!-------------------------------------------------------------------------------------

 ! Fill esmf fields for search test
  call ESMF_FieldScatter(field1_target_grid, field2_search, rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

 call ESMF_FieldBundleGet(bundle_search2,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleGet", rc)

 allocate(field_nums(num_fields))
 field_nums = (/224/)
 input_type="grib2"
 external_model="HRRR"

 !Call the search many routine to test search and replace
  call search_many(num_fields,bundle_search2,1,field_nums,localpet, &
         soilt_climo=soilt_climo,mask=mask_target_search)

 call ESMF_FieldBundleDestroy(bundle_search2,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)
 deallocate(field_nums)

 call ESMF_FieldGather(field1_target_grid, dummy_2d, rootPet=0, rc=rc)

  print*,"Check results for field2_search."

 if (any((abs(dummy_2d - field2_search_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'field2_search SHOULD BE:', field2_search_correct
   print*,'field2_search FROM TEST:', dummy_2d
   stop 2
 endif

!-------------------------------------------------------------------------------------
! DEFAULT TEST 1 CHECKS DEFAULT/CLIMO SST,TERRAIN,SOILTYPE REPLACEMENT
!-------------------------------------------------------------------------------------

 ! Fill esmf fields for default1 test
  call ESMF_FieldScatter(field1_target_grid, field1_default, rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

 call ESMF_FieldScatter(field2_target_grid, field1_default,rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

 call ESMF_FieldScatter(field3_target_grid, field1_default,rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

 call ESMF_FieldBundleGet(bundle_default1,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleGet", rc)

 allocate(field_nums(num_fields))
 field_nums = (/11,7,224/)
 !Call the search many routine to test some branches of default behavior
 call search_many(num_fields,bundle_default1,1,field_nums,localpet, &
      latitude=latitude_default,terrain_land=terrain_land,soilt_climo=soilt_climo,mask=mask_default)

 print*,"Check results for bundle_default1."

 do i = 1,num_fields
   call ESMF_FieldBundleGet(bundle_default1,i,temp_field,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldBundleGet", rc)
   call ESMF_FieldGet(temp_field, name=fname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGet", rc)
   print*, "Check ", trim(fname)
   call ESMF_FieldGather(temp_field,dummy_2d,rootPet=0,tile=1,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)

   field_default_correct = field1_default
   if (i .eq. 1) then
     field_default_correct(2,2) = 273.16
   elseif (i .eq. 2) then
     field_default_correct(2,2) = terrain_land(2,2)
   else
     field_default_correct(2,2) = soilt_climo(2,2)
   endif

   if (any((abs(dummy_2d - field_default_correct)) > EPSILON)) then
     print*,'TEST FAILED '
     print*,trim(fname), ' SHOULD BE:', field_default_correct
     print*,trim(fname), ' FROM TEST:', dummy_2d
     stop 2
   endif
 enddo
  call ESMF_FieldBundleDestroy(bundle_default1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleDestroy", rc)
 deallocate(field_nums)

!---------------------------------------------
! DEFAULT TEST 2 TESTS GFS GRIB2 SOIL TYPE
!---------------------------------------------
  ! Fill esmf fields for default2 test
  call ESMF_FieldScatter(field1_target_grid, field1_default,rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

  call ESMF_FieldBundleGet(bundle_default2,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleGet", rc)

 allocate(field_nums(num_fields))
 field_nums(:) = (/224/)

 input_type="grib2"
 external_model="GFS"
 !Call the search many routine to test behavior for GFS grib2 soil type
 call search_many(num_fields,bundle_default2,1,field_nums,localpet,&
      soilt_climo=soilt_climo,mask=mask_default)

 call ESMF_FieldBundleDestroy(bundle_default2,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleDestroy", rc)
 deallocate(field_nums)

 print*,"Check results for bundle_default2."

 call ESMF_FieldGather(field1_target_grid,dummy_2d,rootPet=0,tile=1,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)

 if (any((abs(dummy_2d - soilt_climo)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'field1_target SHOULD BE:', soilt_climo
   print*,'field1_target FROM TEST:', dummy_2d
   stop 2
 endif
 
!--------------------------------------------------------
! 3D TEST  TESTS REPLACEMENT FOR SOIL TEMPERATURE
!-------------------------------------------------------!  
! Fill esmf fields for default2 test
  call ESMF_FieldScatter(field_3d_target_grid,soil_temp_search,rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

  call ESMF_FieldBundleGet(bundle_3d_search,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleGet", rc)

 allocate(field_nums(num_fields))
 field_nums(:) = (/21/)

 !Call the search many routine to test behavior for GFS grib2 soil type
 call search_many(num_fields,bundle_3d_search,1,field_nums,localpet,mask=mask_target_search)

 call ESMF_FieldBundleDestroy(bundle_3d_search,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleDestroy", rc)
 deallocate(field_nums)

 print*,"Check results for bundle_3d_search."

 call ESMF_FieldGather(field_3d_target_grid,dummy_3d,rootPet=0,tile=1,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)
 do i = 1,lsoil_target
   if (any((abs(dummy_3d(:,:,i) - soil_temp_correct)) > EPSILON)) then
     print*,'TEST FAILED '
     print*,'field_3d_target SHOULD BE:', soil_temp_correct
     print*,'field_3d_target at level ',i,' FROM TEST:', dummy_3d(:,:,i)
     stop 2
   endif
 enddo

 print*,"Tests Passed!"

! Deallocate and destroy
 deallocate(field1_search_correct,field2_search_correct,field1_search,field2_search,mask_target_search)
 deallocate(field1_default,mask_default,latitude_default,dummy_2d,terrain_land,soilt_climo,dummy_3d)
 deallocate(soil_temp_correct,soil_temp_search,field_default_correct)

 call ESMF_FieldDestroy(latitude_target_grid,rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid,rc=rc)
 call ESMF_FieldDestroy(field1_target_grid,rc=rc)
 call ESMF_FieldDestroy(field2_target_grid,rc=rc)
 call ESMF_FieldDestroy(field3_target_grid,rc=rc)
 call ESMF_FieldDestroy(field4_target_grid,rc=rc)
 call ESMF_FieldDestroy(field_3d_target_grid,rc=rc)
 call ESMF_GridDestroy(target_grid, rc=rc)


 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program surface_interp
