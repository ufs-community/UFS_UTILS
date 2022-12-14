 program surface_interp

! Unit test for surface routine interp that regrids surface 
! variables from input to target grid. 
!
! Author: Larissa Reames, OU CIMMS/NOAA NSSL

 use esmf


 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        latitude_input_grid, &
                        longitude_input_grid, &
                        i_target, j_target, &
                        target_grid, num_tiles_target_grid, &
                        latitude_target_grid, &
                        longitude_target_grid

 use sfc_input_data, only: t2m_input_grid, &
                       q2m_input_grid
 
 use surface, only : regrid_many

 use surface_target_data, only : t2m_target_grid, &
                                 q2m_target_grid

 implicit none

 integer, parameter           :: IPTS_INPUT=4
 integer, parameter           :: JPTS_INPUT=3
 integer, parameter           :: IPTS_TARGET=8
 integer, parameter           :: JPTS_TARGET=5

 real, parameter              :: EPSILON=0.0001
 real(esmf_kind_r8)           :: deltalon

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k, num_fields
 integer                      :: isrctermprocessing

 real(esmf_kind_r8), allocatable  :: latitude(:,:), longitude(:,:)
 real(esmf_kind_r8), allocatable  :: q2m_input(:,:), & 
                                     t2m_input(:,:)
 real(esmf_kind_r8), allocatable  :: q2m_correct(:,:), & 
                                     q2m_target(:,:), &
                                     t2m_target(:,:), &
                                     t2m_correct(:,:)
 real(esmf_kind_r8), pointer      :: lon_ptr(:,:), &
                                     lat_ptr(:,:)
 type(esmf_vm)                    :: vm
 type(esmf_polekind_flag)         :: polekindflag(2)
 type(esmf_fieldbundle)           :: bundle_all_target, bundle_all_input
 type(esmf_regridmethod_flag)     :: method
 type(esmf_routehandle)           :: regrid_bl_no_mask
 logical, allocatable             :: dozero(:)

 print*,"Starting test of surface regrid_many."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 !--------------------------------------------------------------------!
 !----------------- Setup Input Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 i_input = IPTS_INPUT
 j_input = JPTS_INPUT

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 input_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_input,j_input/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))  
          
 ! This is a random regional grid. I tried a global grid here but it had an unstable
 ! solution.
 
 deltalon = 2.0_esmf_kind_r8
 do i = 1, i_input
   longitude(i,:) = 90+real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 do j = 1, j_input
   latitude(:,j) = 35.0-real((j-1),kind=esmf_kind_r8) * deltalon
 end do

 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 call ESMF_GridGetCoord(input_grid, &
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
     if (lon_ptr(i,j) > 360.0_esmf_kind_r8) lon_ptr(i,j) = lon_ptr(i,j) - 360.0_esmf_kind_r8
     lat_ptr(i,j) = latitude(i,j)
   enddo
 enddo
 nullify(lat_ptr,lon_ptr)
 

 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)
 deallocate(latitude, longitude) 

 !Initializes input ESMF fields
 t2m_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 q2m_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 !Allocate and fill in the fields on the input grid that we need to create soil type
 allocate(t2m_input(i_input,j_input))
 allocate(q2m_input(i_input,j_input))
 
 t2m_input = reshape((/290.,292.,294.,296., 291.,293.,295.,297., 292.,294.,296.,298./),(/i_input,j_input/))
 q2m_input = reshape((/6.E-4,7.E-4,8.E-4,9.E-4,  7.E-4,8.E-4,9.E-4,10.E-4,   8.E-4,9.E-4,10.E-4,11.E-4/),(/i_input,j_input/))
 
 call ESMF_FieldScatter(t2m_input_grid,t2m_input,rootpet=0,rc=rc)
 call ESMF_FieldScatter(q2m_input_grid,q2m_input,rootpet=0,rc=rc)
 
 deallocate(t2m_input,q2m_input)

 !--------------------------------------------------------------------!
 !---------------- Setup Target Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 i_target = IPTS_TARGET
 j_target = JPTS_TARGET

 num_tiles_target_grid = 1
 target_grid = ESMF_GridCreate1PeriDim(maxIndex=(/i_target,j_target/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))  

 ! Regional grid that fits within the input regional grid but with smaller grid cells
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

 ! Create target t2m and q2m fields
 t2m_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="t2m_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 q2m_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="q2m_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)

 ! Create masks on the target grid and the correct (expected) soil type on the target grid
 ! to check against what returns from interp
 
 allocate(t2m_correct(i_target,j_target))
 allocate(q2m_correct(i_target,j_target))
 allocate(t2m_target(i_target,j_target))
 allocate(q2m_target(i_target,j_target))

   
 !t2m_correct =   reshape((/0., 0., 15.,15.,5., 5., 5., 5., &  
 !                          0., 0., 5., 5., 6., 6., 6., 6., &
 !                          0., 0., 5., 5., 6., 6., 6., 6., &
 !                          0., 0., 5., 5., 6., 6., 0., 0., &
 !                          0., 0., 5., 5., 6., 6., 0., 0. /),(/i_target,j_target/))
 t2m_correct = reshape((/ 292.000000000000,        292.000000000000,&
   292.000000000000,        292.000000000000,        294.000000000000,&
   294.000000000000,        294.000000000000,        294.000000000000,&
   293.000000000000,        293.000000000000,        293.000000000000,&
   293.000000000000,        295.000000000000,        295.000000000000,&
   295.000000000000,        295.000000000000,        293.000000000000,&
   293.000000000000,        293.000000000000,        293.000000000000,&
   295.000000000000,        295.000000000000,        295.000000000000,&
   295.000000000000,        293.000000000000,        293.000000000000,&
   293.000000000000,        293.000000000000,        295.000000000000,&
   295.000000000000,        295.000000000000,        295.000000000000,&
   293.000000000000,        293.000000000000,        293.000000000000,&
   293.000000000000,        295.000000000000,        295.000000000000,&
   295.000000000000,        295.000000000000/),(/i_target,j_target/))
 !q2m_correct =   reshape((/0., 0.,16.,16., 4., 4., 4., 4., &
 !                          0., 0., 3., 3., 5., 5., 5., 5., &
 !                          0., 0., 3., 3., 5., 5., 5., 5., &
 !                          0., 0., 3., 3., 5., 5., 0., 0., &
 !                          0., 0., 3., 3., 5., 5., 0., 0. /),(/i_target,j_target/))
 q2m_correct = reshape((/ 7.000000000000000E-004,  7.000000000000000E-004,&
  7.000000000000000E-004,  7.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  8.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  8.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  9.000000000000000E-004,  9.000000000000000E-004,&
  9.000000000000000E-004,  9.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  8.000000000000000E-004,  8.000000000000000E-004,&
  9.000000000000000E-004,  9.000000000000000E-004,  9.000000000000000E-004,&
  9.000000000000000E-004,  8.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  8.000000000000000E-004,  9.000000000000000E-004,&
  9.000000000000000E-004,  9.000000000000000E-004,  9.000000000000000E-004,&
  8.000000000000000E-004,  8.000000000000000E-004,  8.000000000000000E-004,&
  8.000000000000000E-004,  9.000000000000000E-004,  9.000000000000000E-004,&
  9.000000000000000E-004,  9.000000000000000E-004/),(/i_target,j_target/))


 method=ESMF_REGRIDMETHOD_NEAREST_STOD

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR NON-MASKED BILINEAR INTERPOLATION."
 call ESMF_FieldRegridStore(t2m_input_grid, &
                            t2m_target_grid, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl_no_mask, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 bundle_all_target = ESMF_FieldBundleCreate(name="all points target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_all_input = ESMF_FieldBundleCreate(name="all points input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
      
 call ESMF_FieldBundleAdd(bundle_all_target, (/t2m_target_grid,q2m_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)                        
 call ESMF_FieldBundleAdd(bundle_all_input, (/t2m_input_grid,q2m_input_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)  
                            
 call ESMF_FieldBundleGet(bundle_all_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)  
 
 allocate(dozero(num_fields))
 dozero(:) = .True. 
  
 !Call the routine to unit test.
 call regrid_many(bundle_all_input,bundle_all_target,num_fields,regrid_bl_no_mask,dozero)
 deallocate(dozero) 

 call ESMF_FieldBundleDestroy(bundle_all_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)  
 call ESMF_FieldBundleDestroy(bundle_all_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc) 

 call ESMF_FieldGather(t2m_target_grid, t2m_target, rootPet=0, rc=rc)
 call ESMF_FieldGather(q2m_target_grid, q2m_target, rootPet=0, rc=rc)

 print*,"Check results."

 if (any((abs(t2m_target - t2m_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'T2M SHOULD BE:', t2m_correct
   print*,'T2M FROM TEST:', t2m_target
   stop 2
 endif
 
  if (any((abs(q2m_target - q2m_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'Q2M SHOULD BE:', q2m_correct
   print*,'Q2M FROM TEST:', q2m_target
   stop 2
 endif

 
 print*,"OK"

! Deallocate and destroy
 deallocate(t2m_target,t2m_correct,q2m_target,q2m_correct)
 call ESMF_FieldDestroy(latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(longitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(latitude_target_grid,rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid,rc=rc)
 call ESMF_FieldDestroy(t2m_input_grid,rc=rc)
 call ESMF_FieldDestroy(t2m_input_grid,rc=rc)
 call ESMF_FieldDestroy(q2m_input_grid,rc=rc)
 call ESMF_FieldDestroy(q2m_input_grid,rc=rc)
call ESMF_GridDestroy(input_grid, rc=rc)
 call ESMF_GridDestroy(target_grid, rc=rc)


 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program surface_interp
