 program surface_interp

! Unit test for surface routine interp that regrids surface 
! variables from input to target grid. 
!
! Author: Larissa Reames, OU CIMMS/NOAA NSSL

 use esmf
 
 use program_setup, only : sotyp_from_climo, &
                           vgtyp_from_climo

 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        latitude_input_grid, &
                        longitude_input_grid, &
                        i_target, j_target, &
                        target_grid, num_tiles_target_grid, &
                        latitude_target_grid, &
                        longitude_target_grid, &
                        landmask_target_grid, &
                        seamask_target_grid, &
                        terrain_target_grid, &
                        cleanup_input_target_grid_data

 use static_data, only: create_static_fields, &
                        soil_type_target_grid, &
                        veg_type_target_grid, &
                        cleanup_static_fields
 
 use surface, only : create_surface_esmf_fields, &
                     interp, &
                     cleanup_target_sfc_data
 
 use input_data, only : init_sfc_esmf_fields, &
                        soil_type_input_grid, &
                        veg_type_input_grid, &
                        landsea_mask_input_grid, &
                        terrain_input_grid, &
                        t2m_input_grid, &
                        cleanup_input_sfc_data

 implicit none

 integer, parameter           :: IPTS_INPUT=4
 integer, parameter           :: JPTS_INPUT=3
 integer, parameter           :: IPTS_TARGET=8
 integer, parameter           :: JPTS_TARGET=5

 real, parameter              :: EPSILON=0.0001
 real(esmf_kind_r8)           :: deltalon

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k

 integer(esmf_kind_i8), pointer   :: mask_target_ptr(:,:)
 real(esmf_kind_r8), allocatable  :: latitude(:,:), longitude(:,:)
 real(esmf_kind_r8), allocatable  :: mask_input(:,:)
 integer(esmf_kind_i8), allocatable  :: mask_target(:,:), &
                                     seamask_target(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_input(:,:), & 
                                     vgtyp_input(:,:), &
                                     terrain_input(:,:), &
                                     t2m_input(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_correct(:,:), & 
                                     sotyp_target(:,:), &
                                     vgtyp_target(:,:), &
                                     vgtyp_correct(:,:), &
                                     terrain_correct(:,:)
 real(esmf_kind_r8), pointer      :: lon_ptr(:,:), &
                                     lat_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_corner_ptr(:,:), &
                                     lat_corner_ptr(:,:)
 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)

 print*,"Starting test of surface interp."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 sotyp_from_climo = .False.
 vgtyp_from_climo = .False.

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
 
! Create staggered coordinates for conservative regridding
 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

  print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (i == i_input+1) then
       lon_corner_ptr(i,j) = longitude(i_input,1) + (0.5_esmf_kind_r8*deltalon)
     else
       lon_corner_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
     endif

     if (j == j_input+1) then
       lat_corner_ptr(i,j) = latitude(1,j_input) +(0.5_esmf_kind_r8*deltalon)
     else
       lat_corner_ptr(i,j) = latitude(1,j) -(0.5_esmf_kind_r8*deltalon)
     endif
   enddo
 enddo

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

 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)
 deallocate(latitude, longitude) 
 nullify(lat_corner_ptr, lon_corner_ptr)

 !Initializes surface ESMF fields
 call init_sfc_esmf_fields
 
 !Allocate and fill in the fields on the input grid that we need to create soil type
 allocate(mask_input(i_input,j_input))
 allocate(sotyp_input(i_input,j_input))
 allocate(vgtyp_input(i_input,j_input))
 allocate(terrain_input(i_input,j_input))
 
 mask_input = reshape((/0,1,1,1, 0,1,1,1, 0,1,1,0/),(/i_input,j_input/))
 sotyp_input = reshape((/0.,16.,4.,16.,  0.,3.,5.,16.,   0.,3.,5.,0./),(/i_input,j_input/))
 vgtyp_input = reshape((/0.,15.,5.,15.,  0.,5.,6.,15.,  0.,5.,6.,0./),(/i_input,j_input/))
 terrain_input = reshape((/0.,20.,20.,20.,  0.,20.,20.,20., 0.,20., 20.,0./),(/i_input,j_input/))
 
 call ESMF_FieldScatter(landsea_mask_input_grid,mask_input,rootpet=0,rc=rc)
 call ESMF_FieldScatter(soil_type_input_grid,sotyp_input,rootpet=0,rc=rc)
 call ESMF_FieldScatter(veg_type_input_grid,vgtyp_input,rootpet=0,rc=rc)
 call ESMF_FieldScatter(terrain_input_grid,terrain_input,rootpet=0,rc=rc)
 
 deallocate(mask_input,sotyp_input,vgtyp_input,terrain_input)

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

 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", rc)

 ! Create staggered coordinates for regional target grid
 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (i == i_target+1) then
       lon_corner_ptr(i,j) = longitude(i_input,1) + (0.5_esmf_kind_r8*deltalon)
     else
       lon_corner_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
     endif

     if (j == j_target+1) then
       lat_corner_ptr(i,j) = latitude(1,j_input) +(0.5_esmf_kind_r8*deltalon)
     else
       lat_corner_ptr(i,j) = latitude(1,j) -(0.5_esmf_kind_r8*deltalon)
     endif
   enddo
 enddo

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
 nullify(lat_corner_ptr, lon_corner_ptr)

 ! Create the fields for the target grid land and seamask since these would normally
 ! be created in the appropriate model_grid  subroutine
 
  print*,"- CALL FieldCreate FOR TARGET GRID LANDMASK."
 landmask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_landmask", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate",rc)

 print*,"- CALL FieldCreate FOR TARGET GRID SEAMASK."
 seamask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_seamask", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)
 
 ! Initialize other surface fields on the target grid
 call create_surface_esmf_fields()
 call create_static_fields()

 ! Create masks on the target grid and the correcte (expected) soil type on the target grid
 ! to check against what returns from interp
 
 allocate(mask_target(i_target,j_target))
 allocate(seamask_target(i_target,j_target))
 allocate(sotyp_target(i_target,j_target))
 allocate(sotyp_correct(i_target,j_target))

 mask_target(:,1) = (/0,0,1,1,1,1,1,1/) 
 mask_target(:,2) = (/0,0,1,1,1,1,1,1/)
 mask_target(:,3) = (/0,0,1,1,1,1,1,1/)
 mask_target(:,4) = (/0,0,1,1,1,1,0,0/)
 mask_target(:,5) = (/0,0,1,1,1,1,0,0/)
   
!vgtyp_correct = reshape((/0., 0., 15.,15.,5., 5., 5., 5., &  
!                          0., 0., 5., 5., 6., 6., 6., 6., &
!                          0., 0., 5., 5., 6., 6., 6., 6., &
!                          0., 0., 5., 5., 6., 6., 0., 0., &
!                          0., 0., 5., 5., 6., 6., 0., 0. /),(/i_target,j_target/))
 sotyp_correct = reshape((/0., 0.,16.,16., 4., 4., 4., 4., &
                           0., 0., 3., 3., 5., 5., 5., 5., &
                           0., 0., 3., 3., 5., 5., 5., 5., &
                           0., 0., 3., 3., 5., 5., 0., 0., &
                           0., 0., 3., 3., 5., 5., 0., 0. /),(/i_target,j_target/))
 seamask_target = 0
 where(mask_target .eq. 0) seamask_target = 1

 call ESMF_FieldScatter(landmask_target_grid,mask_target,rootpet=0,rc=rc)
 call ESMF_FieldScatter(seamask_target_grid,seamask_target,rootpet=0,rc=rc)

 !Call the routine to unit test.
 call interp(localpet)

 call ESMF_FieldGather(soil_type_target_grid, sotyp_target, rootPet=0, rc=rc)

 print*,"Check results."

 if (any((abs(sotyp_target - sotyp_correct)) > EPSILON)) then
   print*,'TEST FAILED '
   print*,'ARRAY SHOULD BE:', sotyp_correct
   print*,'ARRAY FROM TEST:', sotyp_target
   stop 2
 endif

 
 print*,"OK"

 deallocate(mask_target, sotyp_target, sotyp_correct,seamask_target)

 call cleanup_static_fields
 call cleanup_input_sfc_data
 call cleanup_target_sfc_data
 call cleanup_input_target_grid_data
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program surface_interp
