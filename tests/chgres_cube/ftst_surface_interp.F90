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
                        target_grid, &
                        latitude_target_grid, &
                        longitude_target_grid, &
                        landmask_target_grid, &
                        terrain_target_grid

 use static_data, only: create_static_fields, &
                        soil_type_target_grid, &
                        veg_type_target_grid
 
 use surface, only : create_surface_esmf_fields, &
                     interp
 
 use input_data, only : init_sfc_esmf_fields, &
                        soil_type_input_grid, &
                        veg_type_input_grid, &
                        landsea_mask_input_grid, &
                        terrain_input_grid

 implicit none

 integer, parameter           :: IPTS_INPUT=4
 integer, parameter           :: JPTS_INPUT=3
 integer, parameter           :: IPTS_TARGET=8
 integer, parameter           :: JPTS_TARGET=5

 real, parameter              :: EPSILON=0.0001

 integer                      :: clb(4), cub(4)
 integer                      :: ierr, localpet, npets, rc
 integer                      :: i, j, k

 real(esmf_kind_r8), allocatable  :: lati_input(:,:), lon_input(:,:)
 real(esmf_kind_i8), allocatable  :: mask_input(:,:), mask_target(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_input(:,:), & 
                                     vgtyp_input(:,:), &
                                     terrain_input(:,:)
 real(esmf_kind_r8), allocatable  :: sotyp_correct(:,:), & 
                                     sotyp_target(:,:), &
                                     vgtyp_correct(:,:), &
                                     terrain_correct(:,:)

 type(esmf_vm)                :: vm
 type(esmf_polekind_flag)     :: polekindflag(2)

 print*,"Starting test of convert_winds."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 i_input = IPTS_INPUT
 j_input = JPTS_INPUT

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
 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))  
                                
 deltalon = 360.0_esmf_kind_r8 / real(i_input,kind=esmf_kind_r8)
 do i = 1, i_input
   longitude(i,:) = real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 allocate(slat(j_input))
 allocate(wlat(j_input))
 call splat(4, j_input, slat, wlat)

 do i = 1, j_input
   latitude(:,i) = 90.0_esmf_kind_r8 - (acos(slat(i))* 180.0_esmf_kind_r8 / &
                  (4.0_esmf_kind_r8*atan(1.0_esmf_kind_r8)))
 enddo   
 
 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)     
 
 deallocate(latitude,longitude)                           
 
 call init_esmf_sfc_fields
 
 allocate(mask(i_input,j_input))
 allocate(sotyp(i_input,j_input))
 allocate(vgtyp(i_input,j_input))
 allocate(terrain(i_input,j_input))
 
 mask_input = reshape((/0,0,0, 1,1,1, 1,1,1, 1,1,0/),(/i_input,j_input/))
 sotyp_input = reshape((/0.,0.,0., 16.,3.,3., 4.,5.,5., 16.,16.,0./),(/i_input,j_input/))
 vgtyp_input = reshape((/0.,0.,0., 15.,5.,5., 5.,6., 6., 15.,15.,0./),(/i_input,j_input/))
 terrain = reshape((/0.,0.,0., 20.,20.,20., 20.,20.,20., 20.,20.,0/),(/i_input,j_input/))
 
 call ESMF_FieldScatter(landsea_mask_input_grid,mask,rootpet=0,rc=rc)
 call ESMF_FieldScatter(sotyp_input_grid,sotyp,rootpet=0,rc=rc)
 call ESMF_FieldScatter(vgtyp_input_grid,vgtyp,rootpet=0,rc=rc)
 call ESMF_FieldScatter(terrain_input_grid,terrain,rootpet=0,rc=rc)
 
 deallocate(mask,sotyp,vgtyp,terrain)

 i_target = IPTS_TARGET
 j_target = JPTS_TARGET

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 target_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_target,j_target/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

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
 
 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))  
                                
 deltalon = 360.0_esmf_kind_r8 / real(i_input,kind=esmf_kind_r8)
 do i = 1, i_input
   longitude(i,:) = real((i-1),kind=esmf_kind_r8) * deltalon
 enddo

 allocate(slat(j_input))
 allocate(wlat(j_input))
 call splat(4, j_input, slat, wlat)

 do i = 1, j_input
   latitude(:,i) = 90.0_esmf_kind_r8 - (acos(slat(i))* 180.0_esmf_kind_r8 / &
                  (4.0_esmf_kind_r8*atan(1.0_esmf_kind_r8)))
 enddo            
 
 call create_surface_esmf_fields
 call create_static_fields
 deallocate(mask_target, sotyp_target, sotyp_correct)

 allocate(mask_target(i_target,j_target))
 allocate(sotyp_target(i_target,j_target))
 allocate(sotyp_correct(i_target,j_target))
 
 mask_target = reshape((/0,0,0,0,0, 0,0,0,0,0, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1,  &
                           1,1,1,1,1, 1,1,1,0,0, 1,1,1,0,0/),(/i_target,j_target/))
 !vgtyp_correct = reshape((/0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 15., 15., 5., 5., 5., &
 !                        15., 15., 6., 5., 5., 15., 15., 15., 6., 6., 15., 15., 15., 6., 6., &
 !                        15., 15., 15., 0., 0., 15., 15., 15., 0., 0./),(/i_target,j_target/))
 sotyp_correct = reshape((/0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 16., 16., 3., 3., 3., &
                         16., 16., 5., 3., 3., 16., 16., 16., 5., 5., 16., 16., 16., 5., 5., &
                         16., 16., 16., 0., 0., 16., 16., 16., 0., 0./),(/i_target,j_target/))
 
 call ESMF_FieldScatter(landmask_target_grid,mask,rootpet=0,rc=rc)
 call ESMF_FieldScatter(sotyp_target_grid,sotyp,rootpet=0,rc=rc)


! Call the routine to unit test.
 call interp(localpet)


 call ESMF_FieldGather(sotyp_target_grid, sotyp_target, rootPet=0, rc=rc)

 if (localpet==0) then
 
   print*,"Check results."

   if any((abs(sotyp_target - sotyp_correct)) > EPSILON) then
     print*,'TEST FAILED '
     print*,'ARRAY SHOULD BE:', soityp_correct
     print*,'ARRAY FROM TEST:', soilt_updated
     stop 2
   endif
   
 endif
 
 print*,"OK"

 deallocate(mask_target, sotyp_target, sotyp_correct)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program winds
