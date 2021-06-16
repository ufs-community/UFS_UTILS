 program surface_nst_landfill 

! Unit test for surface routine interp that regrids surface 
! variables from input to target grid. 
!
! Author: Larissa Reames, OU CIMMS/NOAA NSSL

 use esmf
 
 use model_grid, only : i_target, j_target, &
                        target_grid, num_tiles_target_grid, &
                        landmask_target_grid

 use surface, only : skin_temp_target_grid, &
                     nst_land_fill, &
                     create_nst_esmf_fields, &
                     cleanup_target_nst_data, &
                     c_d_target_grid, &
                     c_0_target_grid, &
                     d_conv_target_grid, &
                     dt_cool_target_grid, &
                     ifd_target_grid, &
                     qrain_target_grid, &
                     tref_target_grid, &
                     w_d_target_grid, &
                     w_0_target_grid, &
                     xs_target_grid, &
                     xt_target_grid, &
                     xu_target_grid, &
                     xv_target_grid, &
                     xz_target_grid, &
                     xtts_target_grid, &
                     xzts_target_grid, &
                     z_c_target_grid, &
                     zm_target_grid
                     
 
 implicit none

 integer, parameter           :: IPTS_TARGET=4
 integer, parameter           :: JPTS_TARGET=4

 real, parameter              :: EPSILON=0.0001
 real(esmf_kind_r8)           :: deltalon, &
                                 init_val = -999.9

 integer                      :: ierr, localpet, npets, rc
 integer                      :: i

 character(len=50)            :: fname

 integer(esmf_kind_i8), pointer     :: mask_target_ptr(:,:)
 integer(esmf_kind_i8), allocatable :: mask(:,:)
 
 real(esmf_kind_r8), pointer        :: skin_temp_ptr(:,:), &
                                       tmp_ptr(:,:)
                    
 real(esmf_kind_r8), allocatable    :: skin_temp(:,:), &
                                       tmp_array_2d(:,:), &
                                       nst_field_correct(:,:), &
                                       tref_correct(:,:),&
                                       xz_correct(:,:)
                                    
 type(esmf_vm)                :: vm
 type(esmf_field)             :: tmp_field
 type(esmf_fieldbundle)       :: nst_bundle

 print*,"Starting test of nst_land_fill."

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)

 !--------------------------------------------------------------------!
 !---------------- Setup Target Grid & Coordinates -------------------!
 !--------------------------------------------------------------------!
 
 i_target = IPTS_TARGET
 j_target = JPTS_TARGET

 num_tiles_target_grid = 1
 target_grid = ESMF_GridCreate1PeriDim(maxIndex=(/i_target,j_target/), &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

  call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN GridAddCoord", rc)

 landmask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_landmask", &
                                   rc=rc)

 skin_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_skin_temp", &
                                   rc=rc)

 call create_nst_esmf_fields

 allocate(mask(i_target,j_target))
 allocate(skin_temp(i_target,j_target))
 allocate(nst_field_correct(i_target,j_target))
 allocate(tref_correct(i_target,j_target))
 allocate(xz_correct(i_target,j_target))

 mask = reshape((/0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0/),(/i_target,j_target/))
 skin_temp = reshape((/280., 280., 280., 280., &
                       280., 290., 290., 280., &
                       280., 290., 290., 280., &
                       280., 280., 280., 280./),(/i_target,j_target/))
 call ESMF_FieldScatter(landmask_target_grid, mask, rootpet=0, rc=rc)
 call ESMF_FieldScatter(skin_temp_target_grid, skin_temp, rootpet=0, rc=rc)
 deallocate(mask, skin_temp)

 nst_field_correct = reshape((/-999.9, -999.9, -999.9, -999.9, &
                               -999.9,     0.,     0., -999.9, &
                               -999.9,     0.,     0., -999.9, &
                               -999.9, -999.9, -999.9, -999.9/),(/i_target,j_target/))

 tref_correct =      reshape((/-999.9, -999.9, -999.9, -999.9, &
                               -999.9,  290.0,  290.0, -999.9, &
                               -999.9,  290.0,  290.0, -999.9, &
                               -999.9, -999.9, -999.9, -999.9/),(/i_target,j_target/))

 xz_correct =        reshape((/-999.9, -999.9, -999.9, -999.9, &
                               -999.9,   30.0,   30.0, -999.9, &
                               -999.9,   30.0,   30.0, -999.9, &
                               -999.9, -999.9, -999.9, -999.9/),(/i_target,j_target/))
 
 nst_bundle = ESMF_FieldBundleCreate(name="nst_bundle", fieldlist= &
                        (/c_d_target_grid,c_0_target_grid,d_conv_target_grid, &
                          dt_cool_target_grid,ifd_target_grid,qrain_target_grid,&
                          w_d_target_grid,w_0_target_grid,xs_target_grid,xt_target_grid,&
                          xu_target_grid,xv_target_grid,xtts_target_grid,xzts_target_grid,&
                          z_c_target_grid, zm_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 !call ESMF_FieldBundleAdd(nst_bundle, (/c_d_target_grid,c_0_target_grid,d_conv_target_grid, &
 !                         dt_cool_target_grid,ifd_target_grid,qrain_target_grid,&
 !                         w_d_target_grid,w_0_target_grid,xs_target_grid,xt_target_grid,&
 !                         xu_target_grid,xv_target_grid,xtts_target_grid,xzts_target_grid, &
 !                         z_c_target_grid, zm_target_grid/), rc=rc)
 !  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
 !     call error_handler("IN FieldBundleAdd", rc) 

 do i = 1,16
   call ESMF_FieldBundleGet(nst_bundle,i,tmp_field,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldBundleGet", rc)
   call ESMF_FieldGet(tmp_field,farrayPtr=tmp_ptr,rc=rc)
   tmp_ptr(:,:) = init_val
 enddo
 call ESMF_FieldGet(tref_target_grid,farrayPtr=tmp_ptr,rc=rc)
   tmp_ptr(:,:) = init_val
 call ESMF_FieldGet(xz_target_grid,farrayPtr=tmp_ptr,rc=rc)
   tmp_ptr(:,:) = init_val

 !Call the routine to unit test.
 call nst_land_fill

 allocate(tmp_array_2d(i_target,j_target))

 print*,"Check results."
 do i = 1, 16
   call ESMF_FieldBundleGet(nst_bundle,i,tmp_field,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldBundleGet", rc)
   call ESMF_FieldGather(tmp_field,tmp_array_2d,rootPet=0,tile=1,rc=rc)
   call ESMF_FieldGet(tmp_field,name=fname,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGet", rc)

   if (any((abs(tmp_array_2d - nst_field_correct)) > EPSILON)) then
     print*,'TEST FAILED '
     print*,'ARRAY ', trim(fname), ' SHOULD BE:', nst_field_correct
     print*,'ARRAY ', trim(fname), ' FROM TEST:', tmp_array_2d
     stop 2
   endif
 enddo
 call ESMF_FieldBundleDestroy(nst_bundle,rc=rc)

 call ESMF_FieldGather(tref_target_grid,tmp_array_2d,rootPet=0,tile=1,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGater", rc)

 if (any((abs(tmp_array_2d - tref_correct)) > EPSILON)) then
     print*,'TEST FAILED '
     print*,'TREF SHOULD BE:', tref_correct
     print*,'TREF FROM TEST:', tmp_array_2d
     stop 2
 endif

  call ESMF_FieldGather(xz_target_grid,tmp_array_2d,rootPet=0,tile=1,rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGater", rc)

 if (any((abs(tmp_array_2d - xz_correct)) > EPSILON)) then
     print*,'TEST FAILED '
     print*,'XZ SHOULD BE:', xz_correct
     print*,'XZ FROM TEST:', tmp_array_2d
     stop 2
 endif

 
 print*,"OK"

 deallocate(tmp_array_2d,nst_field_correct,tref_correct,xz_correct)

 call ESMF_FieldDestroy(landmask_target_grid,rc=rc)
 call ESMF_FieldDestroy(skin_temp_target_grid,rc=rc)
 call cleanup_target_nst_data
 call ESMF_GridDestroy(target_grid,rc=rc)
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)
 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program surface_nst_landfill
