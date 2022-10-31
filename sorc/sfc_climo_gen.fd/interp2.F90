!> @file
!! @brief Read the input source data and interpolate it to the
!! model grid.
!! @author George Gayno @date 2018

!> Read the input source data and interpolate it to the
!! model grid.
!!
!! @param[in] localpet this mpi task
!! @param[in] input_file filename of input source data.
!! @author George Gayno @date 2018
 subroutine interp2(localpet, input_file)

 use esmf
 use netcdf
 use model_grid
 use source_grid
 use utils
 use mpi

 implicit none

 character(len=*), intent(in)       :: input_file

 integer                            :: rc, localpet
 integer                            :: i, j, tile, n, ncid, status
 integer                            :: t
 integer                            :: clb_mdl(3), cub_mdl(3)
 integer                            :: varid, record
 integer                            :: isrctermprocessing
 integer :: category, num_categories

 integer(esmf_kind_i4), allocatable :: mask_mdl_one_tile(:,:)
 integer(esmf_kind_i4), pointer     :: unmapped_ptr(:)

 real(esmf_kind_r4), pointer        :: data_mdl_ptr(:,:,:)
 real(esmf_kind_r4), allocatable    :: data_src_global(:,:)
 real(esmf_kind_r4), allocatable    :: data_src_global2(:,:,:)
 real(esmf_kind_r4), allocatable    :: data_mdl_one_tile(:,:,:)
 real(esmf_kind_r4), allocatable    :: vegt_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lat_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: sum_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lon_mdl_one_tile(:,:)

 type(esmf_regridmethod_flag) :: method
 type(esmf_field)                        :: data_field_src
 type(esmf_field)                        :: data_field_mdl2
 type(esmf_routehandle)                  :: regrid_data
 type(esmf_polemethod_flag)              :: pole
 
! get this from file.
 num_categories = 20

 print*,"- CALL FieldCreate FOR SOURCE GRID DATA."
 data_field_src = ESMF_FieldCreate(grid_src, &
                                  typekind=ESMF_TYPEKIND_R4, &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
                                  ungriddedLBound=(/1/), &
                                  ungriddedUBound=(/num_categories/), &
                                  name="source data", &
                                  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR model GRID veg DATA."
 data_field_mdl2 = ESMF_FieldCreate(grid_mdl, &
                                  typekind=ESMF_TYPEKIND_R4, &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
                                  ungriddedLBound=(/1/), &
                                  ungriddedUBound=(/num_categories/), &
                                  name="mdl data", &
                                  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldGet FOR MODEL GRID DATA."
 nullify(data_mdl_ptr)
 call ESMF_FieldGet(data_field_mdl2, &
                    farrayPtr=data_mdl_ptr,  &
                    computationalLBound=clb_mdl, &
                    computationalUBound=cub_mdl, &
                    rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)


 print*,'got here2 ',localpet,clb_mdl,cub_mdl

 if (localpet == 0) then
   allocate(data_src_global(i_src,j_src))
   allocate(data_src_global2(i_src,j_src,num_categories))
 else
   allocate(data_src_global(0,0))
   allocate(data_src_global2(0,0,0))
 endif

 print*,'- OPEN SOURCE FILE ', trim(input_file)
 status = nf90_open(input_file, nf90_nowrite, ncid)
 call netcdf_err(status, "IN ROUTINE INTERP OPENING SOURCE FILE")

 if (localpet == 0) then
   allocate(data_mdl_one_tile(i_mdl,j_mdl,num_categories))
   allocate(mask_mdl_one_tile(i_mdl,j_mdl))
   allocate(lat_mdl_one_tile(i_mdl,j_mdl))
   allocate(sum_mdl_one_tile(i_mdl,j_mdl))
   allocate(lon_mdl_one_tile(i_mdl,j_mdl))
 else
   allocate(data_mdl_one_tile(0,0,0))
   allocate(mask_mdl_one_tile(0,0))
   allocate(lat_mdl_one_tile(0,0))
   allocate(sum_mdl_one_tile(0,0))
   allocate(lon_mdl_one_tile(0,0))
 endif

 record = 0

 TIME_LOOP : do t = 1, num_time_recs ! loop over each time period

 FIELD_LOOP : do n = 1, num_fields ! loop over each surface field.

   record = record + 1

   if (localpet == 0) then
     status = nf90_inq_varid(ncid, field_names(n), varid)
     call netcdf_err(status, "IN ROUTINE INTERP READING FIELD ID")
     status = nf90_get_var(ncid, varid, data_src_global, start=(/1,1,t/), count=(/i_src,j_src,1/))
     call netcdf_err(status, "IN ROUTINE INTERP READING FIELD")
     data_src_global2 = 0.0
     do j = 1, j_src
     do i = 1, i_src
       category = nint(data_src_global(i,j))
!      if (category < 1) category = 17
       if (category < 1) cycle
       data_src_global2(i,j,category) = 1.0
     enddo
     enddo
   endif

   print*,"- CALL FieldScatter FOR SOURCE GRID DATA."
   call ESMF_FieldScatter(data_field_src, data_src_global2, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter.", rc)

   isrctermprocessing = 1

   if (record == 1) then

     method = ESMF_REGRIDMETHOD_CONSERVE
     pole = ESMF_POLEMETHOD_NONE

     print*,"- CALL FieldRegridStore."
     nullify(unmapped_ptr)
     call ESMF_FieldRegridStore(data_field_src, &
                            data_field_mdl2, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=pole, &                     
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_data, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr,  &
                            rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldRegridStore.", rc)

   endif

   print*,"- CALL Field_Regrid."
   call ESMF_FieldRegrid(data_field_src, &
                         data_field_mdl2, &
                         routehandle=regrid_data, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, &
                         rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid.", rc)
 
!-----------------------------------------------------------------------
! Unmapped points are stored in "unmapped_ptr".  The pointer contains
! "ij" global indices as follows:
!
! tile 1: 1 thru (itile*jtile)
! tile n: (n-1)*(itile*jtile) thru n*(itile*jtile)
! 
! This "ij" index is converted to the tile number and i/j index of the
! field object.  This logic assumes the model grid object was
! created using "GLOBAL" indices. 
!
! Unmapped data points are given the flag value of -9999.9, which
! will be replaced in routine "search".
!-----------------------------------------------------------------------

   OUTPUT_LOOP : do tile = 1, num_tiles

     print*,"- CALL FieldGather FOR MODEL LATITUDE."
     call ESMF_FieldGather(latitude_field_mdl, lat_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather.", rc)

     print*,"- CALL FieldGather FOR MODEL LONGITUDE."
     call ESMF_FieldGather(longitude_field_mdl, lon_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather.", rc)

     print*,"- CALL FieldGather FOR MODEL GRID DATA."
     call ESMF_FieldGather(data_field_mdl2, data_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather.", rc)

     print*,"- CALL FieldGather FOR MODEL GRID LAND MASK."
     call ESMF_FieldGather(mask_field_mdl, mask_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather.", rc)

     select case (trim(field_names(n)))
       case ('substrate_temperature','vegetation_greenness','leaf_area_index','slope_type','soil_type')
         print*,"- CALL FieldGather FOR MODEL GRID VEG TYPE."
         call ESMF_FieldGather(vegt_field_mdl, vegt_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather.", rc)
     end select

     if (localpet == 0) then
       print*,'- CALL SEARCH FOR TILE ',tile
       sum_mdl_one_tile = sum(data_mdl_one_tile, dim=3) ! use unused variable to now.
       do j = 1, j_mdl
       do i = 1, i_mdl
 
         if (mask_mdl_one_tile(i,j) == 1 .and. sum_mdl_one_tile(i,j) == 0.0) then
           data_mdl_one_tile(i,j,:) = -9999.9
         endif

       enddo
       enddo

       call search2 (data_mdl_one_tile, mask_mdl_one_tile, i_mdl, j_mdl, num_categories, tile, field_names(n))
!      where(mask_mdl_one_tile == 0) data_mdl_one_tile = missing
       print*,'after regrid ',data_mdl_one_tile(i_mdl/2,j_mdl/2,:)
       call output2 (data_mdl_one_tile, lat_mdl_one_tile, lon_mdl_one_tile, i_mdl, j_mdl, num_categories, tile, t, n)
     endif

   print*,'after output ', localpet
   call mpi_barrier(mpi_comm_world, rc)
   stop

!    if (field_names(n) == 'vegetation_type') then
!      print*,"- CALL FieldScatter FOR MODEL GRID VEGETATION TYPE."
!      call ESMF_FieldScatter(vegt_field_mdl, data_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
!      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!         call error_handler("IN FieldScatter.", rc)
!    endif

   enddo OUTPUT_LOOP

   if (allocated(vegt_mdl_one_tile)) deallocate(vegt_mdl_one_tile)

 enddo FIELD_LOOP
 enddo TIME_LOOP

 status=nf90_close(ncid)

 deallocate(data_mdl_one_tile, mask_mdl_one_tile)
 deallocate(data_src_global, lat_mdl_one_tile, lon_mdl_one_tile)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_data, rc=rc)

 print*,"- CALL FieldDestroy."
 call ESMF_FieldDestroy(data_field_src, rc=rc)

 end subroutine interp2
