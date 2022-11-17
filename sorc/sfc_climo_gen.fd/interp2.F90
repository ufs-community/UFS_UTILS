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
 use model_grid,  only               : grid_mdl, i_mdl, j_mdl, &
                                       num_tiles, latitude_field_mdl, &
                                       longitude_field_mdl, mask_field_mdl, &
                                       land_frac_field_mdl
 use source_grid
 use utils
 use mpi

 implicit none

 character(len=*), intent(in)       :: input_file

 integer                            :: rc, localpet
 integer                            :: i, j, tile, ncid, status
 integer                            :: varid, water_category
 integer                            :: isrctermprocessing
 integer :: category, num_categories

 integer(esmf_kind_i4), allocatable :: mask_mdl_one_tile(:,:)
 integer(esmf_kind_i4), pointer     :: unmapped_ptr(:)

 real(esmf_kind_r4), allocatable    :: data_src_global(:,:)
 real(esmf_kind_r4), allocatable    :: data_src_global2(:,:,:)
 real(esmf_kind_r4), allocatable    :: data_mdl_one_tile(:,:,:)
 real(esmf_kind_r4), allocatable    :: dom_cat_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lat_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: sum_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lon_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: land_frac_mdl_one_tile(:,:)

 type(esmf_regridmethod_flag) :: method
 type(esmf_field)                        :: data_field_src
 type(esmf_field)                        :: data_field_mdl
 type(esmf_routehandle)                  :: regrid_data
 type(esmf_polemethod_flag)              :: pole
 

 if (localpet == 0) then
   allocate(data_src_global(i_src,j_src))
 else
   allocate(data_src_global(0,0))
 endif

 if (localpet == 0) then
   print*,'- OPEN SOURCE FILE ', trim(input_file)
   status = nf90_open(input_file, nf90_nowrite, ncid)
   call netcdf_err(status, "IN ROUTINE INTERP OPENING SOURCE FILE")
   status = nf90_inq_varid(ncid, field_names(1), varid)
   call netcdf_err(status, "IN ROUTINE INTERP READING FIELD ID")
   status = nf90_get_var(ncid, varid, data_src_global, start=(/1,1,1/), count=(/i_src,j_src,1/))
   call netcdf_err(status, "IN ROUTINE INTERP READING FIELD")
   print*,'number of cats ',maxval(data_src_global)
   num_categories = nint(maxval(data_src_global))
   status = nf90_get_att(ncid, varid, 'water_category', water_category)
   call netcdf_err(status, "IN ROUTINE INTERP READING water_category")
   print*,'water cat ',water_category
 endif

 call mpi_bcast(num_categories,1,MPI_INTEGER,0,MPI_COMM_WORLD,status)

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
 data_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                  typekind=ESMF_TYPEKIND_R4, &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
                                  ungriddedLBound=(/1/), &
                                  ungriddedUBound=(/num_categories/), &
                                  name="mdl data", &
                                  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (localpet == 0) then
   allocate(data_src_global2(i_src,j_src,num_categories))
   allocate(data_mdl_one_tile(i_mdl,j_mdl,num_categories))
   allocate(dom_cat_mdl_one_tile(i_mdl,j_mdl))
   allocate(mask_mdl_one_tile(i_mdl,j_mdl))
   allocate(land_frac_mdl_one_tile(i_mdl,j_mdl))
   allocate(lat_mdl_one_tile(i_mdl,j_mdl))
   allocate(sum_mdl_one_tile(i_mdl,j_mdl))
   allocate(lon_mdl_one_tile(i_mdl,j_mdl))
 else
   allocate(data_src_global2(0,0,0))
   allocate(data_mdl_one_tile(0,0,0))
   allocate(dom_cat_mdl_one_tile(0,0))
   allocate(mask_mdl_one_tile(0,0))
   allocate(land_frac_mdl_one_tile(0,0))
   allocate(lat_mdl_one_tile(0,0))
   allocate(sum_mdl_one_tile(0,0))
   allocate(lon_mdl_one_tile(0,0))
 endif

 if (localpet == 0) then
   data_src_global2 = 0.0
   do j = 1, j_src
   do i = 1, i_src
     category = nint(data_src_global(i,j))
     if (category < 1) cycle
     data_src_global2(i,j,category) = 1.0
   enddo
   enddo
 endif

 deallocate(data_src_global)

 print*,"- CALL FieldScatter FOR SOURCE GRID DATA."
 call ESMF_FieldScatter(data_field_src, data_src_global2, rootPet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter.", rc)

 deallocate(data_src_global2)

 isrctermprocessing = 1

 method = ESMF_REGRIDMETHOD_CONSERVE
 pole = ESMF_POLEMETHOD_NONE

 print*,"- CALL FieldRegridStore."
 nullify(unmapped_ptr)
 call ESMF_FieldRegridStore(data_field_src, &
                            data_field_mdl, &
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

 print*,"- CALL Field_Regrid."
 call ESMF_FieldRegrid(data_field_src, &
                       data_field_mdl, &
                       routehandle=regrid_data, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid.", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_data, rc=rc)

 print*,"- CALL FieldDestroy."
 call ESMF_FieldDestroy(data_field_src, rc=rc)

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
   call ESMF_FieldGather(data_field_mdl, data_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather.", rc)

   print*,"- CALL FieldGather FOR MODEL GRID LAND MASK."
   call ESMF_FieldGather(mask_field_mdl, mask_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather.", rc)

   print*,"- CALL FieldGather FOR MODEL GRID LAND FRACTION."
   call ESMF_FieldGather(land_frac_field_mdl, land_frac_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather.", rc)

   if (localpet == 0) then
     print*,'- CALL SEARCH FOR TILE ',tile

! Where sum is zero, the regridding did not find any input data for the model point
! (ex. and isolated island). Call the search routine at these points.
     sum_mdl_one_tile = sum(data_mdl_one_tile, dim=3)
     do j = 1, j_mdl
     do i = 1, i_mdl
       if (mask_mdl_one_tile(i,j) == 1 .and. sum_mdl_one_tile(i,j) == 0.0) then
         data_mdl_one_tile(i,j,:) = -9999.9 ! flag to tell search routine to search.
       endif
     enddo
     enddo
     call search2 (data_mdl_one_tile, mask_mdl_one_tile, i_mdl, j_mdl, num_categories, tile, field_names(1))
     print*,'after regrid ',data_mdl_one_tile(i_mdl/2,j_mdl/2,:)

! These points are all non-land. Set to 100% of the water category.
     do j = 1, j_mdl
     do i = 1, i_mdl
       if (mask_mdl_one_tile(i,j) == 0) then
         data_mdl_one_tile(i,j,water_category) = 1.0
       endif
     enddo
     enddo

! For fractional grids, need to rescale the category percentages by the
! fraction of land in the model grid cell.

! When running with fractional grids, the land_frac_mdl_one_tile array will 
! contain a fraction between 0 and 1. When not running with fractional
! grids, this array will contain negative fill values.

     if (maxval(land_frac_mdl_one_tile) > 0.0) then
       print*,'before rescale ',land_frac_mdl_one_tile(658,95),data_mdl_one_tile(658,95,:)
       do j = 1, j_mdl
       do i = 1, i_mdl
         if (mask_mdl_one_tile(i,j) == 1) then
           data_mdl_one_tile(i,j,:) = data_mdl_one_tile(i,j,:) * land_frac_mdl_one_tile(i,j)
           data_mdl_one_tile(i,j,water_category) = 1.0 - land_frac_mdl_one_tile(i,j)
         endif
       enddo
       enddo
       print*,'after  rescale ',land_frac_mdl_one_tile(658,95),data_mdl_one_tile(658,95,:)
     endif
! under fractional grids, how do we define dominate category?
     dom_cat_mdl_one_tile = 0.0
     dom_cat_mdl_one_tile = maxloc(data_mdl_one_tile,dim=3)
     call output2 (data_mdl_one_tile, dom_cat_mdl_one_tile, lat_mdl_one_tile, lon_mdl_one_tile, i_mdl, j_mdl, num_categories, tile)
   endif

 enddo OUTPUT_LOOP

 status=nf90_close(ncid)

 deallocate(data_mdl_one_tile, dom_cat_mdl_one_tile, mask_mdl_one_tile)
 deallocate(lat_mdl_one_tile, lon_mdl_one_tile, sum_mdl_one_tile, land_frac_mdl_one_tile)

 print*,"- CALL FieldDestroy."
 call ESMF_FieldDestroy(data_field_mdl, rc=rc)

 call mpi_barrier(mpi_comm_world, rc)

 end subroutine interp2
