!> @file
!! @brief Read the input source data and interpolate it to the
!! model grid.
!! @author George Gayno @date 2018

!> Read the input source data and interpolate it to the
!! model grid.
!!
!! @param[in] localpet this mpi task
!! @param[in] method interpolation method.defined where mask=1
!! @param[in] input_file filename of input source data.
!! @author George Gayno @date 2018
 subroutine interp(localpet, method, input_file)

 use esmf
 use netcdf
 use model_grid
 use program_setup, only : fract_vegsoil_type
 use source_grid
 use utils
 use mpi

 implicit none

 character(len=*), intent(in)       :: input_file

 integer                            :: rc, localpet
 integer                            :: i, j, ij, tile, n, ncid, status
 integer                            :: l(1), u(1), t
 integer                            :: clb_mdl(2), cub_mdl(2)
 integer                            :: varid, record
 integer                            :: tile_num, pt_loc_this_tile
 integer                            :: isrctermprocessing
 double precision                   :: scale
 integer(esmf_kind_i4), allocatable :: mask_mdl_one_tile(:,:)
 integer(esmf_kind_i4), pointer     :: unmapped_ptr(:)

 real(esmf_kind_r4), pointer        :: data_mdl_ptr(:,:)
 real(esmf_kind_r4), pointer        :: data_src_ptr(:,:)
 real(esmf_kind_r4), allocatable    :: data_src_global(:,:)
 real(esmf_kind_r4), allocatable    :: data_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: vegt_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lat_mdl_one_tile(:,:)
 real(esmf_kind_r4), allocatable    :: lon_mdl_one_tile(:,:)

 type(esmf_regridmethod_flag),intent(in) :: method
 type(esmf_field)                        :: data_field_src
 type(esmf_routehandle)                  :: regrid_data
 type(esmf_polemethod_flag)              :: pole
 
 print*,"- CALL FieldCreate FOR SOURCE GRID DATA."
 data_field_src = ESMF_FieldCreate(grid_src, &
                                  typekind=ESMF_TYPEKIND_R4, &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
                                  name="source data", &
                                  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldGet FOR SOURCE GRID DATA."
 nullify(data_src_ptr)
 call ESMF_FieldGet(data_field_src, &
                    farrayPtr=data_src_ptr,  &
                    rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR MODEL GRID DATA."
 nullify(data_mdl_ptr)
 call ESMF_FieldGet(data_field_mdl, &
                    farrayPtr=data_mdl_ptr,  &
                    computationalLBound=clb_mdl, &
                    computationalUBound=cub_mdl, &
                    rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (localpet == 0) then
   allocate(data_src_global(i_src,j_src))
 else
   allocate(data_src_global(0,0))
 endif

 print*,'- OPEN SOURCE FILE ', trim(input_file)
 status = nf90_open(input_file, nf90_nowrite, ncid)
 call netcdf_err(status, "IN ROUTINE INTERP OPENING SOURCE FILE")

 if (localpet == 0) then
   allocate(data_mdl_one_tile(i_mdl,j_mdl))
   allocate(mask_mdl_one_tile(i_mdl,j_mdl))
   allocate(lat_mdl_one_tile(i_mdl,j_mdl))
   allocate(lon_mdl_one_tile(i_mdl,j_mdl))
 else
   allocate(data_mdl_one_tile(0,0))
   allocate(mask_mdl_one_tile(0,0))
   allocate(lat_mdl_one_tile(0,0))
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
     status=nf90_get_att(ncid, varid, 'scale_factor', scale)
      if (status == 0) then
          call scale_data(data_src_global,i_src,j_src,scale)
      endif

   endif

   print*,"- CALL FieldScatter FOR SOURCE GRID DATA."
   call ESMF_FieldScatter(data_field_src, data_src_global, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter.", rc)

   isrctermprocessing = 1

   if (record == 1) then

     if (method == ESMF_REGRIDMETHOD_BILINEAR) then
       pole = ESMF_POLEMETHOD_ALLAVG
     else
       pole = ESMF_POLEMETHOD_NONE
     endif

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

   endif

   print*,"- CALL Field_Regrid."
   call ESMF_FieldRegrid(data_field_src, &
                         data_field_mdl, &
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

   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)
   do ij = l(1), u(1)

     tile_num = ((unmapped_ptr(ij)-1) / (i_mdl*j_mdl)) ! tile number minus 1
     pt_loc_this_tile = unmapped_ptr(ij) - (tile_num * i_mdl * j_mdl) 
                                                  ! "ij" location of point within tile.

     j = (pt_loc_this_tile - 1) / i_mdl + 1
     i = mod(pt_loc_this_tile, i_mdl)
     if (i==0) i = i_mdl
     data_mdl_ptr(i,j) = -9999.9 
                                 
   enddo

! Adjust some fields at permanent land ice points. These points are identified by the
! 'permanent ice' vegetation type category.
!
! When outputting the fraction of each vegetation type, land ice points are
! not defined. So don't do this adjustment.

   if (.not. fract_vegsoil_type) then
     select case (trim(field_names(n)))
       case ('substrate_temperature','vegetation_greenness','leaf_area_index','slope_type','soil_type','soil_color')
       if (localpet == 0) then
         allocate(vegt_mdl_one_tile(i_mdl,j_mdl))
       else
         allocate(vegt_mdl_one_tile(0,0))
       endif
     end select
   endif

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

     if (.not. fract_vegsoil_type) then
       select case (trim(field_names(n)))
         case ('substrate_temperature','vegetation_greenness','leaf_area_index','slope_type','soil_type','soil_color')
           print*,"- CALL FieldGather FOR MODEL GRID VEG TYPE."
           call ESMF_FieldGather(vegt_field_mdl, vegt_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
           if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather.", rc)
       end select
     endif

     if (localpet == 0) then
       print*,'- CALL SEARCH FOR TILE ',tile
       call search (data_mdl_one_tile, mask_mdl_one_tile, i_mdl, j_mdl, tile, field_names(n))
       if (.not. fract_vegsoil_type) then
         select case (field_names(n))
           case ('substrate_temperature','vegetation_greenness','leaf_area_index','slope_type','soil_type','soil_color')
             call adjust_for_landice (data_mdl_one_tile, vegt_mdl_one_tile, i_mdl, j_mdl, field_names(n))
         end select
       endif
       where(mask_mdl_one_tile == 0) data_mdl_one_tile = missing
       call output (data_mdl_one_tile, lat_mdl_one_tile, lon_mdl_one_tile, i_mdl, j_mdl, tile, record, t, n)
     endif

     if (.not. fract_vegsoil_type) then
       if (field_names(n) == 'vegetation_type') then
         print*,"- CALL FieldScatter FOR MODEL GRID VEGETATION TYPE."
         call ESMF_FieldScatter(vegt_field_mdl, data_mdl_one_tile, rootPet=0, tile=tile, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldScatter.", rc)
       endif
     endif

    
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

 end subroutine interp

!> Ensure consistent fields at land ice points.
!! Land ice is vegetation type 15 (variable landice).
!! output is Model field.
!!
!! @param[in] field Model field before adjustments for land ice.
!! @param[in] vegt Vegetation type on the model tile.
!! @param[inout] idim i dimension of model tile.
!! @param[inout] jdim j dimension of model tile.
!! @param[in] field_ch Field name.
!! @author George Gayno NCEP/EMC
 subroutine adjust_for_landice(field, vegt, idim, jdim, field_ch)

 use esmf
 use mpi

 implicit none

 character(len=*), intent(in)      :: field_ch

 integer, intent(in)               :: idim, jdim

 real(esmf_kind_i4), intent(in)    :: vegt(idim,jdim)
 real(esmf_kind_r4), intent(inout) :: field(idim,jdim)

 integer, parameter                :: landice=15

 integer                           :: i, j, ierr

 real                              :: landice_value
 

 select case (field_ch)
   case ('substrate_temperature') ! soil substrate temp
     landice_value = 273.15
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = min(field(i,j), landice_value)
       endif
     enddo
     enddo
   case ('vegetation_greenness') ! vegetation greenness
     landice_value = 0.01 ! 1.0% is bare ground
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = landice_value
       endif
     enddo
     enddo
   case ('leaf_area_index') ! leaf area index
     landice_value = 0.0 ! bare ground
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = landice_value
       endif
     enddo
     enddo
   case ('slope_type') ! slope type
     landice_value = 9.0
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = landice_value
       else
         if (nint(field(i,j)) == nint(landice_value)) field(i,j) = 2.0
       endif
     enddo
     enddo
   case ('soil_type') ! soil type
     landice_value = 16.0
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = landice_value
       else
         if (nint(field(i,j)) == nint(landice_value)) field(i,j) = 6.0
       endif
     enddo
     enddo
   case ('soil_color') ! soil color
     landice_value = 10.0
     do j = 1, jdim
     do i = 1, idim
       if (nint(vegt(i,j)) == landice) then
         field(i,j) = landice_value
       endif
     enddo
     enddo
   case default
     print*,'- FATAL ERROR IN ROUTINE ADJUST_FOR_LANDICE.  UNIDENTIFIED FIELD : ', field_ch
     call mpi_abort(mpi_comm_world, 57, ierr)
 end select

 end subroutine adjust_for_landice
 
 
 
 
!> use Scale to fix the data to the correct value
!! 
!! @param[inout] idim i dimension of model tile.
!! @param[inout] jdim j dimension of model tile.
!! @param[in] field_ch Field name.
!! @author George Gayno NCEP/EMC
!! @author Sanath Kumar NCEP/EMC
!!
subroutine scale_data(field,idim,jdim,scale)

 use esmf
 use mpi

 implicit none
 integer, intent(in)               :: idim, jdim
 integer                           :: i, j
 real(esmf_kind_r4), intent(inout) :: field(idim,jdim)
 double precision                  :: scale

     do j = 1, jdim
      do i = 1, idim
         field(i,j) = field(i,j)*scale
      
      enddo
     enddo
 end subroutine scale_data

 
 
 
 
