!> @file
!! @brief Defines the model grid.
!! @author George Gayno @date 2018

!> This module defines the model grid.
!!
!! Variables named with '_mdl' refer to the model grid.
!!
!! @author George Gayno @date 2018
 module model_grid

 use esmf

 implicit none

 private

 character(len=5), allocatable, public  :: grid_tiles(:) !< Array of model grid tile names.

 integer, public               :: i_mdl !< i dimension of model tile.
 integer, public               :: j_mdl !< j dimension of model tile.
 integer, public               :: ij_mdl !< Total number of points on a model tile.
 integer, public               :: num_tiles !< Total number of model grid tiles.

 real(kind=4), public          :: missing = -999. !<Value assigned to undefined points
                                                  !! (i.e., ocean points for a land field).

 type(esmf_grid),  public      :: grid_mdl !< ESMF grid object for the model grid.
 type(esmf_field), public      :: data_field_mdl !< ESMF field object that holds the
                                                 !! data interpolated to model grid.
 type(esmf_field), public      :: land_frac_field_mdl !< ESMF field object that holds the
                                                     !! model land fraction. When running
                                                     !! with fractional grids, will be
                                                     !! between zero and one. For non-
                                                     !! fractional grids, will contain a
                                                     !! fill value.
 type(esmf_field), public      :: mask_field_mdl !< ESMF field object that holds the
                                                 !! model land mask. Equal to '1' if
                                                 !! point is partial or all land. Equal
                                                 !! to zero is point is all non-land.
 type(esmf_field), public      :: latitude_field_mdl !< ESMF field object that holds the
                                                     !! model grid latitude.
 type(esmf_field), public      :: longitude_field_mdl !< ESMF field object that holds the
                                                      !! model grid longitude.
 type(esmf_field), public      :: vegt_field_mdl !< ESMF field object that holds the
                                                 !! vegetation type on the model grid.

 public                        :: define_model_grid
 public                        :: model_grid_cleanup

 contains

!> Define model grid.
!!
!! Define the model grid from the mosaic and orography
!! files. Create the ESMF grid object for the model grid.
!!
!! @param[in] localpet this mpi task      
!! @param[in] npets total number of mpi tasks      
!! @author George Gayno @date 2018
 subroutine define_model_grid(localpet, npets)

 use esmf
 use netcdf
 use program_setup 
 use utils
 use mpi

 implicit none

 integer, intent(in)              :: localpet, npets

 character(len=500)               :: the_file

 integer                          :: error, id_dim, id_tiles, ncid
 integer                          :: id_grid_tiles, ierr
 integer                          :: extra, rc, tile
 integer, allocatable             :: decomptile(:,:)

 integer(esmf_kind_i4), allocatable :: mask_mdl_one_tile(:,:)
 integer(esmf_kind_i4), pointer   :: mask_field_mdl_ptr(:,:)
 integer(esmf_kind_i4), pointer   :: mask_mdl_ptr(:,:)

 real(esmf_kind_r4), allocatable  :: latitude_one_tile(:,:)
 real(esmf_kind_r4), allocatable  :: longitude_one_tile(:,:)
 real(esmf_kind_r4), allocatable  :: land_frac_one_tile(:,:)

!-----------------------------------------------------------------------
! Get the number of tiles from the mosaic file.
!-----------------------------------------------------------------------

 print*,'- OPEN MODEL GRID MOSAIC FILE: ',trim(mosaic_file_mdl)
 error=nf90_open(trim(mosaic_file_mdl),nf90_nowrite,ncid)
 call netcdf_err(error, "OPENING MODEL GRID MOSAIC FILE")

 print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 call netcdf_err(error, "READING NTILES ID")
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles)
 call netcdf_err(error, "READING NTILES")
 error=nf90_inq_varid(ncid, 'gridtiles', id_grid_tiles)
 call netcdf_err(error, "READING GRIDTILES ID")
 allocate(grid_tiles(num_tiles))
 grid_tiles="NULL"
 print*,"- READ TILE NAMES"
 error=nf90_get_var(ncid, id_grid_tiles, grid_tiles)
 call netcdf_err(error, "READING GRIDTILES")

 error = nf90_close(ncid)

 print*,'- NUMBER OF TILES, MODEL GRID IS ', num_tiles

 if (mod(npets,num_tiles) /= 0) then
   print*,'- FATAL ERROR: MUST RUN THIS PROGRAM WITH A TASK COUNT THAT'
   print*,'- IS A MULTIPLE OF THE NUMBER OF TILES.'
   call mpi_abort(mpi_comm_world, 44, ierr)
 endif

!-----------------------------------------------------------------------
! Get the model grid specs and land mask from the orography files.
!-----------------------------------------------------------------------

 orog_dir_mdl = trim(orog_dir_mdl) // '/'
 the_file = trim(orog_dir_mdl) // trim(orog_files_mdl(1))

 print*,'- OPEN FIRST MODEL GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, "OPENING MODEL GRID OROGRAPHY FILE")
 print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, "READING MODEL LON ID")
 error=nf90_inquire_dimension(ncid,id_dim,len=i_mdl)
 call netcdf_err(error, "READING MODEL LON")
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, "READING MODEL LAT ID")
 error=nf90_inquire_dimension(ncid,id_dim,len=j_mdl)
 call netcdf_err(error, "READING MODEL LAT")
 error = nf90_close(ncid)

 print*,"- I/J DIMENSIONS OF THE MODEL GRID TILES ", i_mdl, j_mdl

 ij_mdl = i_mdl * j_mdl

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles

 allocate(decomptile(2,num_tiles))

 do tile = 1, num_tiles
   decomptile(:,tile)=(/1,extra/)
 enddo

 print*,"- CALL GridCreateMosaic FOR MODEL GRID"
 grid_mdl = ESMF_GridCreateMosaic(filename=trim(mosaic_file_mdl), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, &
                                                   ESMF_STAGGERLOC_CORNER/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_mdl), &
                                  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", rc)

 print*,"- CALL FieldCreate FOR DATA INTERPOLATED TO MODEL GRID."
 data_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="data on model grid", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (.not. fract_vegsoil_type) then
   print*,"- CALL FieldCreate FOR VEGETATION TYPE INTERPOLATED TO MODEL GRID."
   vegt_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="veg type on model grid", &
                                   rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)
 endif

 print*,"- CALL FieldCreate FOR MODEL GRID LATITUDE."
 latitude_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="latitude on model grid", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR MODEL GRID LONGITUDE."
 longitude_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="longitude on model grid", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

!-----------------------------------------------------------------------
! Set model land mask.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR MODEL GRID LANDMASK."
 mask_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_I4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="model grid land mask", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldGet FOR MODEL GRID LANDMASK."
 nullify(mask_field_mdl_ptr)
 call ESMF_FieldGet(mask_field_mdl, &
                    farrayPtr=mask_field_mdl_ptr,  &
                    rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldCreate FOR MODEL GRID LAND FRACTION."
 land_frac_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="model grid land fraction", &
                                   rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (localpet == 0) then
   allocate(mask_mdl_one_tile(i_mdl,j_mdl))
   allocate(land_frac_one_tile(i_mdl,j_mdl))
   allocate(latitude_one_tile(i_mdl,j_mdl))
   allocate(longitude_one_tile(i_mdl,j_mdl))
 else
   allocate(mask_mdl_one_tile(0,0))
   allocate(land_frac_one_tile(0,0))
   allocate(latitude_one_tile(0,0))
   allocate(longitude_one_tile(0,0))
 endif

 do tile = 1, num_tiles
   if (localpet == 0) then
     the_file = trim(orog_dir_mdl) // trim(orog_files_mdl(tile))
     call get_model_info(trim(the_file), mask_mdl_one_tile, land_frac_one_tile, & 
                         latitude_one_tile, longitude_one_tile, i_mdl, j_mdl)
   endif

   print*,"- CALL FieldScatter FOR MODEL GRID MASK. TILE IS: ", tile
   call ESMF_FieldScatter(mask_field_mdl, mask_mdl_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

   print*,"- CALL FieldScatter FOR MODEL GRID LAND FRACTION. TILE IS: ", tile
   call ESMF_FieldScatter(land_frac_field_mdl, land_frac_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

   print*,"- CALL FieldScatter FOR MODEL LATITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_field_mdl, latitude_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

   print*,"- CALL FieldScatter FOR MODEL LONGITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_field_mdl, longitude_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo

 deallocate(mask_mdl_one_tile, latitude_one_tile, longitude_one_tile, land_frac_one_tile)

 print*,"- CALL GridAddItem FOR MODEL GRID."
 call ESMF_GridAddItem(grid_mdl, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddItem", rc)

 print*,"- CALL GridGetItem FOR MODEL GRID."
 nullify(mask_mdl_ptr)
 call ESMF_GridGetItem(grid_mdl, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_mdl_ptr, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetItem", rc)

 mask_mdl_ptr = mask_field_mdl_ptr

 end subroutine define_model_grid

!> Get model information
!!
!! Read model land/sea mask, land fraction and lat/lon from the orography file.
!!
!! @param[in] orog_file the orography file
!! @param[out] mask land/sea mask 0-all non-land; 1-some land.
!! @param[out] land_frac land fraction between 0 and 1.
!! @param[out] lat2d latitude
!! @param[out] lon2d longitude
!! @param[in] idim i dimension of the model tile
!! @param[in] jdim j dimension of the model tile
!! @author George Gayno @date 2018
 subroutine get_model_info(orog_file, mask, land_frac, lat2d, lon2d, idim, jdim)

 use esmf
 use netcdf
 use utils

 implicit none

 character(len=*), intent(in)       :: orog_file

 integer, intent(in)                :: idim, jdim
 integer(esmf_kind_i4), intent(out) :: mask(idim,jdim)

 real(esmf_kind_r4), intent(out)    :: lat2d(idim,jdim)
 real(esmf_kind_r4), intent(out)    :: lon2d(idim,jdim)
 real(esmf_kind_r4), intent(out)    :: land_frac(idim,jdim)

 integer                            :: error, lat, lon, i, j
 integer                            :: ncid, id_dim, id_var

 real(kind=4), allocatable          :: dummy(:,:)

 print*,"- READ MODEL OROGRAPHY FILE"

 print*,'- OPEN FILE: ', orog_file
 error=nf90_open(orog_file,nf90_nowrite,ncid)
 call netcdf_err(error, "OPENING MODEL OROGRAPHY FILE")

 print*,"- READ I-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, "READING LON ID")
 error=nf90_inquire_dimension(ncid,id_dim,len=lon)
 call netcdf_err(error, "READING LON")

 print*,"- READ J-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, "READING LAT ID")
 error=nf90_inquire_dimension(ncid,id_dim,len=lat)
 call netcdf_err(error, "READING LAT")

 print*,"- I/J DIMENSIONS: ", lon, lat

 if ((lon /= idim) .or. (lat /= jdim)) then
   call error_handler("MISMATCH IN DIMENSIONS.")
 endif

 allocate(dummy(idim,jdim))

!-----------------------------------------------------------------------
! If the lake maker was used, we are running with a fractional
! land/non-land grid and there will be a 'lake_frac' record.
! In that case, land/non-land is determined by 'land_frac'.
!
! If the lake maker was not used, use 'slmsk', which is defined
! as the nint(land_frac).
!
! In summary, if 'mask' is one, the point is all land or
! partial land and surface data will be mapped to it. Otherwise,
! when 'mask' is zero, then the point is all non-land and
! surface data will not be mapped to it.
!-----------------------------------------------------------------------

!error=nf90_inq_varid(ncid, 'lake_frac', id_var)
!if (error /= 0) then
!  print*,"- READ LAND MASK (SLMSK)"
!  error=nf90_inq_varid(ncid, 'slmsk', id_var)
!  call netcdf_err(error, "READING SLMSK ID")
!  error=nf90_get_var(ncid, id_var, dummy)
!  call netcdf_err(error, "READING SLMSK")
!  mask = nint(dummy)
!  land_frac = -999.
!else
   print*,"- READ LAND FRACTION"
   error=nf90_inq_varid(ncid, 'land_frac', id_var)
   call netcdf_err(error, "READING LAND_FRAC ID")
   error=nf90_get_var(ncid, id_var, land_frac)
   call netcdf_err(error, "READING LAND_FRAC")
   mask = 0
   do j = 1, lat
   do i = 1, lon
     if (land_frac(i,j) > 0.0) then
       mask(i,j) = 1
     endif
   enddo
   enddo
!endif

 print*,"- READ LATITUDE"
 error=nf90_inq_varid(ncid, 'geolat', id_var)
 call netcdf_err(error, "READING GEOLAT ID")
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, "READING GEOLAT")
 lat2d=dummy

 print*,"- READ LONGITUDE"
 error=nf90_inq_varid(ncid, 'geolon', id_var)
 call netcdf_err(error, "READING GEOLON ID")
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, "READING GEOLON")
 lon2d=dummy

 error = nf90_close(ncid)

 deallocate (dummy)

 end subroutine get_model_info

!> Model grid cleanup.
!!
!! Free up memory associated with this module.
!!
!! @author George Gayno @date 2018
 subroutine model_grid_cleanup

 implicit none

 integer         :: rc

 print*,"- CALL GridDestroy FOR MODEL GRID."
 call ESMF_GridDestroy(grid_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID LAND MASK."
 call ESMF_FieldDestroy(mask_field_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID LAND FRACTION."
 call ESMF_FieldDestroy(land_frac_field_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID DATA FIELD."
 call ESMF_FieldDestroy(data_field_mdl,rc=rc)

 if (ESMF_FieldIsCreated(vegt_field_mdl)) then
   print*,"- CALL FieldDestroy FOR MODEL GRID VEGETATION TYPE."
   call ESMF_FieldDestroy(vegt_field_mdl,rc=rc)
 endif

 print*,"- CALL FieldDestroy FOR MODEL GRID LATITUDE."
 call ESMF_FieldDestroy(latitude_field_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID LONGITUDE."
 call ESMF_FieldDestroy(longitude_field_mdl,rc=rc)

 end subroutine model_grid_cleanup

 end module model_grid
