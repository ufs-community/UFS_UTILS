 module model_grid

!--------------------------------------------------------------------------
! module documentation block
!
! Module: model grid
!   pgrmmr: gayno           org: w/np2           date: 2018
!
! Abstract: Defines the model grid.
!
! Usage:  use model_grid 
!
! Public Subroutines:
! -------------------
! define_model_grid            Defines esmf grid object for the
!                              model grid.
! model_grid_cleanup           Free up memory used in this module.
!
! Public variables:
! -----------------
!
! Variables named with 'mdl' refer to the model grid.
!
! data_field_mdl               ESMF field object that holds the
!                              data interpolated to model grid.
! grid_mdl                     ESMF grid object for the model grid.
! grid_tiles                   Array of model grid tile names.
! i/j_mdl                      i/j dimensions of model tile.
! mdl_field_mdl                ESMF field object that holds the
!                              model land mask.
! missing                      Value assigned to undefined points
!                              (i.e., ocean points for a land
!                              field).
! num_tiles                    Total number of model grid tiles.
! vegt_field_mdl               ESMF field object that holds the
!                              vegetation type on the model grid.
! 
!-----------------------------------------------------------------------

 use esmf

 implicit none

 private

 character(len=5), allocatable, public  :: grid_tiles(:)

 integer, public               :: i_mdl, j_mdl, ij_mdl, num_tiles

 real(kind=4), public          :: missing = -999.

 type(esmf_grid),  public      :: grid_mdl
 type(esmf_field), public      :: data_field_mdl, mask_field_mdl
 type(esmf_field), public      :: vegt_field_mdl

 public                        :: define_model_grid
 public                        :: model_grid_cleanup

 contains

 subroutine define_model_grid(localpet, npets)

!-----------------------------------------------------------------------
!  subroutine documentation block
!
! Subroutine: define model grid
!   prgmmr: gayno          org: w/np2           date: 2018
!
! Abstract: Define the model grid from the mosaic and orography
!   files.  Create the ESMF grid object for the model grid.
!
! Usage:  define_model_grid(localpet, npets)
!
!   input argument list:
!     localpet               this mpi task      
!     npets                  total number of mpi tasks      
!-----------------------------------------------------------------------

 use esmf
 use netcdf
 use program_setup 
 use utils

 implicit none

 integer, intent(in)              :: localpet, npets

 character(len=500)               :: the_file

 integer                          :: error, id_dim, id_tiles, ncid
 integer                          :: id_grid_tiles
 integer                          :: extra, rc, tile
 integer, allocatable             :: decomptile(:,:)

 integer, allocatable             :: mask_mdl_one_tile(:,:)
 integer(esmf_kind_i4), pointer   :: mask_field_mdl_ptr(:,:)
 integer(esmf_kind_i4), pointer   :: mask_mdl_ptr(:,:)

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
   call mpi_abort
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

 print*,"- CALL FieldCreate FOR VEGETATION TYPE INTERPOLATED TO MODEL GRID."
 vegt_field_mdl = ESMF_FieldCreate(grid_mdl, &
                                   typekind=ESMF_TYPEKIND_R4, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="veg type on model grid", &
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

 if (localpet == 0) then
   allocate(mask_mdl_one_tile(i_mdl,j_mdl))
 else
   allocate(mask_mdl_one_tile(0,0))
 endif

 do tile = 1, num_tiles
   if (localpet == 0) then
     the_file = trim(orog_dir_mdl) // trim(orog_files_mdl(tile))
     call get_model_mask(trim(the_file), mask_mdl_one_tile, i_mdl, j_mdl)
   endif
   print*,"- CALL FieldScatter FOR MODEL GRID MASK. TILE IS: ", tile
   call ESMF_FieldScatter(mask_field_mdl, mask_mdl_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 deallocate(mask_mdl_one_tile)

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

 subroutine get_model_mask(orog_file, mask, idim, jdim)

!-----------------------------------------------------------------------
!  subroutine documentation block
!
! Subroutine: get model mask
!   prgmmr: gayno          org: w/np2           date: 2018
!
! Abstract: Read model land/sea mask from the orography file.
!
! Usage:  call get_model_mask(orog_file, mask, idim, jdim)
!
!   input argument list:
!     orog_file              the orography file
!     i/jdim                 i/j dimension of the model tile
!
!   output argument list:
!     mask                   land/sea mask
!
!-----------------------------------------------------------------------

 use esmf
 use netcdf
 use utils

 implicit none

 character(len=*), intent(in)       :: orog_file

 integer, intent(in)                :: idim, jdim
 integer(esmf_kind_i4), intent(out) :: mask(idim,jdim)

 integer                            :: error, lat, lon
 integer                            :: ncid, id_dim, id_var

 real(kind=4), allocatable          :: dummy(:,:)

 print*,"- READ MODEL LAND MASK FILE"

 print*,'- OPEN LAND MASK FILE: ', orog_file
 error=nf90_open(orog_file,nf90_nowrite,ncid)
 call netcdf_err(error, "OPENING MODEL LAND MASK FILE")

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

 print*,"- READ LAND MASK"
 error=nf90_inq_varid(ncid, 'slmsk', id_var)
 call netcdf_err(error, "READING SLMSK ID")
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, "READING SLMSK")

 error = nf90_close(ncid)

 mask = nint(dummy)

 deallocate (dummy)

 end subroutine get_model_mask

 subroutine model_grid_cleanup

!-----------------------------------------------------------------------
!  subroutine documentation block
!
! Subroutine: model grid cleanup
!   prgmmr: gayno          org: w/np2           date: 2018
!
! Abstract: Free up memory associated with this module.
!
! Usage:  call model_grid_cleanup
!
!-----------------------------------------------------------------------

 implicit none

 integer         :: rc

 print*,"- CALL GridDestroy FOR MODEL GRID."
 call ESMF_GridDestroy(grid_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID LAND MASK."
 call ESMF_FieldDestroy(mask_field_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID DATA FIELD."
 call ESMF_FieldDestroy(data_field_mdl,rc=rc)

 print*,"- CALL FieldDestroy FOR MODEL GRID VEGETATION TYPE."
 call ESMF_FieldDestroy(vegt_field_mdl,rc=rc)

 end subroutine model_grid_cleanup

 end module model_grid
