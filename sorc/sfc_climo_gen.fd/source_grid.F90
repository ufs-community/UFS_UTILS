!> @file
!! @brief Read grid specs, date information and land/sea mask for
!! the source data that will be interpolated to the model grid.
!! @author George Gayno @date 2018

!> Read grid specs, date information and land/sea mask for
!! the source data that will be interpolated to the model grid.
!! Also, sets up the ESMF grid object for the source grid.
!! Source grid is assumed to be global lat/lon.
!!
!! @author George Gayno @date 2018
module source_grid

 use esmf
 use utils

 implicit none

 private 

 character(len=50), allocatable, public :: field_names(:) !< Names of fields to be processed.
 character(len=75), public              :: source !< Original source of the data.

 integer, public               :: i_src !< i dimension of the source grid.
 integer, public               :: j_src !< j dimension of the source grid.
 integer, public               :: num_records !< Number of fields times time records.
 integer, public               :: num_time_recs !< Number of time records.
 integer, public               :: num_fields !< Number of fields in the file. Some
                                             !! files have more than one (ex: 
                                             !! the 4-component albedo).
 integer, allocatable, public  :: day_of_rec(:) !< Day of each time record with
                                                !! respect to Jan 1.

 type(esmf_grid), public       :: grid_src !< ESMF grid object for the source grid.

 public                        :: define_source_grid
 public                        :: source_grid_cleanup

 contains

 !> Defines esmf grid object for source grid. Retrieves date and field
 !! information from source file.
 !!
 !! Read date information from input source data file.
 !! Create esmf grid object for the source grid.
 !!
 !! @param[in] localpet mpi task number
 !! @param[in] npets total number mpi tasks
 !! @param[in] input_file file containing the source grid data.
 !! @author George Gayno @date 2018
 subroutine define_source_grid(localpet, npets, input_file)

 use mpi
 use netcdf

 implicit none

 character(len=*), intent(in)     :: input_file

 integer, intent(in)              :: localpet, npets
 
 character(len=50)                :: vname

 integer                          :: dimid, dims(1), ncid, status
 integer                          :: count, num_vars, n
 integer                          :: rc, varid, i, j, i_src_corner
 integer(esmf_kind_i4), pointer   :: mask_ptr(:,:)
 integer                          :: clb(2), cub(2)
 integer                          :: clb_corner(2), cub_corner(2)

 real(esmf_kind_r4), allocatable  :: mask_global(:,:)
 real(esmf_kind_r8), allocatable  :: lat_global(:)
 real(esmf_kind_r8), allocatable  :: lon_global(:)
 real(esmf_kind_r8), allocatable  :: lat_corner_global(:)
 real(esmf_kind_r8), allocatable  :: lon_corner_global(:)
 real(esmf_kind_r4), pointer      :: mask_field_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lat_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lat_corner_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_corner_ptr(:,:)
 real                             :: lon_extent
 real(esmf_kind_r4)               :: missing

 type(esmf_field)                 :: mask_field
 type(esmf_polekind_flag)         :: polekindflag(2)

 print*,"- OPEN FILE: ", trim(input_file)
 status = nf90_open(input_file, nf90_nowrite, ncid)
 call netcdf_err(status, "OPENING INPUT SOURCE FILE")

 status = nf90_get_att(ncid, nf90_global, 'source', source)
 if (status /= nf90_noerr) source="unknown"

 if(localpet == 0) print*,'- SOURCE OF DATA IS: ', trim(source)

 status = nf90_inq_dimid(ncid, 'idim', dimid)
 call netcdf_err(status, "READING I DIMENSION ID OF SOURCE DATA")
 status = nf90_inquire_dimension(ncid, dimid, len=i_src)
 call netcdf_err(status, "READING I DIMENSION OF SOURCE DATA")

 status = nf90_inq_dimid(ncid, 'jdim', dimid)
 call netcdf_err(status, "READING J DIMENSION ID OF SOURCE DATA")
 status = nf90_inquire_dimension(ncid, dimid, len=j_src)
 call netcdf_err(status, "READING J DIMENSION OF SOURCE DATA")

 if(localpet == 0) print*,'- I/J DIMENSION OF SOURCE DATA: ', i_src, j_src

 allocate(lat_global(j_src))
 status = nf90_inq_varid(ncid, 'lat', varid)
 call netcdf_err(status, "READING SOURCE DATA LATITUDE ID")
 status = nf90_get_var(ncid, varid, lat_global)
 call netcdf_err(status, "READING SOURCE DATA LATITUDES")

 allocate(lon_global(i_src))
 status = nf90_inq_varid(ncid, 'lon', varid)
 call netcdf_err(status, "READING SOURCE DATA LONGITUDE ID")
 status = nf90_get_var(ncid, varid, lon_global)
 call netcdf_err(status, "READING SOURCE DATA LONGITUDE")

 allocate(lat_corner_global(j_src+1))
 status = nf90_inq_varid(ncid, 'lat_corner', varid)
 call netcdf_err(status, "READING SOURCE DATA CORNER LATITUDE ID")
 status = nf90_get_var(ncid, varid, lat_corner_global)
 call netcdf_err(status, "READING SOURCE DATA CORNER LATITUDE")

!-----------------------------------------------------------------------
! Dimension of lon_corner depends on whether input data is periodic
! (global) or regional.
!-----------------------------------------------------------------------

 status = nf90_inq_varid(ncid, 'lon_corner', varid)
 call netcdf_err(status, "READING SOURCE DATA CORNER LONGITUDE ID")
 status = nf90_inquire_variable(ncid, varid, dimids=dims)
 call netcdf_err(status, "READING SOURCE DATA CORNER LONGITUDE DIMIDS")
 status = nf90_inquire_dimension(ncid, dims(1), len=i_src_corner)
 call netcdf_err(status, "READING SOURCE DATA CORNER LONGITUDE DIMS")
 allocate(lon_corner_global(i_src_corner))
 status = nf90_get_var(ncid, varid, lon_corner_global)
 call netcdf_err(status, "READING SOURCE DATA CORNER LONGITUDE")

 status = nf90_inq_dimid(ncid, 'time', dimid)
 call netcdf_err(status, "READING SOURCE DATA TIME ID")
 status = nf90_inquire_dimension(ncid, dimid, len=num_time_recs)
 call netcdf_err(status, "READING SOURCE DATA NUM TIME PERIODS")

 allocate(day_of_rec(num_time_recs))
 status = nf90_inq_varid(ncid, 'time', varid)
 call netcdf_err(status, "READING SOURCE DATA TIME VARID")
 status = nf90_get_var(ncid, varid, day_of_rec)
 call netcdf_err(status, "READING SOURCE DATA DAY OF RECORD")

 print*,'- SOURCE DATA DAYS OF RECORD: ', day_of_rec

 status = nf90_inquire(ncid, nVariables=num_vars)
 call netcdf_err(status, "READING NUMBER OF VARIABLES SOURCE DATA")

!-----------------------------------------------------------------------
! Assumes files contain records for time, lat, lon, lat_corner, 
! lon_corner.  So number of fields processed will be the total
! number of variables minus 5.  
!-----------------------------------------------------------------------

 num_fields = num_vars - 5
 num_records = num_vars * num_time_recs

 allocate(field_names(num_fields))

 count = 0
 do n = 1, num_vars

  status = nf90_inquire_variable(ncid, n, name=vname)
  call netcdf_err(status, "READING FIELD NAMES")

  if (trim(vname) == 'time') cycle
  if (trim(vname) == 'lon') cycle
  if (trim(vname) == 'lon_corner') cycle
  if (trim(vname) == 'lat') cycle
  if (trim(vname) == 'lat_corner') cycle

  count = count + 1
  field_names(count) = vname

 enddo

 if(localpet==0) print*,'- FIELDS TO BE PROCESSED: ', field_names

 if (localpet == 0) then
   allocate(mask_global(i_src,j_src))
   status = nf90_inq_varid(ncid, field_names(1), varid)
   call netcdf_err(status, "READING FIELD ID")
   status = nf90_get_var(ncid, varid, mask_global)
   call netcdf_err(status, "READING FIELD")
 else
   allocate(mask_global(0,0))
 endif

!--------------------------------------------------------------------------
! Read in missing value.  This is used to mask out data at non-land
! points.
!--------------------------------------------------------------------------

 status = nf90_inq_varid(ncid, field_names(1), varid)
 call netcdf_err(status, "READING FIELD 1 ID")
 status=nf90_get_att(ncid, varid, 'missing_value', missing)
 call netcdf_err(status, "READING MISSING VALUE")

 status = nf90_close(ncid)

!--------------------------------------------------------------------------
! Create ESMF grid object for the source data grid.  Check if 
! data is periodic in the east/west direction.
!
! Note: When using regional data, there is always the chance of
! the target grid being located outside the input grid.  ESMF
! support recommends passing back the dstField (esmf_typekind_i4) from
! all calls to FieldRegridStore and checking for the
! ESMF_REGRIDSTATUS_OUTSIDE flag.
!--------------------------------------------------------------------------

 lon_extent = mod((lon_global(i_src)-lon_global(1))-1+3600,360.)+1.0
 
 if (lon_extent < 350.0) then

   print*,"- CALL GridCreateNoPeriDim FOR SOURCE GRID"
   grid_src = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_src,j_src/), &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    name="source_grid", &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridCreateNoPeriDim.", rc)

 else

   polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

   print*,"- CALL GridCreate1PeriDim FOR SOURCE GRID"
   grid_src = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_src,j_src/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    name="source_grid", &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridCreate1PeriDim.", rc)

 endif

 print*,"- CALL FieldCreate FOR SOURCE GRID LANDMASK."
 mask_field = ESMF_FieldCreate(grid_src, &
                               typekind=ESMF_TYPEKIND_R4, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, &
                               name="source grid land mask", &
                               rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate.", rc)

 print*,"- CALL FieldScatter FOR SOURCE GRID MASK."
 call ESMF_FieldScatter(mask_field, mask_global, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter.", rc)

 print*,"- CALL GridAddItem FOR SOURCE GRID MASK."
 call ESMF_GridAddItem(grid_src, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddItem.", rc)

 print*,"- CALL GridGetItem FOR SOURCE GRID MASK."
 nullify(mask_ptr)
 call ESMF_GridGetItem(grid_src, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_ptr, &
                       totalLBound=clb,  &
                       totalUBound=cub, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetItem", rc)
 
 print*,"- CALL FieldGet FOR SOURCE GRID LANDMASK."
 nullify(mask_field_ptr)
 call ESMF_FieldGet(mask_field, &
                    farrayPtr=mask_field_ptr,  &
                    rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet.", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if ( abs(mask_field_ptr(i,j)-missing) < 0.001) then
       mask_ptr(i,j) = 0
     else
       mask_ptr(i,j) = 1
     endif
   enddo
 enddo

 deallocate(mask_global)

 print*,"- CALL FieldDestroy FOR SOURCE GRID LAND MASK."
 call ESMF_FieldDestroy(mask_field,rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldDestroy.", rc)

! Set lat/lons of grid points

 print*,"- CALL GridAddCoord FOR SOURCE GRID CENTER LOCATION."
 call ESMF_GridAddCoord(grid_src, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord.", rc)

 print*,"- CALL GridGetCoord FOR SOURCE GRID CENTER LONGITUDE."
 nullify(lon_ptr)
 call ESMF_GridGetCoord(grid_src, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord.", rc)

 print*,"- CALL GridGetCoord FOR SOURCE GRID CENTER LATITUDE."
 nullify(lat_ptr)
 call ESMF_GridGetCoord(grid_src, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        farrayPtr=lat_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord.", rc)

 do j = clb(2), cub(2)
   lat_ptr(:,j) = lat_global(j)
 enddo

 do i = clb(1), cub(1)
   lon_ptr(i,:) = lon_global(i)
 enddo

 print*,"- CALL GridAddCoord FOR SOURCE GRID CORNER LOCATION."
 call ESMF_GridAddCoord(grid_src, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, &
                        rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord.", rc)

 print*,"- CALL GridGetCoord FOR SOURCE GRID CORNER LONGITUDE."
 nullify(lon_corner_ptr)
 call ESMF_GridGetCoord(grid_src, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_ptr, &
                        rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord.", rc)

 print*,"- CALL GridGetCoord FOR SOURCE GRID CORNER LATITUDE."
 nullify(lat_corner_ptr)
 call ESMF_GridGetCoord(grid_src, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb_corner, &
                        computationalUBound=cub_corner, &
                        farrayPtr=lat_corner_ptr, &
                        rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord.", rc)

 do j = clb_corner(2), cub_corner(2)
   lat_corner_ptr(:,j) = lat_corner_global(j)
 enddo

 do i = clb_corner(1), cub_corner(1)
   lon_corner_ptr(i,:) = lon_corner_global(i)
 enddo

 deallocate(lat_global)
 deallocate(lon_global)
 deallocate(lat_corner_global)
 deallocate(lon_corner_global)

 end subroutine define_source_grid

 !> Free up memory associated with this module.
 !!
 !! @author George Gayno @date 2018 
 subroutine source_grid_cleanup

 implicit none

 integer  :: rc

 print*,"- CALL GridDestroy FOR SOURCE GRID."
 call ESMF_GridDestroy(grid_src,rc=rc)

 deallocate(day_of_rec)
 deallocate(field_names)

 end subroutine source_grid_cleanup

 end module source_grid
