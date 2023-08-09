!> @file
!! @brief Specify input and target model grids.
!! @author George Gayno NCEP/EMC

!> Sets up the ESMF grid objects for the input data grid and target
!! FV3 grid.
!!
!! @author George Gayno NCEP/EMC
 module model_grid

 use esmf
 use ESMF_LogPublicMod

 use utilities, only                    : error_handler, netcdf_err
 implicit none

 private

 character(len=5), allocatable, public  :: tiles_target_grid(:)
                                           !< Tile names of target grid.
 character(len=50), public              :: input_grid_type = "latlon"
                                           !< map projection of input grid

 ! Made lsoil_target non-parameter to allow for RAP land surface initiation
 integer, public                        :: lsoil_target = 4 ! # soil layers
                                           !< Number of soil layers, target grid.
 integer, public                        :: i_input
                                           !< i-dimension of input grid
                                           !! (or of each global tile)
 integer, public                        :: j_input
                                           !< j-dimension of input grid
                                           !! (or of each global tile)
 integer, public                        :: ip1_input
                                           !< i_input plus 1
 integer, public                        :: jp1_input
                                           !< j_input plus 1
 integer, public                        :: i_target
                                           !< i dimension of each global tile, 
                                           !! or of a nest, target grid.
 integer, public                        :: j_target
                                           !< j dimension of each global tile,
                                           !! or of a nest, target grid.
 integer, public                        :: ip1_target
                                           !< ip1_target plus 1
 integer, public                        :: jp1_target
                                           !< jp1_target plus 1
 integer, public                        :: num_tiles_input_grid
                                           !< Number of tiles, input grid
 integer, public                        :: num_tiles_target_grid
                                           !< Number of tiles, target grid

 type(esmf_grid),  public               :: input_grid
                                           !< input grid esmf grid object
 type(esmf_grid),  public               :: target_grid
                                           !< target grid esmf grid object.

 type(esmf_field),  public              :: latitude_input_grid
                                           !< latitude of grid center, input grid
 type(esmf_field),  public              :: longitude_input_grid
                                           !< longitude of grid center, input grid
 type(esmf_field),  public              :: latitude_s_input_grid
                                           !< latitude of 'south' edge of grid
                                           !! box, input grid
 type(esmf_field),  public              :: longitude_s_input_grid
                                           !< longitude of 'south' edge of grid
                                           !! box, input grid
 type(esmf_field),  public              :: latitude_w_input_grid
                                           !< latitude of 'west' edge of grid
                                           !! box, input grid
 type(esmf_field),  public              :: longitude_w_input_grid
                                           !< longitude of 'west' edge of grid
                                           !! box, input grid

 type(esmf_field),  public              :: landmask_target_grid
                                           !< land mask target grid - '1' land;
                                           !! '0' non-land
 type(esmf_field),  public              :: land_frac_target_grid
                                           !< land fraction, target grid
 type(esmf_field),  public              :: latitude_target_grid
                                           !< latitude of grid center, target grid
 type(esmf_field),  public              :: latitude_s_target_grid
                                           !< latitude of 'south' edge of grid
                                           !! box, target grid
 type(esmf_field),  public              :: latitude_w_target_grid
                                           !< latitude of 'west' edge of grid
                                           !! box, target grid
 type(esmf_field),  public              :: longitude_target_grid
                                           !< longitude of grid center, target grid
 type(esmf_field),  public              :: longitude_s_target_grid
                                           !< longitude of 'south' edge of grid
                                           !! box, target grid
 type(esmf_field),  public              :: longitude_w_target_grid
                                           !< longitude of 'west' edge of grid
                                           !! box, target grid
 type(esmf_field),  public              :: seamask_target_grid
                                           !< sea mask target grid - '1' non-land;
                                           !! '0' land
 type(esmf_field),  public              :: terrain_target_grid
                                           !< terrain height target grid

 public :: define_target_grid
 public :: define_input_grid
 public :: cleanup_input_target_grid_data

 contains

!> Driver routine to setup the esmf grid object for the input grid.
!!
!! If the input source is tiled fv3 restart or history data, the grid
!! is created by reading the mosaic and grid files.  If the input
!! source is fv3 global gaussian nemsio, spectral gfs global gaussian
!! nemsio, or spectral gfs global gaussian sigio/sfcio, the grid is
!! setup by computing lat/lons using the sp library.
!!
!! @param [in] localpet ESMF local persistent execution thread 
!! @param [in] npets  Number of persistent execution threads
!! @author George Gayno NCEP/EMC   
 subroutine define_input_grid(localpet, npets)

 use program_setup, only       : input_type

 implicit none

 integer, intent(in)          :: localpet, npets

 if (trim(input_type) == "gaussian_nemsio" .or. &
     trim(input_type) == "gfs_gaussian_nemsio" .or. &
     trim(input_type) == "gfs_sigio" .or. &
     trim(input_type) == "gaussian_netcdf") then
   call define_input_grid_gaussian(npets)
 elseif (trim(input_type) == "grib2") then
   call define_input_grid_grib2(npets)
 else
   call define_input_grid_mosaic(localpet, npets)
 endif

 end subroutine define_input_grid

!> Define grid object for input data on global gaussian grids.
!!
!! Recognized file formats: 
!!  - fv3gfs nemsio
!!  - spectral gfs nemsio (starting July 19, 2017)
!!  - spectral gfs sigio  (prior to July 19, 2017)
!!  - spectral gfs sfcio  (prior to July 19, 2017)
!!
!! @param [in] npets  Number of  persistent execution threads.
!! @author George Gayno NCEP/EMC   
 subroutine define_input_grid_gaussian(npets)

 use nemsio_module

 use program_setup, only       : data_dir_input_grid, &
                                 atm_files_input_grid, &
                                 sfc_files_input_grid, &
                                 input_type, &
                                 convert_atm, convert_sfc

 use sfcio_module
 use sigio_module
 use netcdf

 implicit none

 integer, intent(in)              :: npets

 character(len=250)               :: the_file

 integer                          :: i, j, rc, clb(2), cub(2), ncid, id_grid
 integer(sfcio_intkind)           :: rc2
 integer(sigio_intkind)           :: rc3

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)
 real(esmf_kind_r8), pointer      :: lat_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lat_corner_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_corner_src_ptr(:,:)
 real(esmf_kind_r8)               :: deltalon
 real(esmf_kind_r8), allocatable  :: slat(:), wlat(:)

 type(nemsio_gfile)               :: gfile
 type(esmf_polekind_flag)         :: polekindflag(2)
 type(sfcio_head)                 :: sfchead
 type(sigio_head)                 :: sighead

 print*,"- DEFINE INPUT GRID OBJECT FOR GAUSSIAN DATA."

 num_tiles_input_grid = 1

 if (convert_sfc) then
   the_file=trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))
 elseif (convert_atm) then
   the_file=trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))
 endif

 if (trim(input_type) == "gfs_sigio") then  ! sigio/sfcio format, used by
                                               ! spectral gfs prior to 7/19/2017.

   if (convert_sfc) then   ! sfcio format
     print*,"- OPEN AND READ ", trim(the_file)
     call sfcio_sropen(21, trim(the_file), rc2)
     if (rc2 /= 0) call error_handler("OPENING FILE", rc2)
     call sfcio_srhead(21, sfchead, rc2)
     if (rc2 /= 0) call error_handler("READING FILE", rc2)
     call sfcio_sclose(21, rc2)
     i_input = sfchead%lonb
     j_input = sfchead%latb
   elseif (convert_atm) then ! sigio format
     print*,"- OPEN AND READ ", trim(the_file)
     call sigio_sropen(21, trim(the_file), rc3)
     if (rc3 /= 0) call error_handler("OPENING FILE", rc3)
     call sigio_srhead(21, sighead, rc3)
     if (rc3 /= 0) call error_handler("READING FILE", rc3)
     call sigio_sclose(21, rc3)
     i_input = sighead%lonb
     j_input = sighead%latb
   endif

 elseif (trim(input_type) == "gaussian_netcdf") then

   print*,'- OPEN AND READ: ',trim(the_file)
   rc=nf90_open(trim(the_file),nf90_nowrite,ncid)
   call netcdf_err(rc, 'opening file')

   print*,"- READ grid_xt"
   rc=nf90_inq_dimid(ncid, 'grid_xt', id_grid)
   call netcdf_err(rc, 'reading grid_xt id')
   rc=nf90_inquire_dimension(ncid,id_grid,len=i_input)
   call netcdf_err(rc, 'reading grid_xt')

   print*,"- READ grid_yt"
   rc=nf90_inq_dimid(ncid, 'grid_yt', id_grid)
   call netcdf_err(rc, 'reading grid_yt id')
   rc=nf90_inquire_dimension(ncid,id_grid,len=j_input)
   call netcdf_err(rc, 'reading grid_yt')

   rc = nf90_close(ncid)

 else ! nemsio format

   call nemsio_init(iret=rc)

   print*,"- OPEN AND READ ", trim(the_file)
   call nemsio_open(gfile, the_file, "read", iret=rc)
   if (rc /= 0) call error_handler("OPENING FILE", rc)

   call nemsio_getfilehead(gfile, iret=rc, dimx=i_input, dimy=j_input)
   if (rc /= 0) call error_handler("READING FILE", rc)

   call nemsio_close(gfile)
 
 endif

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_input,j_input/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN GridCreate1PeriDim", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 allocate(longitude(i_input,j_input))
 allocate(latitude(i_input,j_input))

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

 deallocate(slat, wlat)

 print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE."
 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE."
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_src_ptr(i,j) = longitude(i,j)
     if (lon_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_src_ptr(i,j) = lon_src_ptr(i,j) - 360.0_esmf_kind_r8
     lat_src_ptr(i,j) = latitude(i,j)
   enddo
 enddo

 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_corner_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_corner_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_corner_src_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
     if (lon_corner_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_corner_src_ptr(i,j) = lon_corner_src_ptr(i,j) - 360.0_esmf_kind_r8
     if (j == 1) then 
       lat_corner_src_ptr(i,j) = 90.0_esmf_kind_r8
       cycle
     endif
     if (j == jp1_input) then
       lat_corner_src_ptr(i,j) = -90.0_esmf_kind_r8
       cycle
     endif
     lat_corner_src_ptr(i,j) = 0.5_esmf_kind_r8 * (latitude(i,j-1)+ latitude(i,j))
   enddo
 enddo

 deallocate(latitude,longitude)

 end subroutine define_input_grid_gaussian

!> Define input grid for tiled data using the 'mosaic',
!! 'grid' and orography files.
!!
!! @param localpet ESMF local persistent execution thread 
!! @param npets Total number of persistent execution threads
!! @author George Gayno NCEP/EMC   
 subroutine define_input_grid_mosaic(localpet, npets)

 use netcdf
 use program_setup, only       : mosaic_file_input_grid,  &
                                 orog_dir_input_grid, &
                                 orog_files_input_grid

 implicit none

 character(len=500)           :: the_file

 integer, intent(in)          :: localpet, npets

 integer                      :: id_tiles, id_dim, tile
 integer                      :: extra, error, ncid
 integer, allocatable         :: decomptile(:,:)

 real(esmf_kind_r8), allocatable       :: latitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_w_one_tile(:,:)

 print*,'- OPEN INPUT GRID MOSAIC FILE: ',trim(mosaic_file_input_grid)
 error=nf90_open(trim(mosaic_file_input_grid),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening grid mosaic file')

 print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 call netcdf_err(error, 'reading ntiles id')
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles_input_grid)
 call netcdf_err(error, 'reading ntiles')

 error = nf90_close(ncid)

 print*,'- NUMBER OF TILES, INPUT MODEL GRID IS ', num_tiles_input_grid

 if (mod(npets,num_tiles_input_grid) /= 0) then
   call error_handler("MUST RUN WITH A TASK COUNT THAT IS A MULTIPLE OF 6.", 1)
 endif

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles_input_grid

 allocate(decomptile(2,num_tiles_input_grid))

 do tile = 1, num_tiles_input_grid
   decomptile(:,tile)=(/1,extra/)
 enddo

 print*,"- CALL GridCreateMosaic FOR INPUT MODEL GRID"
 input_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file_input_grid), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_input_grid), &
                                  rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", error)

!-----------------------------------------------------------------------
! Read the mask and lat/lons.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE_S."
 latitude_s_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="input_grid_latitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE_S."
 longitude_s_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="input_grid_longitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE_W."
 latitude_w_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="input_grid_latitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE_W."
 longitude_w_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="input_grid_longitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 the_file = trim(orog_dir_input_grid) // trim(orog_files_input_grid(1))

 print*,'- OPEN FIRST INPUT GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening ororgraphy file')
 print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, 'reading lon id')
 error=nf90_inquire_dimension(ncid,id_dim,len=i_input)
 call netcdf_err(error, 'reading lon')
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, 'reading lat id')
 error=nf90_inquire_dimension(ncid,id_dim,len=j_input)
 call netcdf_err(error, 'reading lat')
 error = nf90_close(ncid)

 print*,"- I/J DIMENSIONS OF THE INPUT GRID TILES ", i_input, j_input

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 if (localpet == 0) then
   allocate(longitude_one_tile(i_input,j_input))
   allocate(longitude_s_one_tile(i_input,jp1_input))
   allocate(longitude_w_one_tile(ip1_input,j_input))
   allocate(latitude_one_tile(i_input,j_input))
   allocate(latitude_s_one_tile(i_input,jp1_input))
   allocate(latitude_w_one_tile(ip1_input,j_input))
 else
   allocate(longitude_one_tile(0,0))
   allocate(longitude_s_one_tile(0,0))
   allocate(longitude_w_one_tile(0,0))
   allocate(latitude_one_tile(0,0))
   allocate(latitude_s_one_tile(0,0))
   allocate(latitude_w_one_tile(0,0))
 endif

 do tile = 1, num_tiles_input_grid
   if (localpet == 0) then
     call get_model_latlons(mosaic_file_input_grid, orog_dir_input_grid, num_tiles_input_grid, tile, &
                            i_input, j_input, ip1_input, jp1_input, latitude_one_tile, &
                            latitude_s_one_tile, latitude_w_one_tile, longitude_one_tile, &
                            longitude_s_one_tile, longitude_w_one_tile)
   endif
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_input_grid, latitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_input_grid, longitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_s_input_grid, latitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_s_input_grid, longitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_w_input_grid, latitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_w_input_grid, longitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", error)
 enddo

 deallocate(longitude_one_tile)
 deallocate(longitude_s_one_tile)
 deallocate(longitude_w_one_tile)
 deallocate(latitude_one_tile)
 deallocate(latitude_s_one_tile)
 deallocate(latitude_w_one_tile)

 end subroutine define_input_grid_mosaic

!> Define input grid object for grib2 input data.
!!
!! @param [in] npets  Number of persistent execution threads
!! @author Larissa Reames
!! @author Jeff Beck
!! @author George Gayno
 subroutine define_input_grid_grib2(npets)

 use grib_mod
 use gdswzd_mod
 use program_setup, only       : grib2_file_input_grid, data_dir_input_grid

 implicit none

 integer, intent(in)              :: npets

 character(len=500)               :: the_file

 integer                          :: i, j, k, jdisc, jgdtn, jpdtn, lugb, lugi
 integer                          :: jids(200), jgdt(200), jpdt(200), rc
 integer                          :: kgds(200), nret, clb(2), cub(2)

 logical                          :: unpack

 real                             :: res
 real, allocatable                :: rlon(:,:),rlat(:,:),xpts(:,:),ypts(:,:)
 real, allocatable                :: rlon_corner(:,:),rlat_corner(:,:)
 real, allocatable                :: rlon_diff(:,:),rlat_diff(:,:)
 real, allocatable                :: xpts_corner(:,:),ypts_corner(:,:)
 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)
 real(esmf_kind_r8), allocatable  :: latitude_corner(:,:)
 real(esmf_kind_r8), allocatable  :: longitude_corner(:,:)
 real(esmf_kind_r8), pointer      :: lat_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lat_corner_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_src_ptr(:,:)
 real(esmf_kind_r8), pointer      :: lon_corner_src_ptr(:,:)

 type(esmf_polekind_flag)         :: polekindflag(2)

 type(gribfield)                  :: gfld

 the_file = trim(data_dir_input_grid) // "/" // grib2_file_input_grid

 lugb=12

 print*,"- OPEN AND READ INPUT DATA GRIB2 FILE: ", trim(the_file)
 call baopenr(lugb,the_file,rc)
 if (rc /= 0) call error_handler("OPENING FILE", rc)

! Read the first record and get the grid definition template.

 j       = 0      ! Search at beginning of file
 lugi    = 0      ! No grib index file
 jdisc   = -1     ! Search for any discipline
 jpdtn   = -1     ! Search for any product definition template number
 jgdtn   = -1     ! Search for any grid definition template number
 jids    = -9999  ! Array of values in identification section, set to wildcard.
 jgdt    = -9999  ! Array of values in grid definition template, set to wildcard.
 jpdt    = -9999  ! Array of values in product definition template, set to wildcard.
 unpack  = .false. ! unpack data
   
 call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
 if (rc /= 0) call error_handler("DEGRIBBING INPUT FILE.", rc)

 call baclose(lugb,rc)

 if (gfld%igdtnum == 0) then
   print*,"- INPUT DATA ON LAT/LON GRID."
   input_grid_type = 'latlon'
 elseif (gfld%igdtnum == 30) then
   print*,"- INPUT DATA ON LAMBERT CONFORMAL GRID."
   input_grid_type = 'lambert'
 elseif (gfld%igdtnum == 32769) then
   print*,"- INPUT DATA ON ROTATED LAT/LON GRID."
   input_grid_type = 'rotated_latlon'
 else
   call error_handler("INPUT GRID TEMPLATE NOT SUPPORTED.", 2)
 endif

 kgds = 0
 call gdt_to_gds(gfld%igdtnum, gfld%igdtmpl, gfld%igdtlen, kgds, i_input, j_input, res)

 ip1_input = i_input + 1
 jp1_input = j_input + 1

 allocate(rlat(i_input,j_input))
 allocate(rlon(i_input,j_input))
 allocate(rlat_diff(i_input,j_input))
 allocate(rlon_diff(i_input,j_input))
 allocate(xpts(i_input,j_input))
 allocate(ypts(i_input,j_input))
 allocate(rlat_corner(ip1_input,jp1_input))
 allocate(rlon_corner(ip1_input,jp1_input))
 allocate(xpts_corner(ip1_input,jp1_input))
 allocate(ypts_corner(ip1_input,jp1_input))

 do j = 1, j_input
 do i = 1, i_input
   xpts(i,j) = float(i)
   ypts(i,j) = float(j)
 enddo
 enddo

 print*,"- COMPUTE GRID CELL CENTER COORDINATES."
 call gdswzd(kgds,1,(i_input*j_input),-9999.,xpts,ypts,rlon,rlat,nret)

 if (nret /= (i_input*j_input)) then
   call error_handler("GDSWZD RETURNED WRONG NUMBER OF POINTS.", 2)
 endif

 deallocate(xpts, ypts)

 do j = 1, jp1_input
 do i = 1, ip1_input
   xpts_corner(i,j) = float(i) - 0.5
   ypts_corner(i,j) = float(j) - 0.5
 enddo
 enddo

 print*,"- COMPUTE GRID CELL CORNER COORDINATES."
 call gdswzd(kgds,1,(ip1_input*jp1_input),-9999.,xpts_corner,ypts_corner,rlon_corner,rlat_corner,nret)

 if (nret /= (ip1_input*jp1_input)) then
   call error_handler("GDSWZD RETURNED WRONG NUMBER OF POINTS.", 2)
 endif

 deallocate(xpts_corner, ypts_corner)

 if (gfld%igdtnum == 0) then ! gfs lat/lon data

   print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."

   polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

   input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_input,j_input/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 else

   print*,"- CALL GridCreateNoPeriDim FOR INPUT GRID."

   input_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_input,j_input/), & 
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  rc=rc)

 endif

 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN GridCreate1PeriDim", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID LATITUDE."
 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID LONGITUDE."
 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))

 latitude = rlat
 longitude = rlon

 deallocate (rlat, rlon)

 print*,"- CALL FieldScatter FOR INPUT GRID LONGITUDE."
 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL FieldScatter FOR INPUT GRID LATITUDE."
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_src_ptr(i,j) = longitude(i,j)
     if (lon_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_src_ptr(i,j) = lon_src_ptr(i,j) - 360.0_esmf_kind_r8
     lat_src_ptr(i,j) = latitude(i,j)
   enddo
 enddo

 deallocate(latitude, longitude)

 allocate(latitude_corner(ip1_input,jp1_input))
 allocate(longitude_corner(ip1_input,jp1_input))

 latitude_corner = rlat_corner
 longitude_corner = rlon_corner

 deallocate (rlat_corner, rlon_corner)

 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(input_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_corner_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_corner_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_corner_src_ptr)
 call ESMF_GridGetCoord(input_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_corner_src_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     lon_corner_src_ptr(i,j) = longitude_corner(i,j)
     if (lon_corner_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_corner_src_ptr(i,j) = lon_corner_src_ptr(i,j) - 360.0_esmf_kind_r8
     lat_corner_src_ptr(i,j) = latitude_corner(i,j)
   enddo
 enddo

 deallocate(latitude_corner, longitude_corner)

 end subroutine define_input_grid_grib2

!> Setup the esmf grid object for the target grid.
!!
!! @param [in] localpet ESMF local persistent execution thread 
!! @param [in] npets Number of persistent execution threads
!! @author George Gayno NCEP/EMC   
 subroutine define_target_grid(localpet, npets)

 use netcdf
 use program_setup, only       : mosaic_file_target_grid, &
                                 orog_dir_target_grid,    &
                                 orog_files_target_grid,  &
                                 nsoill_out

 implicit none

 integer, intent(in)                   :: localpet, npets

 character(len=500)                    :: the_file

 integer                               :: error, ncid, extra
 integer                               :: id_tiles
 integer                               :: id_dim, id_grid_tiles
 integer                               :: tile
 integer, allocatable                  :: decomptile(:,:)
 integer(esmf_kind_i8), allocatable    :: landmask_one_tile(:,:)
 integer(esmf_kind_i8), allocatable    :: seamask_one_tile(:,:)
 
 real(esmf_kind_r8), allocatable       :: land_frac_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: latitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_s_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: longitude_w_one_tile(:,:)
 real(esmf_kind_r8), allocatable       :: terrain_one_tile(:,:)

 lsoil_target = nsoill_out
 
 print*,'- OPEN TARGET GRID MOSAIC FILE: ',trim(mosaic_file_target_grid)
 error=nf90_open(trim(mosaic_file_target_grid),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening grid mosaic file')

 print*,"- READ NUMBER OF TILES"
 error=nf90_inq_dimid(ncid, 'ntiles', id_tiles)
 call netcdf_err(error, 'reading ntile id')
 error=nf90_inquire_dimension(ncid,id_tiles,len=num_tiles_target_grid)
 call netcdf_err(error, 'reading ntiles')
 error=nf90_inq_varid(ncid, 'gridtiles', id_grid_tiles)
 call netcdf_err(error, 'reading gridtiles id')
 allocate(tiles_target_grid(num_tiles_target_grid))
 tiles_target_grid="NULL"
 print*,"- READ TILE NAMES"
 error=nf90_get_var(ncid, id_grid_tiles, tiles_target_grid)
 call netcdf_err(error, 'reading gridtiles')

 error = nf90_close(ncid)

 print*,'- NUMBER OF TILES, TARGET MODEL GRID IS ', num_tiles_target_grid

 if (mod(npets,num_tiles_target_grid) /= 0) then
   call error_handler("MUST RUN WITH TASK COUNT THAT IS A MULTIPLE OF # OF TILES.", 1)
 endif

!-----------------------------------------------------------------------
! Get the model grid specs and land mask from the orography files.
!-----------------------------------------------------------------------

 the_file = trim(orog_dir_target_grid) // trim(orog_files_target_grid(1))

 print*,'- OPEN FIRST TARGET GRID OROGRAPHY FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening orography file')
 print*,"- READ GRID DIMENSIONS"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, 'reading lon id')
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 call netcdf_err(error, 'reading lon')
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, 'reading lat id')
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 call netcdf_err(error, 'reading lat')
 error = nf90_close(ncid)

 print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 extra = npets / num_tiles_target_grid

 allocate(decomptile(2,num_tiles_target_grid))

 do tile = 1, num_tiles_target_grid
   decomptile(:,tile)=(/1,extra/)
 enddo

 print*,"- CALL GridCreateMosaic FOR TARGET GRID"
 target_grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file_target_grid), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(orog_dir_target_grid), rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", error)

!-----------------------------------------------------------------------
! Set target model landmask (1 - land, 0 - not land) and 
! seamask (1 - non-land, 0 -land).  Read lat/lon on target grid.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR TARGET GRID LAND FRACTION."
 land_frac_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_landmask", rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LANDMASK."
 landmask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_landmask", rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID SEAMASK."
 seamask_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_I8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_seamask", rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE_S."
 latitude_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_latitude_s", rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE_W."
 latitude_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_latitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE_S."
 longitude_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_longitude_s", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE_W."
 longitude_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_longitude_w", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID TERRAIN."
 terrain_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_terrain", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet == 0) then
   allocate(land_frac_one_tile(i_target,j_target))
   allocate(landmask_one_tile(i_target,j_target))
   allocate(seamask_one_tile(i_target,j_target))
   allocate(latitude_one_tile(i_target,j_target))
   allocate(latitude_s_one_tile(i_target,jp1_target))
   allocate(latitude_w_one_tile(ip1_target,j_target))
   allocate(longitude_one_tile(i_target,j_target))
   allocate(longitude_s_one_tile(i_target,jp1_target))
   allocate(longitude_w_one_tile(ip1_target,j_target))
   allocate(terrain_one_tile(i_target,j_target))
 else
   allocate(land_frac_one_tile(0,0))
   allocate(landmask_one_tile(0,0))
   allocate(seamask_one_tile(0,0))
   allocate(longitude_one_tile(0,0))
   allocate(longitude_s_one_tile(0,0))
   allocate(longitude_w_one_tile(0,0))
   allocate(latitude_one_tile(0,0))
   allocate(latitude_s_one_tile(0,0))
   allocate(latitude_w_one_tile(0,0))
   allocate(terrain_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid
   if (localpet == 0) then
     the_file = trim(orog_dir_target_grid) // trim(orog_files_target_grid(tile))
     call get_model_mask_terrain(trim(the_file), i_target, j_target, landmask_one_tile, &
                                 terrain_one_tile, land_frac_one_tile)
     
     seamask_one_tile = 0  ! all land
     where(floor(land_frac_one_tile) == 0) seamask_one_tile = 1  ! at least some non-land.
     landmask_one_tile = 0 ! all non-land
     where(ceiling(land_frac_one_tile) == 1) landmask_one_tile = 1  ! at least some land

     call get_model_latlons(mosaic_file_target_grid, orog_dir_target_grid, num_tiles_target_grid, tile, &
                            i_target, j_target, ip1_target, jp1_target, latitude_one_tile, &
                            latitude_s_one_tile, latitude_w_one_tile, longitude_one_tile, &
                            longitude_s_one_tile, longitude_w_one_tile)
   endif
   print*,"- CALL FieldScatter FOR TARGET GRID LAND FRACTION. TILE IS: ", tile
   call ESMF_FieldScatter(land_frac_target_grid, land_frac_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LANDMASK. TILE IS: ", tile
   call ESMF_FieldScatter(landmask_target_grid, landmask_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID SEAMASK. TILE IS: ", tile
   call ESMF_FieldScatter(seamask_target_grid, seamask_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_target_grid, longitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_s_target_grid, longitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(longitude_w_target_grid, longitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_target_grid, latitude_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE_S. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_s_target_grid, latitude_s_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE_W. TILE IS: ", tile
   call ESMF_FieldScatter(latitude_w_target_grid, latitude_w_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
   print*,"- CALL FieldScatter FOR TARGET GRID TERRAIN. TILE IS: ", tile
   call ESMF_FieldScatter(terrain_target_grid, terrain_one_tile, rootpet=0, tile=tile, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", error)
 enddo

 deallocate(land_frac_one_tile)
 deallocate(landmask_one_tile)
 deallocate(seamask_one_tile)
 deallocate(longitude_one_tile)
 deallocate(longitude_s_one_tile)
 deallocate(longitude_w_one_tile)
 deallocate(latitude_one_tile)
 deallocate(latitude_s_one_tile)
 deallocate(latitude_w_one_tile)
 deallocate(terrain_one_tile)

 end subroutine define_target_grid

!> Read model lat/lons for a single tile from the "grid" 
!! specificaton file.
!!
!! @param [in] mosaic_file The mosaic file associated with the 'grid' files.
!! @param [in] orog_dir  Directory containing the 'grid' and orography files.
!! @param [in] num_tiles  Total number of tiles
!! @param [in] tile  Tile number to be read
!! @param [in] i_tile "i" dimension of the tile
!! @param [in] j_tile "j" dimension of the tile
!! @param [in] ip1_tile "i" dimension of the tile plus 1
!! @param [in] jp1_tile "j" dimension of the tile plus 1
!! @param [out] latitude  grid box center latitude
!! @param [out] latitude_s  latitude of 'south' edge of grid box
!! @param [out] latitude_w  latitude of 'west' edge of grid box
!! @param [out] longitude  grid box center longitude
!! @param [out] longitude_s  longitude of 'south' edge of grid box
!! @param [out] longitude_w  longitude of 'west' edge of grid box
!! @author George Gayno NCEP/EMC   
 subroutine get_model_latlons(mosaic_file, orog_dir, num_tiles, tile, &
                              i_tile, j_tile, ip1_tile, jp1_tile,  &
                              latitude, latitude_s, latitude_w, &
                              longitude, longitude_s, longitude_w)

 use netcdf

 implicit none

 character(len=*), intent(in)      :: mosaic_file, orog_dir

 integer, intent(in)               :: num_tiles, tile
 integer, intent(in)               :: i_tile, j_tile
 integer, intent(in)               :: ip1_tile, jp1_tile

 real(esmf_kind_r8), intent(out)   :: latitude(i_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: latitude_s(i_tile, jp1_tile)
 real(esmf_kind_r8), intent(out)   :: latitude_w(ip1_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: longitude(i_tile, j_tile)
 real(esmf_kind_r8), intent(out)   :: longitude_s(i_tile, jp1_tile)
 real(esmf_kind_r8), intent(out)   :: longitude_w(ip1_tile, j_tile)

 character(len=50)                 :: grid_files(num_tiles)
 character(len=255)                :: grid_file

 integer                           :: error, id_var, ncid
 integer                           :: id_dim, nxp, nyp, i, j, ii, jj

 real(esmf_kind_r8), allocatable   :: tmpvar(:,:)

 print*,"- READ MODEL GRID FILE"

 print*,'- OPEN MOSAIC FILE: ', trim(mosaic_file)
 error=nf90_open(trim(mosaic_file), nf90_nowrite, ncid)
 call netcdf_err(error, 'opening mosaic file')

 print*,"- READ GRID FILE NAMES"
 error=nf90_inq_varid(ncid, 'gridfiles', id_var)
 call netcdf_err(error, 'reading gridfiles id')
 error=nf90_get_var(ncid, id_var, grid_files)
 call netcdf_err(error, 'reading gridfiles')

 error = nf90_close(ncid)

 grid_file = trim(orog_dir) // trim(grid_files(tile))

 print*,'- OPEN GRID FILE: ', trim(grid_file)
 error=nf90_open(trim(grid_file), nf90_nowrite, ncid)
 call netcdf_err(error, 'opening grid file')

 print*,'- READ NXP ID'
 error=nf90_inq_dimid(ncid, 'nxp', id_dim)
 call netcdf_err(error, 'reading nxp id')

 print*,'- READ NXP'
 error=nf90_inquire_dimension(ncid,id_dim,len=nxp)
 call netcdf_err(error, 'reading nxp')

 print*,'- READ NYP ID'
 error=nf90_inq_dimid(ncid, 'nyp', id_dim)
 call netcdf_err(error, 'reading nyp id')

 print*,'- READ NYP'
 error=nf90_inquire_dimension(ncid,id_dim,len=nyp)
 call netcdf_err(error, 'reading nyp')

 if ((nxp/2 /= i_tile) .or. (nyp/2 /= j_tile)) then
   call error_handler("DIMENSION MISMATCH IN GRID FILE.", 1)
 endif

 allocate(tmpvar(nxp,nyp))

 print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'x', id_var)
 call netcdf_err(error, 'reading longitude id')

 print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, tmpvar)
 call netcdf_err(error, 'reading longitude')

 do j = 1, j_tile
 do i = 1, i_tile
   ii = 2*i
   jj = 2*j
   longitude(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, jp1_tile
 do i = 1, i_tile
   ii = 2*i
   jj = (2*j) - 1
   longitude_s(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, j_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = 2*j
   longitude_w(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'y', id_var)
 call netcdf_err(error, 'reading latitude id')

 print*,'- READ LATIITUDE'
 error=nf90_get_var(ncid, id_var, tmpvar)
 call netcdf_err(error, 'reading latitude')

 do j = 1, j_tile
 do i = 1, i_tile
   ii = 2*i
   jj = 2*j
   latitude(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, jp1_tile
 do i = 1, i_tile
   ii = 2*i
   jj = (2*j) - 1
   latitude_s(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 do j = 1, j_tile
 do i = 1, ip1_tile
   ii = (2*i) - 1
   jj = 2*j
   latitude_w(i,j) = tmpvar(ii,jj)
 enddo
 enddo

 deallocate(tmpvar)

 error = nf90_close(ncid)

 end subroutine get_model_latlons
 
!> Read the model land mask and terrain for a single tile
!! from the orography file.
!!
!! @param [in] orog_file  Path/name of orography file
!! @param [in] idim  "i" dimension of tile
!! @param [in] jdim  "j" dimension of tile
!! @param [out] mask  land mask of tile
!! @param [out] terrain  terrain height of tile
!! @param [out] land_frac  The fraction of the grid point that is land.
!! @author George Gayno NCEP/EMC   
 subroutine get_model_mask_terrain(orog_file, idim, jdim, mask, terrain, land_frac)

 use netcdf

 implicit none

 character(len=*), intent(in)       :: orog_file

 integer, intent(in)                :: idim, jdim
 integer(esmf_kind_i8), intent(out) :: mask(idim,jdim)

 real(esmf_kind_i8), intent(out)    :: terrain(idim,jdim)
 real(esmf_kind_i8), intent(out)    :: land_frac(idim,jdim)

 integer                            :: error, lat, lon
 integer                            :: ncid, id_dim, id_var

 real(kind=4), allocatable          :: dummy(:,:)

 print*,"- READ MODEL LAND MASK FILE"

 print*,'- OPEN LAND MASK FILE: ', orog_file
 error=nf90_open(orog_file,nf90_nowrite,ncid)
 call netcdf_err(error, 'opening land mask file')

 print*,"- READ I-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, 'reading idim id')
 error=nf90_inquire_dimension(ncid,id_dim,len=lon)
 call netcdf_err(error, 'reading idim')

 print*,"- READ J-DIMENSION"
 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, 'reading jdim id')
 error=nf90_inquire_dimension(ncid,id_dim,len=lat)
 call netcdf_err(error, 'reading jdim')

 print*,"- I/J DIMENSIONS: ", lon, lat

 if ((lon /= idim) .or. (lat /= jdim)) then
   call error_handler("MISMATCH IN DIMENSIONS.", 1)
 endif

 allocate(dummy(idim,jdim))

 print*,"- READ LAND MASK"
 error=nf90_inq_varid(ncid, 'slmsk', id_var)
 call netcdf_err(error, 'reading slmsk id')
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, 'reading slmsk')
 mask = nint(dummy)

 print*,"- READ RAW OROGRAPHY."
 error=nf90_inq_varid(ncid, 'orog_raw', id_var)
 call netcdf_err(error, 'reading orog_raw id')
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, 'reading orog_raw')
 terrain = dummy

 print*,"- READ LAND FRACTION."
 error=nf90_inq_varid(ncid, 'land_frac', id_var)
 call netcdf_err(error, 'reading land_frac id')
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, 'reading orog_raw')
 land_frac = dummy

!print*,'land frac ',dummy(idim/2,:)
 error = nf90_close(ncid)

 deallocate (dummy)

 end subroutine get_model_mask_terrain

!> Deallocate all esmf grid objects.
!!
!! @author George Gayno NCEP/EMC   
 subroutine cleanup_input_target_grid_data

 implicit none

 integer                                :: rc

 print*,"- DESTROY MODEL DATA."
 
 call ESMF_FieldDestroy(latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(longitude_input_grid,rc=rc)
 if (ESMF_FieldIsCreated(latitude_s_input_grid)) then
   call ESMF_FieldDestroy(latitude_s_input_grid, rc=rc)
 endif
 if (ESMF_FieldIsCreated(latitude_w_input_grid)) then
   call ESMF_FieldDestroy(latitude_w_input_grid, rc=rc)
 endif
 if (ESMF_FieldIsCreated(longitude_s_input_grid)) then
   call ESMF_FieldDestroy(longitude_s_input_grid, rc=rc)
 endif
 if (ESMF_FieldIsCreated(longitude_w_input_grid)) then
   call ESMF_FieldDestroy(longitude_w_input_grid, rc=rc)
 endif
 call ESMF_FieldDestroy(landmask_target_grid, rc=rc)
 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(latitude_s_target_grid)) then
   call ESMF_FieldDestroy(latitude_s_target_grid, rc=rc)
 endif
 if (ESMF_FieldIsCreated(latitude_w_target_grid)) then
   call ESMF_FieldDestroy(latitude_w_target_grid, rc=rc)
 endif
 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(longitude_s_target_grid)) then
   call ESMF_FieldDestroy(longitude_s_target_grid, rc=rc)
 endif
 if (ESMF_FieldIsCreated(longitude_w_target_grid)) then
   call ESMF_FieldDestroy(longitude_w_target_grid, rc=rc)
 endif
 call ESMF_FieldDestroy(seamask_target_grid, rc=rc)
 call ESMF_FieldDestroy(terrain_target_grid, rc=rc)
 call ESMF_GridDestroy(input_grid, rc=rc)
 call ESMF_GridDestroy(target_grid, rc=rc)

 end subroutine cleanup_input_target_grid_data

!> Convert the GRIB2 grid description template to
!! to the GRIB1 grid description section.
!!
!! @param [in] igdtnum GRIB2 grid description template number.
!! @param [in] igdstmpl Length of grib2 grid description template.
!! @param [in] igdtlen Array of GRIB2 grid description template octets.
!! @param [out] kgds Array of GRIB1 grid description octets.
!! @param [out] ni I-dimension of grid.
!! @param [out] nj J-dimension of grid.
!! @param [out] res Resolution of grid in km.
!! @author George Gayno NCEP/EMC   
 subroutine gdt_to_gds(igdtnum, igdstmpl, igdtlen, kgds, ni, nj, res)

 implicit none

 integer, intent(in   )  :: igdtnum, igdtlen, igdstmpl(igdtlen)
 integer, intent(  out)  :: kgds(200), ni, nj
 integer                 :: iscale

 real,    intent(  out)  :: res

 kgds=0

 if (igdtnum.eq.32769) then        ! rot lat/lon b grid

     iscale=igdstmpl(10)*igdstmpl(11)
     if (iscale == 0) iscale = 1e6
     kgds(1)=205                    ! oct 6,     rotated lat/lon for Non-E
                                    !            Stagger grid
     kgds(2)=igdstmpl(8)            ! octs 7-8,  Ni
     ni = kgds(2)
     kgds(3)=igdstmpl(9)            ! octs 9-10, Nj
     nj = kgds(3)
     kgds(4)=nint(float(igdstmpl(12))/float(iscale)*1000.)  ! octs 11-13, Lat of
                                                            ! 1st grid point
     kgds(5)=nint(float(igdstmpl(13))/float(iscale)*1000.)  ! octs 14-16, Lon of
                                                            ! 1st grid point

     kgds(6)=0                      ! oct 17, resolution and component flags
     if (igdstmpl(1)==2 ) kgds(6)=64
     if ( btest(igdstmpl(14),4).OR.btest(igdstmpl(14),5) )  kgds(6)=kgds(6)+128
     if ( btest(igdstmpl(14),3) ) kgds(6)=kgds(6)+8

     kgds(7)=nint(float(igdstmpl(15))/float(iscale)*1000.)       ! octs 18-20,
                                                                 ! Lat of cent of rotation
     kgds(8)=nint(float(igdstmpl(16))/float(iscale)*1000.)       ! octs 21-23,
                                                                 ! Lon of cent of rotation
     kgds(9)=nint(float(igdstmpl(17))/float(iscale)*1000.)       ! octs 24-25,
                                                                 ! Di
     kgds(10)=nint(float(igdstmpl(18))/float(iscale)*1000.)      ! octs 26-27,
                                                                 ! Dj

     kgds(11) = 0                   ! oct 28, scan mode
     if (btest(igdstmpl(19),7)) kgds(11) = 128
     if (btest(igdstmpl(19),6)) kgds(11) = kgds(11) +  64
     if (btest(igdstmpl(19),5)) kgds(11) = kgds(11) +  32

     kgds(12)=nint(float(igdstmpl(20))/float(iscale)*1000.) ! octs 29-31, Lat of
                                                            ! last grid point
     kgds(13)=nint(float(igdstmpl(21))/float(iscale)*1000.) ! octs 32-34, Lon of
                                                            ! last grid point

     kgds(19)=0    ! oct 4, # vert coordinate parameters
     kgds(20)=255  ! oct 5, used for thinned grids, set to 255

     res = ((float(kgds(9)) / 1000.0) + (float(kgds(10)) / 1000.0)) &
             * 0.5 * 111.0

   elseif(igdtnum==30) then

     kgds(1)=3                      ! oct 6,     lambert conformal
     kgds(2)=igdstmpl(8)            ! octs 7-8,  Ni
     ni = kgds(2)
     kgds(3)=igdstmpl(9)            ! octs 9-10, Nj
     nj = kgds(3)

     iscale = 1e6
     kgds(4) = nint(float(igdstmpl(10))/1000.0)
     kgds(5) = nint(float(igdstmpl(11))/1000.0)

     kgds(6)=0                      ! oct 17, resolution and component flags
     if (igdstmpl(1)==2 ) kgds(6)=64
     if ( btest(igdstmpl(12),4).OR.btest(igdstmpl(12),5) )  kgds(6)=kgds(6)+128
     if ( btest(igdstmpl(12),3) ) kgds(6)=kgds(6)+8

     kgds(7) = nint(float(igdstmpl(14))/1000.0)
     kgds(8) = nint(float(igdstmpl(15))/1000.0)
     kgds(9) = nint(float(igdstmpl(16))/1000.0)
     kgds(10) = 0 

     kgds(11) = 0                   ! oct 28, scan mode
     if (btest(igdstmpl(18),7)) kgds(11) = 128
     if (btest(igdstmpl(18),6)) kgds(11) = kgds(11) +  64
     if (btest(igdstmpl(18),5)) kgds(11) = kgds(11) +  32

     kgds(12) = nint(float(igdstmpl(19))/1000.0)
     kgds(13) = nint(float(igdstmpl(20))/1000.0)
     kgds(14) = -90
     kgds(15) = 0

   elseif(igdtnum==0) then  ! lat/lon grid

     iscale=igdstmpl(10)*igdstmpl(11)
     if (iscale == 0) iscale = 1e6
     kgds(1)=0                   ! oct 6, data representation type.
     kgds(2)=igdstmpl(8)         ! octs 7-8, Ni
     ni = kgds(2)
     kgds(3)=igdstmpl(9)         ! octs 9-10, Nj
     nj = kgds(3)
     kgds(4)=nint(float(igdstmpl(12))/float(iscale)*1000.)  ! octs 11-13, Lat of 1st grid point
     kgds(5)=nint(float(igdstmpl(13))/float(iscale)*1000.)  ! octs 14-16, Lon of 1st grid point

     kgds(6)=0                   ! oct 17, resolution and component flags
     if (igdstmpl(1)==2 ) kgds(6)=64
     if ( btest(igdstmpl(14),4).OR.btest(igdstmpl(14),5) ) kgds(6)=kgds(6)+128
     if ( btest(igdstmpl(14),3) ) kgds(6)=kgds(6)+8

     kgds(7)=nint(float(igdstmpl(15))/float(iscale)*1000.)  ! octs 18-20, Lat of last grid point
     kgds(8)=nint(float(igdstmpl(16))/float(iscale)*1000.)  ! octs 21-23, Lon of last grid point
     kgds(9)=nint(float(igdstmpl(17))/float(iscale)*1000.)  ! octs 24-25, "i" resolution.
     kgds(10)=nint(float(igdstmpl(18))/float(iscale)*1000.) ! octs 26-27, "j" resolution.

     kgds(11) = 0              ! oct 28, scan mode
     if (btest(igdstmpl(19),7)) kgds(11) = 128
     if (btest(igdstmpl(19),6)) kgds(11) = kgds(11) +  64
     if (btest(igdstmpl(19),5)) kgds(11) = kgds(11) +  32

     kgds(12)=0      ! octs 29-32, reserved
     kgds(19)=0      ! oct 4, # vert coordinate parameters
     kgds(20)=255    ! oct 5, used for thinned grids, set to 255

   else
      
     call error_handler("UNRECOGNIZED INPUT GRID TYPE ", 1)

 endif

 end subroutine gdt_to_gds

 end module model_grid
