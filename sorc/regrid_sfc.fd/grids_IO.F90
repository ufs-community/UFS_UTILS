!> @file
!! @brief ESMF grid-specific routines for GFS regridding program, including IO.
!! @author Clara Draper, Aug 2024.
 module grids_IO

 use esmf
 use netcdf
 use ESMF_LogPublicMod
 use utilities, only     : error_handler, netcdf_err

 implicit none

 private

 integer, public, parameter  :: n_tiles=6 !< number tiles in fv3 grid
 ! mask values for land / ocean mask built from veg type
 integer, public, parameter  :: vtype_water=0, & !< non-land
                                vtype_landice=15 !< land ice
 ! mask values for soilsnow_mask calculated in the GSI EnKF
 integer, public, parameter  :: mtype_water=0, & !< water
                                mtype_snow=2     !< snow
 type, public  :: grid_setup_type !< kl
        character(7)   :: descriptor       !< options: gau_inc fv3_rst 
        character(100) :: fname !< n
        character(100) :: dir !< l 
        character(15)  :: mask_variable(1) !< j l
        character(100) :: fname_mask !< jlk
        character(100) :: dir_mask !< j
        character(100) :: fname_coord !< j 
        character(100) :: dir_coord !< l 
        integer        :: ires !< j 
        integer        :: jres !< l 
 end type

 public :: setup_grid, &
           read_into_fields, &
           write_from_fields

 contains

!> Create ESMF grid objects, with mask if requested

 subroutine setup_grid(localpet, npets, grid_setup, mod_grid )

 implicit none

 ! INTENT IN
 type(grid_setup_type), intent(in)    :: grid_setup
 integer, intent(in)            :: localpet, npets

 ! INTENT OUT
 type(esmf_grid)                :: mod_grid

 ! LOCAL
 type(esmf_field)               :: mask_field(1,1)
 real(esmf_kind_r8), pointer    :: ptr_maskvar(:,:)
 integer(esmf_kind_i4), pointer :: ptr_mask(:,:)

 integer                        :: ierr, ncid, tile

!--------------------------
! Create grid object, and set up pet distribution

 select case (grid_setup%descriptor)
 case ('fv3_rst')
     call create_grid_fv3(grid_setup%ires, trim(grid_setup%dir_coord), npets, localpet ,mod_grid)
 case ('gau_inc')
     call create_grid_gauss(grid_setup, npets, localpet,  mod_grid)
 case default
     call error_handler("unknown grid_setup%descriptor in setup_grid", 1)
 end select

!--------------------------
! Calculate and add the mask

 mask_field(1,1) = ESMF_FieldCreate(mod_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input variable for mask", &
                                   rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate, mask_variable", ierr)

 call read_into_fields(localpet, grid_setup%ires, grid_setup%jres, trim(grid_setup%fname_mask), &
                         trim(grid_setup%dir_mask), grid_setup, 1, &
                         grid_setup%mask_variable(1), mask_field(1,1))

! get pointer to mask
 call ESMF_FieldGet(mask_field(1,1), &
                    farrayPtr=ptr_maskvar, &
                    rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", ierr)

! create and get pointer to the mask
 call ESMF_GridAddItem(mod_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("in GridAddItem mask", ierr)

 call ESMF_GridGetItem(mod_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=ptr_mask, &
                       rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("in GridGetItem mask", ierr)

! calculate the mask
 ptr_mask = 1 ! initialize land everywhere
 select case (trim(grid_setup%mask_variable(1)))
 case("vegetation_type") ! removing non-land and glaciers using veg class
     where (nint(ptr_maskvar) == vtype_water )   ptr_mask = 0 ! exclude water
     where (nint(ptr_maskvar) == vtype_landice ) ptr_mask = 0 ! exclude glaciers
 case("soilsnow_mask") ! removing snow and non-land using pre-computed mask
     where (nint(ptr_maskvar) == mtype_water )   ptr_mask = 0 ! exclude non-soil
     where (nint(ptr_maskvar) == mtype_snow ) ptr_mask = 0 ! exclude snow
 case default
    call error_handler("unknown mask_variable", 1)
 end select

! destroy mask field
 call ESMF_FieldDestroy(mask_field(1,1),rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING FIELD", ierr)

 end subroutine setup_grid

!> read variables from fv3 netcdf restart file into ESMF Fields
 subroutine read_into_fields(localpet, i_dim, j_dim , fname_read, dir_read, &
                               grid_setup, n_vars, variable_list, fields)

 implicit none

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read
 type(grid_setup_type), intent(in)        :: grid_setup

 character(len=15), dimension(n_vars), intent(in)   :: variable_list

 ! INTENT OUT
 type(esmf_field), dimension(1,n_vars), intent(inout) :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, v
 integer                         :: n_files
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_in(:,:,:)

 if (localpet==0) then
     allocate(array_in(n_vars,i_dim, j_dim))
     allocate(array2D(i_dim, j_dim))
 else
     allocate(array_in(0,0,0))
     allocate(array2D(0,0))
 end if

 select case (grid_setup%descriptor)
 case ('fv3_rst')
     n_files=n_tiles
 case ('gau_inc')
     n_files=1
 case default
     call error_handler("unknown grid_setup%descriptor in read into fields", 1)
 end select

 do tt = 1, n_files

      ! read from restart
      if (localpet == 0) then

         if ( n_files > 1) then
             write(tchar,'(i1)') tt
             fname = dir_read//"/"//fname_read//".tile"//tchar//".nc"
         else
             fname = dir_read//"/"//fname_read
         endif

         print *, 'Reading ', trim(fname)

         ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
         call netcdf_err(ierr, 'opening: '//trim(fname) )

         do v =1, n_vars
             print *, 'Reading ', trim(variable_list(v))
             ierr=nf90_inq_varid(ncid, trim(variable_list(v)), id_var)
             call netcdf_err(ierr, 'reading variable id' )

             ierr=nf90_get_var(ncid, id_var, array_in(v,:,:))
             call netcdf_err(ierr, 'reading variable' )
         enddo
         ierr = nf90_close(ncid)

      endif

      ! scatter
      do v =1, n_vars
          array2D=array_in(v,:,:) ! scatter misbehaves if given indexed 3D array.
          call ESMF_FieldScatter(fields(1,v), array2D, rootpet=0, tile=tt, rc=ierr)
          if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldScatter", ierr)

      enddo

 enddo

 ! clean up
 deallocate(array_in)

 end subroutine read_into_fields


!> read lat and lon from SCRIP file, for use in Gaussian grid
 subroutine lonlat_read_into_fields(localpet, i_dim, j_dim, fname_read, dir_read, &
                               grid_setup, gauss_grid, lon_fields, lat_fields)

 implicit none

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read
 type(grid_setup_type), intent(in)        :: grid_setup
 type(esmf_grid), intent(in)     :: gauss_grid


 ! INTENT OUT
 type(esmf_field), dimension(2), intent(out) :: lon_fields
 type(esmf_field), dimension(2), intent(out) :: lat_fields

 ! LOCAL
 integer                         :: id_var, ncid, ierr
 integer                         :: i,j,v
 character(len=500)              :: fname
 real(esmf_kind_r8)              :: vec1(i_dim*j_dim), vec4(4,i_dim*j_dim)
 real(esmf_kind_r8), allocatable :: array_lat(:,:,:), array_lon(:,:,:), array2d(:,:)
 real(esmf_kind_r8)              :: array_tmp(i_dim,j_dim)
 character(len=30)               :: vname
 character(len=15), dimension(2) :: lonvar_list
 character(len=15), dimension(2) :: latvar_list

 ! center coord names
 lonvar_list(1) = 'grid_center_lon'
 latvar_list(1) = 'grid_center_lat'

 ! corner coord names
 lonvar_list(2) = 'grid_corner_lon'
 latvar_list(2) = 'grid_corner_lat'

 ! Create the fields
 do v=1,2
        vname="gauss_grid_" // trim(lonvar_list(v))
        lon_fields(v) = ESMF_FieldCreate(gauss_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name=vname, rc=ierr)
        if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate", ierr)

        vname="gauss_grid_" // trim(latvar_list(v))
        lat_fields(v) = ESMF_FieldCreate(gauss_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name=vname, rc=ierr)
        if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate", ierr)
 enddo

 if (localpet==0) then
     allocate(array_lat(2, i_dim, j_dim))
     allocate(array_lon(2, i_dim, j_dim))
     allocate(array2d( i_dim, j_dim))
 else
     allocate(array_lat(2,0,0))
     allocate(array_lon(2,0,0))
     allocate(array2d(0,0))
 endif

 ! read coordinates from restart
 if (localpet == 0) then

     fname = dir_read//"/"//fname_read

     print *, 'Reading ', trim(fname)

     ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
     call netcdf_err(ierr, 'opening: '//trim(fname) )

     ! reading the center points

     print *, 'Reading ', trim(lonvar_list(1))
     ierr=nf90_inq_varid(ncid, trim(lonvar_list(1)), id_var)
     call netcdf_err(ierr, 'reading variable id' )

     ierr=nf90_get_var(ncid, id_var, vec1)
     call netcdf_err(ierr, 'reading variable' )

     array_lon(1,:,:) = reshape(vec1,(/i_dim,j_dim/))

     print *, 'Reading ', trim(latvar_list(1))
     ierr=nf90_inq_varid(ncid, trim(latvar_list(1)), id_var)
     call netcdf_err(ierr, 'reading variable id' )

     ierr=nf90_get_var(ncid, id_var, vec1)
     call netcdf_err(ierr, 'reading variable' )


     ! increment files are S->N, need to flip input lats from N->S
     if ( grid_setup%descriptor == 'gau_inc') then
        do j =1,j_dim
                array_tmp = reshape(vec1,(/i_dim,j_dim/))
                array_lat(1,:,j) = array_tmp(:,j_dim+1-j)
        enddo
     else
        array_lat(1,:,:) = reshape(vec1,(/i_dim,j_dim/))
     endif

     ! read in the corners

     print *, 'Reading ', trim(lonvar_list(2))
     ierr=nf90_inq_varid(ncid, trim(lonvar_list(2)), id_var)
     call netcdf_err(ierr, 'reading variable id' )

     ierr=nf90_get_var(ncid, id_var, vec4)
     call netcdf_err(ierr, 'reading variable' )

     ! 2nd element of vec4 is the NW corner
     array_lon(2,:,:) = reshape(vec4(2,:),(/i_dim,j_dim/))

     print *, 'Reading ', trim(latvar_list(2))
     ierr=nf90_inq_varid(ncid, trim(latvar_list(2)), id_var)
     call netcdf_err(ierr, 'reading variable id' )

     ierr=nf90_get_var(ncid, id_var, vec4)
     call netcdf_err(ierr, 'reading variable' )

     ! increment files are S->N, need to flip input lats from N->S
     if ( grid_setup%descriptor == 'gau_inc') then
        do j =1,j_dim
                ! 2nd element of vec4 is the NW corner
                array_tmp = reshape(vec4(2,:),(/i_dim,j_dim/))
                array_lat(2,:,j) = array_tmp(:,j_dim+1-j)
        enddo
     else
        ! 2nd element of vec4 is the NW corner
        array_lat(2,:,:) = reshape(vec4(2,:),(/i_dim,j_dim/))
     endif

     ierr = nf90_close(ncid)
  endif

  ! scatter the arrays into the fields
  do v =1,2

      ! scatter longitudes
      if (localpet==0) print *, 'scattering lon', lonvar_list(v)
      array2d=array_lon(v,:,:)
      call ESMF_FieldScatter(lon_fields(v), array2d, rootpet=0,  rc=ierr)
      if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldScatter", ierr)

      ! scatter latitude
      if (localpet==0) print *, 'scattering lat', latvar_list(v)
      array2d=array_lat(v,:,:)
      call ESMF_FieldScatter(lat_fields(v), array2d, rootpet=0,  rc=ierr)
      if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldScatter", ierr)
 enddo

 ! clean up
 deallocate(array_lat)
 deallocate(array_lon)
 deallocate(array2d)

 end subroutine lonlat_read_into_fields

!> write variables from ESMF Fields into netcdf restart-like file

 subroutine write_from_fields(localpet, i_dim, j_dim , fname_out, dir_out, &
                                n_vars, n_tims, variable_list, fields)

 implicit none

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim,  n_vars, n_tims
 character(*), intent(in)        :: fname_out
 character(*), intent(in)        :: dir_out
 character(15), dimension(n_vars), intent(in)     :: variable_list
 type(esmf_field), dimension(n_tims,n_vars), intent(in)  :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, &
                                    id_x, id_y, id_t, v, t
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_out(:,:,:,:)

 do v = 1, n_vars
        if (localpet == 0)  print *, 'Writing ', trim(variable_list(v)), ' into field'
 enddo

 if (localpet==0) then
     allocate(array_out(n_vars, i_dim, j_dim, n_tims))
     allocate(array2D(i_dim, j_dim))
 else
     allocate(array_out(0,0,0,0))
     allocate(array2D(0,0))
 end if

 do tt = 1, n_tiles

      ! fetch data (all PETs)
      do t =1 , n_tims
          do v = 1, n_vars
              call ESMF_FieldGather(fields(t,v), array2D, rootPet=0, tile=tt, rc=ierr)
              if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                 call error_handler("IN FieldGather", ierr)
              array_out(v,:,:,t) = array2D
          enddo
      enddo

      ! write to netcdf
      if (localpet == 0) then

         ! open file, set dimensions
         write(tchar,'(i1)') tt
         fname = dir_out//"/"//fname_out//".tile"//tchar//".nc"

         ierr = nf90_create(trim(fname), NF90_NETCDF4, ncid)
         call netcdf_err(ierr, 'creating file='//trim(fname) )

         if (n_tims>1) then ! UFS_UTILS expects input with no time dim
                           ! GFS (for IAU) expects a time dimension
                           ! later: update GFS to not expect a time dimension
             ierr = nf90_def_dim(ncid, 'Time', n_tims, id_t)
             call netcdf_err(ierr, 'defining taxis dimension' )
         endif

         ierr = nf90_def_dim(ncid, 'xaxis_1', i_dim, id_x)
         call netcdf_err(ierr, 'defining xaxis dimension' )

         ierr = nf90_def_dim(ncid, 'yaxis_1', j_dim, id_y)
         call netcdf_err(ierr, 'defining yaxis dimension' )


         do v=1, n_vars

             if (n_tims>1) then
                 ! UFS model code to read in the increments is expecting
                 ! dimensions: time, y, x (in the ncdump read out - which reverses fortran indexes)
                 ! need dimensions to be x,y,t below.
                 ierr = nf90_def_var(ncid, trim(variable_list(v)), NF90_DOUBLE, &
                                     (/id_x, id_y, id_t/) , id_var)

                 call netcdf_err(ierr, 'defining '//variable_list(v) )
             else
                 ierr = nf90_def_var(ncid, trim(variable_list(v)), NF90_DOUBLE, &
                                     (/id_x, id_y/) , id_var)
             endif

             call netcdf_err(ierr, 'defining '//variable_list(v) )

             ierr = nf90_put_var( ncid, id_var, array_out(v,:,:,:) )
             call netcdf_err(ierr, 'writing '//variable_list(v) )

         enddo

         ierr = nf90_close(ncid)

      endif

 enddo

 ! clean up
 deallocate(array_out)

 end subroutine write_from_fields

!> subroutine to create grid object for fv3 grid
!!  also sets distribution across procs

 subroutine create_grid_fv3(res_atm, dir_fix, npets, localpet, fv3_grid)

! INTENT IN
 integer, intent(in)    :: npets, localpet
 integer, intent(in)    :: res_atm
 character(*), intent(in)       :: dir_fix

 ! INTENT OUT
 type(esmf_grid)                :: fv3_grid

 integer                :: ierr, extra, tile
 integer                :: decomptile(2,n_tiles)

 character(len=5)       :: rchar
 character(len=200)     :: fname

 if (localpet == 0) print*," creating fv3 grid for ", res_atm

! pet distribution
 extra = npets / n_tiles
 do tile = 1, n_tiles
   decomptile(:,tile)=(/1,extra/)
 enddo

 ! mosaic file
 write(rchar,'(i5)') res_atm
 fname = trim(dir_fix)//"/C"//trim(adjustl(rchar))// "_mosaic.nc"

! create the grid
 fv3_grid = ESMF_GridCreateMosaic(filename=trim(fname), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(dir_fix), &
                                  rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", ierr)

 end subroutine create_grid_fv3

!> subroutine to create grid object for gaussian grids
!!  also sets distribution across procs

 subroutine create_grid_gauss(grid_setup, npets, localpet, gauss_grid)

 ! INTENT IN
 type(grid_setup_type), intent(in) :: grid_setup
 integer, intent(in)   :: npets, localpet

 ! INTENT OUT
 type(esmf_grid)                   :: gauss_grid

 type(esmf_polekind_flag)          :: polekindflag(2)
 type(esmf_staggerloc)             :: staggerloc(2)

 type(esmf_field)                  :: lon_fields(2)
 type(esmf_field)                  :: lat_fields(2)
 real(esmf_kind_r8), pointer       :: lon_ptr_field(:,:), lon_ptr_coord(:,:)
 real(esmf_kind_r8), pointer       :: lat_ptr_field(:,:), lat_ptr_coord(:,:)

 integer :: ierr, i_dim, j_dim, v

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE


 i_dim = grid_setup%ires
 j_dim = grid_setup%jres

 if (localpet == 0) print*," creating gauss grid for ", i_dim, j_dim
 gauss_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_dim,j_dim/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreate1PeriDim", ierr)

 ! read lat lon coordinates into fields
 call lonlat_read_into_fields(localpet, i_dim, j_dim,  &
                         trim(grid_setup%fname_coord), trim(grid_setup%dir_coord), &
                         grid_setup, gauss_grid, lon_fields, lat_fields)

 staggerloc = (/ ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER /)
 do v = 1,2
     ! add coordinates to the grid
     call ESMF_GridAddCoord(gauss_grid, &
                            staggerloc=staggerloc(v), rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN GridAddCoord", ierr)

     nullify(lon_ptr_coord)
     nullify(lon_ptr_field)

     ! get pointer to lon coord
     if (localpet == 0) print*," getting GridCoord for long"
     call ESMF_GridGetCoord(gauss_grid, &
                            staggerLoc=staggerloc(v), &
                            coordDim=1, &
                            farrayPtr=lon_ptr_coord, rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN GridGetCoord longitude ", ierr)

     ! fetch lon from field into pointer
     call ESMF_FieldGet(lon_fields(v), &
                        farrayPtr=lon_ptr_field, &
                        rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGet longitude", ierr)

     ! set coord pointe to field pointer
     lon_ptr_coord = lon_ptr_field

     nullify(lat_ptr_coord)
     nullify(lat_ptr_field)

     ! get pointer to lat coord. Need bounds?
     call ESMF_GridGetCoord(gauss_grid, &
                            staggerLoc=staggerloc(v), &
                            coordDim=2, &
                            farrayPtr=lat_ptr_coord, rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN GridGetCoord latitude", ierr)

     ! fetch lat from field into pointer
     call ESMF_FieldGet(lat_fields(v), &
                        farrayPtr=lat_ptr_field, &
                        rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGet latitude", ierr)

     ! set lat coord to field values
     lat_ptr_coord =  lat_ptr_field
 enddo

 ! clean up
 do v = 1,2
     call ESMF_FieldDestroy(lat_fields(v),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING FIELD", ierr)
     call ESMF_FieldDestroy(lon_fields(v),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING FIELD", ierr)
 enddo

 end subroutine create_grid_gauss

end  module grids_IO
