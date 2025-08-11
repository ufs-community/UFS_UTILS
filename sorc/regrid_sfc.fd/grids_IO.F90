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
 integer, public, parameter  :: vtype_nonland=0, & !< non-land
                                vtype_water=17, & !< water
                                vtype_landice=15 !< land ice
 ! mask values for soilsnow_mask calculated in the GSI EnKF
 integer, public, parameter  :: mtype_water=0, & !< water
                                mtype_snow=2     !< snow
 type, public  :: grid_setup_type
        character(7)   :: descriptor       !< options: gau_inc fv3_rst 
        character(100) :: fname            !< file name
        character(100) :: dir              !< directory
        character(15)  :: mask_variable(1) !< name of variables used for mask
        character(100) :: fname_mask       !< file name for reading in mask
        character(100) :: dir_mask         !< directory name for reading in mask
        logical        :: mask_from_input  !< read mask from input file
        character(100) :: fname_coord      !< file name with coordinate info
        character(100) :: dir_coord        !< directory name for coordinate info
        integer        :: ires             !< latitudinal dimension
        integer        :: jres             !< longitudinal dimension
 end type

 public :: setup_grid, &
           read_into_fields, &
           write_from_fields

 contains

!> Create ESMF grid objects, with mask if requested
!! @param[in] localpet          local pet
!! @param[in] npets             total number of pets
!! @param[in] grid_setup        data structure with grid details 
!! @param[out] mod_grid         output esmf_grid structure 
!! @param[in] timestamp      timestep of input file

 subroutine setup_grid(localpet, npets, grid_setup, mod_grid, timestamp )

 implicit none

 ! INTENT IN
 type(grid_setup_type), intent(in)    :: grid_setup
 integer, intent(in)            :: localpet, npets
 integer, intent(in), optional  :: timestamp

 ! INTENT OUT
 type(esmf_grid), intent(out)   :: mod_grid


 ! LOCAL
 type(esmf_field)               :: mask_field(1,1)
 real(esmf_kind_r8), pointer    :: ptr_maskvar(:,:)
 integer(esmf_kind_i4), pointer :: ptr_mask(:,:)

 integer                        :: ierr, ncid, tile
 character(len=128)             :: fname_mask
 character(len=3)               :: tstr

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

 if (present(timestamp)) then 
    write(tstr,"(I3.3)") timestamp
    fname_mask = trim(grid_setup%fname_mask)//tstr//".nc"
 else
    fname_mask = trim(grid_setup%fname_mask)
 endif

 call read_into_fields(localpet, grid_setup%ires, grid_setup%jres, trim(fname_mask), &
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
 case("vegetation_type") ! removing non-land, water, and glaciers using veg class
     where (nint(ptr_maskvar) == vtype_nonland)   ptr_mask = 0 ! exclude non-land
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
!! @param[in] localpet          local pet
!! @param[in] i_dim             longitudinal dimension
!! @param[in] j_dim             latitudinal dimension
!! @param[in] fname_read        file name to read in
!! @param[in] dir_read          directory of file name to read
!! @param[in] grid_setup        grid details
!! @param[in] n_vars            number of variables to read in 
!! @param[in] variable_list     variables to read in
!! @param[inout] fields         fields to read variables into

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
 integer                         :: tt, id_var, ncid, ierr, v, j
 integer                         :: n_files
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_in(:,:,:)
 real(esmf_kind_r8), allocatable :: temp_array(:,:,:)

 allocate(array_in(n_vars,i_dim, j_dim))
 allocate(array2D(i_dim, j_dim))

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

         ! increment files are S->N, ESMF expects N->S
         if  ( grid_setup%descriptor == 'gau_inc') then 
            allocate(temp_array(n_vars,i_dim, j_dim))
            temp_array = array_in
            do j=1,j_dim
                array_in(:,:,j) = temp_array(:,:,j_dim-j+1)
            enddo
            deallocate(temp_array)
         endif

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
 deallocate(array2D)

 end subroutine read_into_fields

!> write variables from ESMF Fields into netcdf restart-like file
!! @param[in] localpet          local pet
!! @param[in] i_dim             longitudinal dimension
!! @param[in] j_dim             latitudinal dimension
!! @param[in] fname_out         file name to write to
!! @param[in] dir_out           directory of file name to write to
!! @param[in] n_vars            number of variables to read in 
!! @param[in] n_tims            number of times to write out
!! @param[in] variable_list     variables to read in
!! @param[in] fields         fields to read variables into
!! @param[in] add_time_dim      specify whether output file has time dimension

 subroutine write_from_fields(localpet, i_dim, j_dim , fname_out, dir_out, &
                                n_vars, n_tims, variable_list, fields, add_time_dim)

 implicit none

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim,  n_vars, n_tims
 character(*), intent(in)        :: fname_out
 character(*), intent(in)        :: dir_out
 character(15), dimension(n_vars), intent(in)     :: variable_list
 type(esmf_field), dimension(n_tims,n_vars), intent(in)  :: fields
 logical,      intent(in)        :: add_time_dim

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

         if (add_time_dim) then ! UFS_UTILS expects input with no time dim
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

             if (add_time_dim) then
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
!! @param[in] res_atm           resolution of grid
!! @param[in] dir_fix           orog fix directory
!! @param[in] localpet          local pet
!! @param[in] npets             total number of pets
!! @param[out] fv3_grid         output ESMF grid 


 subroutine create_grid_fv3(res_atm, dir_fix, npets, localpet, fv3_grid)

! INTENT IN
 integer, intent(in)    :: npets, localpet
 integer, intent(in)    :: res_atm
 character(*), intent(in)       :: dir_fix

 ! INTENT OUT
 type(esmf_grid), intent(out)   :: fv3_grid

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
!! @param[in] grid_setup        data structure with grid details 
!! @param[in] npets             total number of pets
!! @param[in] localpet          local pet
!! @param[out] gauss_grid       output ESMF grid 

 subroutine create_grid_gauss(grid_setup, npets, localpet, gauss_grid)

 ! INTENT IN
 type(grid_setup_type), intent(in) :: grid_setup
 integer, intent(in)   :: npets, localpet

 ! INTENT OUT
 type(esmf_grid)                   :: gauss_grid

 integer :: ierr, fac
 character(len=200)     :: fname

 fname = trim(grid_setup%dir_coord)//trim(grid_setup%fname_coord)

 if (localpet == 0) print*," creating gauss grid for ", trim(fname)

 fac = npets / n_tiles
 gauss_grid = ESMF_GridCreate(filename=trim(fname),  &
              fileFormat=ESMF_FILEFORMAT_SCRIP,  &
              regDecomp=(/n_tiles,fac/), addCornerStagger=.true., rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN Gauss GridCreate", ierr)

 end subroutine create_grid_gauss

end  module grids_IO
