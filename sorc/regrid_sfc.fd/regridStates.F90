 program regridStates
! Program to re-grid a list of FV3 variables.
! Intended for use in DA applications (regridding of restarts for recentering, regridding increments).
!
! Clara Draper, and George Gayno  Aug, 2024.

! TO DO - add option to pre-compute, and read in regridding route.

 use mpi_f08
 use esmf

 use grids_IO, only     : setup_grid, &
                          write_from_fields, &
                          read_into_fields, &
                          n_tiles, &
                          grid_setup_type

 use utilities, only    : error_handler

 implicit none

 integer, parameter             :: max_vars = 10  ! increase if wish to specify more variables

 ! namelist inputs
 character(len=15)              :: variable_list(max_vars)
 character(len=2)               :: time_list(9)
 integer                        :: n_vars, n_tims, extrap_levs
 real(esmf_kind_r8)             :: missing_value ! value given to unmapped cells in the output grid

 type(grid_setup_type)          :: grid_setup_in, grid_setup_out

 integer                        :: ierr, localpet, npets
 integer                        :: v, t, SRCTERM

 character(100)                 :: fname_time

 type(esmf_vm)                  :: vm
 type(esmf_grid)                :: grid_in, grid_out
 type(esmf_field), allocatable  :: fields_in(:,:)
 type(esmf_field), allocatable  :: fields_out(:,:)
 type(esmf_routehandle)         :: regrid_route
 real(esmf_kind_r8), pointer    :: ptr_out(:,:)

 integer :: ut

 real :: t1, t2, t3, t4

 ! see README for details of namelist variables.
 namelist /config/ n_vars, n_tims, time_list, variable_list, missing_value, extrap_levs

! INITIALIZE
!-------------------------------------------------------------------------

 call cpu_time(t1)

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", ierr)

!-------------------------------------------------------------------------
! RUN
!-------------------------------------------------------------------------

 print*,'** pets: local, total: ',localpet, npets

 ! checks

 if (mod(npets,n_tiles) /= 0) then
   call error_handler("must run with a task count that is a multiple of 6", 1)
 endif

!------------------------
! read in namelist

 ! defaults
 missing_value=-999.
 extrap_levs=2
 n_tims=1

 open(newunit=ut, file='regrid.nml', iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING regrid NAMELIST.", ierr)
 read(ut, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING config NAMELIST.", ierr)
 call readin_setup(ut,"input",grid_setup_in)
 call readin_setup(ut,"output",grid_setup_out)
 close (ut)


!------------------------
! Create esmf grid objects for input and output grids, and add land masks

! TO DO - can we make the number of tasks more flexible for fv3

 if (localpet==0) print*,'** Setting up grids'
 call setup_grid(localpet, npets, grid_setup_in, grid_in )

 call setup_grid(localpet, npets, grid_setup_out, grid_out )

!------------------------
! Create input and output fields

 if (localpet==0) print*,'** Creating/Reading fields'

! input
 allocate(fields_in(n_tims,n_vars))

 do t = 1, n_tims
     do v = 1, n_vars

        fields_in(t,v)  = ESMF_FieldCreate(grid_in, &
                            typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            name="input for regridding", &
                            rc=ierr)

        if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
           call error_handler("in FieldCreate "//trim(variable_list(v)), ierr)
     end do
 end do

! output
 allocate(fields_out(n_tims,n_vars))

 do t = 1, n_tims
     do v = 1, n_vars

         fields_out(t,v)  = ESMF_FieldCreate(grid_out, &
                                           typekind=ESMF_TYPEKIND_R8, &
                                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                                           name="output of regridding", &
                                           rc=ierr)
         if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("in FieldCreate, field_out", ierr)


         ! set the default output value (for non-mapped cells)
         call ESMF_FieldGet(fields_out(t,v), &
                            farrayPtr=ptr_out, &
                            rc=ierr)
         if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", ierr)

         ptr_out=missing_value

     enddo
 enddo

!------------------------
! read data into input fields

 do t = 1, n_tims

        if (n_tims>1) then
                fname_time = trim(grid_setup_in%fname)//"."//time_list(t)
        else
                fname_time = trim(grid_setup_in%fname)
        endif
        write(6,*) 'reading into ', trim(fname_time)
        call read_into_fields(localpet, grid_setup_in%ires, grid_setup_in%jres, &
                                 trim(fname_time), trim(grid_setup_in%dir), &
                                 grid_setup_in, n_vars, variable_list(1:n_vars), fields_in(t,:))
 enddo

 call cpu_time(t2)
!------------------------
! regrid the input fields to the output grid

 if (localpet==0) print*,'** Performing regridding'

 SRCTERM=1
 ! get regriding route for a field (only uses the grid info in the field)
 ! to turn off masking, remove [src/dstMaskVales] argumemnts
 call ESMF_FieldRegridStore(srcField=fields_in(1,1), srcMaskValues=(/0/), &
                            dstField=fields_out(1,1), dstMaskValues=(/0/), &
                            ! allow unmapped grid cells, without returning error
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            ! fill un-mapped grid cells with a neighbour
                            extrapMethod=ESMF_EXTRAPMETHOD_CREEP, &
                            ! number of "levels" of neighbours to search for a value
                            extrapNumLevels=extrap_levs, &
                            ! needed for reproducibility
                            ! (combined with ESMF_TERMORDER_SRCSEQ below)
                            srctermprocessing=SRCTERM, &
                            routehandle=regrid_route, &
                            ! use bilinear interp (slightly better results than PATCH)
                            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", ierr)

! do the re-gridding

 call cpu_time(t3)

 do t=1, n_tims
     do v=1, n_vars
         call ESMF_FieldRegrid(fields_in(t,v), &
                               fields_out(t,v), &
                               routehandle=regrid_route, &
                               zeroregion=ESMF_REGION_SELECT, & ! initialize output with missing_value
                               termorderflag=ESMF_TERMORDER_SRCSEQ, rc=ierr)
         if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegrid", ierr)
     enddo
 enddo

! TO-DO: terrain-correct temperatures (all layers?)

! write out fields on destination grid. All times into same file.

 if (localpet==0) print*,'** Writing out regridded fields'

 call write_from_fields(localpet, grid_setup_out%ires, grid_setup_out%jres,     &
                          trim(grid_setup_out%fname), trim(grid_setup_out%dir), &
                          n_vars, n_tims, variable_list(1:n_vars), fields_out)


! clean up

 call ESMF_FieldRegridRelease(routehandle=regrid_route, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", ierr)

 do t = 1, n_tims
     do v = 1, n_vars
         call ESMF_FieldDestroy(fields_in(t,v),rc=ierr)
         if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                call error_handler("DESTROYING FIELD", ierr)

         call ESMF_FieldDestroy(fields_out(t,v),rc=ierr)
         if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("DESTROYING FIELD", ierr)
     enddo
 enddo

 call ESMF_GridDestroy(grid_in,rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING GRID", ierr)

 call ESMF_GridDestroy(grid_out,rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING GRID", ierr)

!-------------------------------------------------------------------------
! FINALIZE
!-------------------------------------------------------------------------

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)
 call mpi_finalize(ierr)

 call cpu_time(t4)
 if (localpet==0) print*, '** time in tile2tile', t4 - t1
 if (localpet==0) print*, '** time in RegridStore', t3 - t2

 print*,"** DONE.", localpet

 end program regridStates

 subroutine readin_setup(unt,namel,grid_setup)
! subroutine to read in namelists, and convert
! values into setupgrid.
! Also fills in some values, and tests have all
! needed vals, according to the selected grid type.

 use grids_IO, only     : grid_setup_type
 use utilities, only    : error_handler

 implicit none

 ! INPUTS
 integer, intent(in) :: unt
 character(*), intent(in) :: namel
 ! OUTPUTS
 type(grid_setup_type), intent(out) :: grid_setup

 character(len=7)   :: gridtype
 character(len=100) :: fname, fname_mask, fname_coord
 character(len=100) :: dir, dir_mask, dir_coord
 character(len=4)   :: default_str="NULL"
 integer            :: ires, jres
 integer            :: ierr

 namelist /input/  fname, dir, &
                   gridtype, &
                   fname_mask, dir_mask, &
                   fname_coord, dir_coord, &
                   ires, jres

 namelist /output/  fname, dir, &
                    gridtype, &
                    fname_mask, dir_mask, &
                    fname_coord, dir_coord,&
                    ires, jres

 ! set defaults
 fname = default_str
 dir = default_str
 fname_mask = default_str
 dir_mask = default_str
 fname_coord = default_str
 dir_coord = default_str
 ires = 0
 jres = 0

 select case (namel)
 case ("input")
     read(unt, nml=input, iostat=ierr)
     if (ierr /= 0) call error_handler("READING input NAMELIST.", ierr)
 case ("output")
     read(unt, nml=output, iostat=ierr)
     if (ierr /= 0) call error_handler("READING output NAMELIST.", ierr)
 case default
     call error_handler("unknown namel in readin_setup", 1)
 end select

 grid_setup%descriptor = gridtype

 grid_setup%dir = dir
 grid_setup%fname = fname

 ! set-up mask details, based on file type
 select case (gridtype)
 case ("fv3_rst") ! for history file and restarts, use veg type
     grid_setup%dir_mask = dir_mask ! get from a fix file
     grid_setup%fname_mask = fname_mask
     grid_setup%mask_variable(1) =  "vegetation_type" ! if getting from fix file
 case ("gau_inc") ! gsi-output incr files only, use calculated mask
     if (trim(fname_mask) == default_str) then ! if not specified, use input file
         grid_setup%dir_mask = dir
         grid_setup%fname_mask = fname
     else
         grid_setup%dir_mask = dir_mask
         grid_setup%fname_mask = fname_mask
     endif
     grid_setup%mask_variable(1) =  "soilsnow_mask  "
 case default
     call error_handler("unknown gridtype in readin_setup", 1)
 end select

 grid_setup%dir_coord = dir_coord
 grid_setup%fname_coord = fname_coord
 grid_setup%ires = ires
 grid_setup%jres = jres

 end subroutine
!-------------------------------------------------------------------------
