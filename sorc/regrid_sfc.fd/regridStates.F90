!> @file
!! @brief Program to re-grid a list of FV3 variables.
!! @author Clara Draper, and George Gayno  Aug, 2024.

!> Program to re-grid a list of FV3 variables.
!! Intended for use in DA applications (regridding of restarts for recentering, regridding increments).
!! @return 0 for success, error code otherwise.

 program regridStates

 use mpi
 use esmf

 use grids_IO, only     : setup_grid, &
                          write_from_fields, &
                          read_into_fields, &
                          n_tiles, &
                          grid_setup_type

 use utilities, only    : error_handler

 implicit none

 integer, parameter             :: max_vars = 10  !< increase if wish to specify more variables

 ! namelist inputs
 character(len=15)              :: variable_list(max_vars)
 integer                        :: n_vars, n_tims, extrap_levs
 integer                        :: time_list(10)               !< increment forecast hours
 logical                        :: add_time_dim                !< specify whether the output increment has time dimension 
 real(esmf_kind_r8)             :: missing_value ! value given to unmapped cells in the output grid
 integer                        :: nmem_ens

 type(grid_setup_type)          :: grid_setup_in, grid_setup_out

 integer                        :: ierr, localpet, npets, localcomm, subpet, imem_ens
 integer                        :: v, t, SRCTERM

 character(100)                 :: fname_time

 type(esmf_vm)                  :: vmlocal
 type(esmf_grid), allocatable   :: grid_in(:)
 type(esmf_grid)                :: grid_out
 type(esmf_field), allocatable  :: fields_in(:,:)
 type(esmf_field), allocatable  :: fields_out(:,:)
 type(esmf_routehandle)         :: regrid_route
 real(esmf_kind_r8), pointer    :: ptr_out(:,:)

 integer :: ut

 real :: t1, t2, t3, t4
 character(len=3)               :: tstr

 ! see README for details of namelist variables.
 namelist /config/ n_vars, variable_list, missing_value, extrap_levs, time_list, add_time_dim, nmem_ens

! INITIALIZE
!-------------------------------------------------------------------------

 call cpu_time(t1)


 ! intialize mpi

 call mpi_init(ierr)
 if (ierr .ne. MPI_SUCCESS) call error_handler("mpi_init", ierr)

 call mpi_comm_rank(MPI_COMM_WORLD, localpet, ierr)
 if (ierr .ne. MPI_SUCCESS) call error_handler("mpi_comm_rank", ierr)

 call mpi_comm_size(MPI_COMM_WORLD, npets, ierr)
 if (ierr .ne. MPI_SUCCESS) call error_handler("mpi_comm_size", ierr)

 if (mod(npets,n_tiles) /= 0) then
   call error_handler("must run with a task count that is a multiple of 6", 1)
 endif

 imem_ens = localpet/n_tiles + 1

 call mpi_comm_split(MPI_COMM_WORLD, imem_ens-1, localpet, localcomm, ierr)
 if (ierr .ne. MPI_SUCCESS) call error_handler("mpi_comm_split", ierr)

 call mpi_comm_rank(localcomm, subpet, ierr)
 if (ierr .ne. MPI_SUCCESS) call error_handler("mpi_comm_rank(localcomm)", ierr)

 ! initialize esmf

 call ESMF_Initialize(rc=ierr, mpiCommunicator=localcomm, logkindflag=ESMF_LOGKIND_MULTI)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

!-------------------------------------------------------------------------
! RUN
!-------------------------------------------------------------------------

 print*,'** pets: local, total: ',localpet, npets

!------------------------
! read in namelist

 ! defaults
 missing_value=-999.
 extrap_levs=2
 time_list=-1

 open(newunit=ut, file='regrid.nml', iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING regrid NAMELIST.", ierr)
 read(ut, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING config NAMELIST.", ierr)
 if(npets/n_tiles /= nmem_ens) then
   call error_handler("number of processor divided by number of tiles must equal number of ensemble members", 1)
 endif
 call readin_setup(ut,"input",nmem_ens,imem_ens,grid_setup_in)
 call readin_setup(ut,"output",nmem_ens,imem_ens,grid_setup_out)
 close (ut)

 n_tims = 0
 do t=1,10
   if (time_list(t) .lt. 0) exit
   n_tims = n_tims + 1
 enddo
 if (n_tims < 1) then
   call error_handler("n_tims < 1. must have at least one valid increment hour in time_list", 1)
 endif

!------------------------
! Create esmf grid objects for input and output grids, and add land masks

! TO DO - can we make the number of tasks more flexible for fv3

 if (subpet==0) print*,'** Setting up grids for ensemble member ', imem_ens
 allocate(grid_in(n_tims))
 do t = 1, n_tims
   if (grid_setup_in%mask_from_input) then
     call setup_grid(subpet, n_tiles, imem_ens, grid_setup_in, grid_in(t), time_list(t) )
   else
     call setup_grid(subpet, n_tiles, imem_ens, grid_setup_in, grid_in(t))
   endif
 enddo
 call setup_grid(subpet, n_tiles, imem_ens, grid_setup_out, grid_out )

!------------------------
! Create input and output fields

 if (subpet==0) print*,'** Creating/Reading fields for ensemble member ', imem_ens

! input
 allocate(fields_in(n_tims,n_vars))

 do t = 1, n_tims
     do v = 1, n_vars

        fields_in(t,v)  = ESMF_FieldCreate(grid_in(t), &
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
        
        write(tstr,"(I3.3)")time_list(t)
        fname_time = trim(grid_setup_in%fname)//tstr//".nc"
        write(6,*) 'reading into ', trim(fname_time)
        call read_into_fields(subpet, grid_setup_in%ires, grid_setup_in%jres, &
                                 trim(fname_time), trim(grid_setup_in%dir), &
                                 grid_setup_in, n_vars, variable_list(1:n_vars), fields_in(t,:))
 enddo

 call cpu_time(t2)
!------------------------
! regrid the input fields to the output grid

 if (subpet==0) print*,'** Performing regridding for ensemble member', imem_ens

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

 if (subpet==0) print*,'** Writing out regridded fields for ensemble member ', imem_ens

 call write_from_fields(subpet, imem_ens, grid_setup_out%ires, grid_setup_out%jres,     &
                          trim(grid_setup_out%fname), trim(grid_setup_out%dir), &
                          n_vars, n_tims, variable_list(1:n_vars), fields_out, add_time_dim)


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
 
     call ESMF_GridDestroy(grid_in(t), rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING GRID", ierr)
 enddo

 call ESMF_GridDestroy(grid_out,rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING GRID", ierr)

!-------------------------------------------------------------------------
! FINALIZE
!-------------------------------------------------------------------------

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)
 call mpi_finalize(ierr)

 call cpu_time(t4)
 if (subpet==0) print*, '** time in tile2tile', t4 - t1, 'for ensemble member ', imem_ens
 if (subpet==0) print*, '** time in RegridStore', t3 - t2, 'for ensemble member ', imem_ens

 print*,"** DONE.", localpet

 end program regridStates
