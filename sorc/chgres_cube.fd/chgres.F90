!> @file
!! @brief Initialize an FV3 model run.
!!
!! @author George Gayno NOAA/EMC

!> Initialize an FV3 model run.
!!
!! Initialize an FV3 run using history or restart data from another
!! FV3 run, the spectral GFS, and a few other models. Converts
!! atmospheric, surface and nst data.
!!
!! This file reads a configuration namelist.
!!
!! Link the configuration namelist to ./fort.41. Then run the program
!! with a multiple of six mpi tasks (an ESMF library requirement for
!! fv3 cubed sphere grids).
!!
!! @note For variable names “input” refers to the data input to the
!! program (i.e., GRIB2, NEMSIO, NetCDF). “Target” refers to the
!! target or FV3 model grid.
!!
!! @author George Gayno NOAA/EMC
!! @return 0 for success, error code otherwise.
 program chgres

 use mpi
 use esmf

 use atmosphere, only          : atmosphere_driver

 use program_setup, only       : read_setup_namelist, &
                                 read_varmap, &
                                 convert_atm, &
                                 convert_sfc 

 use model_grid, only          : define_target_grid,  &
                                 define_input_grid, &
                                 cleanup_input_target_grid_data

 use surface, only             : surface_driver

 use utilities, only           : error_handler
 implicit none

 integer                      :: ierr, localpet, npets

 type(esmf_vm)                :: vm

!-------------------------------------------------------------------------
! Initialize mpi and esmf environment.
!-------------------------------------------------------------------------

 call mpi_init(ierr)

 print*,"- INITIALIZE ESMF"
 call ESMF_Initialize(rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", ierr)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", ierr)

 print*,'- NPETS IS  ',npets
 print*,'- LOCAL PET ',localpet

!-------------------------------------------------------------------------
! Read program configuration namelist.
!-------------------------------------------------------------------------

 call read_setup_namelist
 
!-------------------------------------------------------------------------
! Read variable mapping file (used for grib2 input data only).
!-------------------------------------------------------------------------

 call read_varmap

!-------------------------------------------------------------------------
! Create esmf grid objects for input and target grids.
!-------------------------------------------------------------------------

 call define_target_grid(localpet, npets)
 
 call define_input_grid(localpet, npets)

!-------------------------------------------------------------------------
! Convert atmospheric fields
!-------------------------------------------------------------------------

 if (convert_atm) then

   call atmosphere_driver(localpet)

 end if

!-------------------------------------------------------------------------
! Convert surface/nsst fields
!-------------------------------------------------------------------------

 if (convert_sfc) then

   call surface_driver(localpet)

 end if

 call cleanup_input_target_grid_data

 print*,"- CALL ESMF_finalize"
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

 call mpi_finalize(ierr)

 print*,"- DONE."

 end program chgres
