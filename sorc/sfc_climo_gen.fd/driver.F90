!> @file
!! @brief Interpolates static surface data from lat/lon grid,
!! to FV3 model grid.
!! @author George Gayno

!> Reads static surface data on a global lat/lon grid,
!! interpolates the data to the fv3 model grid, and outputs the
!! result in netcdf format.
!!
!! Program execution is controlled by variables defined in the
!! program configuration namelist (see module program_setup for
!! details).
!!
!! Requires the following input files:
!!   1) Model mosaic file (netcdf)
!!   2) Model orography files (netcdf)
!!   3) Model grid files (netcdf)
!!   4) Source data file on global lat/lon grid (netcdf)
!!
!! Outputs surface data on the model tiles in netcdf format.
!!
!! @return 0 for success, error code otherwise.
!! @author George Gayno
 program driver  

 use model_grid
 use source_grid
 use program_setup
 use utils
 use esmf

 implicit none

 integer                      :: rc
 integer                      :: localpet
 integer                      :: npets

 type(esmf_vm)                :: vm
 type(esmf_regridmethod_flag) :: method

!-------------------------------------------------------------------------
! Initialize MPI and ESMF environments.
!-------------------------------------------------------------------------

 call mpi_init(rc)

 print*,"- INITIALIZE ESMF"
 call ESMF_Initialize(rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", rc)

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal.", rc)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet.", rc)

 print*,'- NPETS IS  ',npets
 print*,'- LOCAL PET ',localpet

!-------------------------------------------------------------------------
! Program set up step.
!-------------------------------------------------------------------------

 call read_setup_namelist(localpet)

!-------------------------------------------------------------------------
! Define fv3 model grid.
!-------------------------------------------------------------------------

 call define_model_grid(localpet, npets)

!-------------------------------------------------------------------------
! Interpolate fields.  Vegetation type must be done first as
! it defines which points are permanent land ice.
!-------------------------------------------------------------------------

 call define_source_grid(localpet, npets, input_vegetation_type_file)
 if (fract_vegsoil_type) then
   print*,'- WILL OUTPUT VEGETATION TYPE FRACTION.'
   call interp_frac_cats(localpet, input_vegetation_type_file)
 else
   print*,'- WILL OUTPUT DOMINATE VEGETATION TYPE.'
   method=ESMF_REGRIDMETHOD_NEAREST_STOD
   call interp(localpet, method, input_vegetation_type_file)
 endif
 call source_grid_cleanup

! Snow free albedo

 if (trim(input_snowfree_albedo_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_snowfree_albedo_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   if (trim(snowfree_albedo_method)=="conserve") method=ESMF_REGRIDMETHOD_CONSERVE
   call interp(localpet, method, input_snowfree_albedo_file)
   call source_grid_cleanup
 endif

! Maximum snow albedo

 if (trim(input_maximum_snow_albedo_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_maximum_snow_albedo_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   if (trim(maximum_snow_albedo_method)=="conserve") method=ESMF_REGRIDMETHOD_CONSERVE
   call interp(localpet, method, input_maximum_snow_albedo_file)
   call source_grid_cleanup
 endif

! FACSF - fractional coverage for strong zenith angle
! dependant albedo.

 if (trim(input_facsf_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_facsf_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   call interp(localpet, method, input_facsf_file)
   call source_grid_cleanup
 endif

! Soil substrate temperature

 if (trim(input_substrate_temperature_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_substrate_temperature_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   call interp(localpet, method, input_substrate_temperature_file)
   call source_grid_cleanup
 endif

! Slope type

 if (trim(input_slope_type_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_slope_type_file)
   method=ESMF_REGRIDMETHOD_NEAREST_STOD
   call interp(localpet, method, input_slope_type_file)
   call source_grid_cleanup
 endif

! Soil type

 if (trim(input_soil_type_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_soil_type_file)
   if (fract_vegsoil_type) then
     print*,'- WILL OUTPUT SOIL TYPE FRACTION.'
     call interp_frac_cats(localpet, input_soil_type_file)
   else
     print*,'- WILL OUTPUT DOMINATE SOIL TYPE.'
     method=ESMF_REGRIDMETHOD_NEAREST_STOD
     call interp(localpet, method, input_soil_type_file)
   endif
   call source_grid_cleanup
 endif

! Vegetation greenness

 if (trim(input_vegetation_greenness_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_vegetation_greenness_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   if (trim(vegetation_greenness_method)=="conserve") method=ESMF_REGRIDMETHOD_CONSERVE
   call interp(localpet, method, input_vegetation_greenness_file)
   call source_grid_cleanup
 endif

! Leaf Area Index

 if (trim(input_leaf_area_index_file) /= "NULL") then
   call define_source_grid(localpet, npets, input_leaf_area_index_file)
   method=ESMF_REGRIDMETHOD_BILINEAR
   if (trim(leaf_area_index_method)=="conserve") method=ESMF_REGRIDMETHOD_CONSERVE
   call interp(localpet, method, input_leaf_area_index_file)
   call source_grid_cleanup
 endif

 call model_grid_cleanup

 print*,"- CALL ESMF_finalize"
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=rc)

 call mpi_finalize(rc)

 print*,'- DONE.'
 stop

 end program driver
