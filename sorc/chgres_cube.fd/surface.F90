!> @file
!! @brief Process land, sea/lake ice, open water/Near Sea Surface
!! Temperature (NSST) fields.
!! @author George Gayno NCEP/EMC

!> Process surface and nst fields. Interpolates fields from the input
!! to target grids. Adjusts soil temperature according to differences
!! in input and target grid terrain. Rescales soil moisture for soil
!! type differences between input and target grid. Computes frozen
!! portion of total soil moisture.
!!
!! Assumes the input land data are Noah LSM-based, and the fv3 run
!! will use the Noah LSM. NSST fields are not available when using
!! GRIB2 input data.
!!
!! Public variables are defined below. "target" indicates field
!! associated with the target grid. "input" indicates field associated
!! with the input grid.
!!
!! @author George Gayno NCEP/EMC
 module surface

 use esmf

 use surface_target_data, only : canopy_mc_target_grid, t2m_target_grid, &
                                 q2m_target_grid, tprcp_target_grid, &
                                 f10m_target_grid, seaice_fract_target_grid, &
                                 ffmm_target_grid, ustar_target_grid, &
                                 srflag_target_grid, soil_temp_target_grid, &
                                 seaice_depth_target_grid, snow_liq_equiv_target_grid, &
                                 seaice_skin_temp_target_grid, skin_temp_target_grid, &
                                 snow_depth_target_grid, z0_target_grid, &
                                 c_d_target_grid, c_0_target_grid, &
                                 d_conv_target_grid, dt_cool_target_grid, &
                                 ifd_target_grid, qrain_target_grid, &
                                 tref_target_grid, w_d_target_grid, &
                                 w_0_target_grid, xs_target_grid, &
                                 xt_target_grid, xu_target_grid, &
                                 xv_target_grid, xz_target_grid, &
                                 xtts_target_grid, xzts_target_grid, &
                                 z_c_target_grid, zm_target_grid, &
                                 soilm_tot_target_grid, lai_target_grid, &
                                 soilm_liq_target_grid, ice_temp_target_grid, &
                                 snow_depth_at_ice_target_grid, alvsf_nl_target_grid, &
                                 alnsf_nl_target_grid, alvwf_nl_target_grid, &
                                 alnwf_nl_target_grid, z0_water_target_grid, &
                                 z0_ice_target_grid, sst_target_grid, &
                                 seaice_substrate_temp_target_grid, &
                                 snow_liq_equiv_at_ice_target_grid

 use write_data, only : write_fv3_sfc_data_netcdf

 use utilities, only  : error_handler

 implicit none

 private

 integer, parameter                 :: veg_type_landice_target = 15
                                       !< Vegetation type category that
                                       !< defines permanent land ice points.
                                       !< The Noah LSM land ice physics
                                       !< are applied at these points.

 type(esmf_field)                   :: soil_type_from_input_grid
                                       !< soil type interpolated from
                                       !< input grid
 type(esmf_field)                   :: terrain_from_input_grid
                                       !< terrain height interpolated
                                       !< from input grid
 type(esmf_field)                   :: terrain_from_input_grid_land
                                       !< terrain height interpolated
                                       !< from input grid at all land points 

 real, parameter, private           :: blim        = 5.5
                                       !< soil 'b' parameter limit
 real, parameter, private           :: frz_h2o     = 273.15
                                       !< melting pt water
 real, parameter, private           :: frz_ice     = 271.21
                                       !< melting pt sea ice
 real, parameter, private           :: grav        = 9.81
                                       !< gravity
 real, parameter, private           :: hlice       = 3.335E5
                                       !< latent heat of fusion
                                       

 type realptr_2d
   real(esmf_kind_r8), pointer :: p(:,:)
                                       !< array of 2d pointers
 end type realptr_2d
                                       !< pointer to hold array of 2d pointers 
  type realptr_3d
   real(esmf_kind_r8), pointer :: p(:,:,:)
                                       !< array of 3d pointers
 end type realptr_3d
                                       !< pointer to hold array of 3d pointers

 public :: surface_driver
 public :: create_nst_esmf_fields
 public :: interp
 public :: create_surface_esmf_fields
 public :: nst_land_fill
 public :: regrid_many
 public :: search_many

 contains

!> Driver routine to process surface/nst data
!!
!! @param[in] localpet  ESMF local persistent execution thread
!!
!! @author George Gayno NCEP/EMC
 subroutine surface_driver(localpet)

 use sfc_input_data, only            : cleanup_input_sfc_data, &
                                       read_input_sfc_data

 use nst_input_data, only            : cleanup_input_nst_data, &
                                       read_input_nst_data

 use program_setup, only             : calc_soil_params_driver, &
                                       convert_nst
                                  
 use static_data, only               :  get_static_fields, &
                                       cleanup_static_fields

 use surface_target_data, only       : cleanup_target_nst_data

 use utilities, only                 : error_handler

 implicit none

 integer, intent(in)                :: localpet

!-----------------------------------------------------------------------
! Compute soil-based parameters.
!-----------------------------------------------------------------------

 call calc_soil_params_driver(localpet)

!-----------------------------------------------------------------------
! Get static data (like vegetation type) on the target grid.
!-----------------------------------------------------------------------

 call get_static_fields(localpet)

!-----------------------------------------------------------------------
! Read surface data on input grid.
!-----------------------------------------------------------------------

 call read_input_sfc_data(localpet)

!-----------------------------------------------------------------------
! Read nst data on input grid.
!-----------------------------------------------------------------------

 if (convert_nst) call read_input_nst_data(localpet)

!-----------------------------------------------------------------------
! Create surface field objects for target grid.
!-----------------------------------------------------------------------

 call create_surface_esmf_fields

!-----------------------------------------------------------------------
! Create nst field objects for target grid.
!-----------------------------------------------------------------------

 if (convert_nst) call create_nst_esmf_fields
 
!-----------------------------------------------------------------------
! Adjust soil levels of input grid !! not implemented yet
!-----------------------------------------------------------------------

 call adjust_soil_levels(localpet)

!-----------------------------------------------------------------------
! Horizontally interpolate fields.
!-----------------------------------------------------------------------

 call interp(localpet)
 
!---------------------------------------------------------------------------------------------
! Adjust soil/landice column temperatures for any change in elevation between  the
! input and target grids.
!---------------------------------------------------------------------------------------------

 call adjust_soilt_for_terrain
 
!---------------------------------------------------------------------------------------------
! Rescale soil moisture for changes in soil type between the input and target grids.
!---------------------------------------------------------------------------------------------

 call rescale_soil_moisture
 
!---------------------------------------------------------------------------------------------
! Compute liquid portion of total soil moisture.
!---------------------------------------------------------------------------------------------

 call calc_liq_soil_moisture

!---------------------------------------------------------------------------------------------
! Set z0 at land and sea ice.
!---------------------------------------------------------------------------------------------

 call roughness

!---------------------------------------------------------------------------------------------
! Perform some final qc checks.
!---------------------------------------------------------------------------------------------

 call qc_check

!---------------------------------------------------------------------------------------------
! Set flag values at land for nst fields.
!---------------------------------------------------------------------------------------------

 if (convert_nst) call nst_land_fill

!---------------------------------------------------------------------------------------------
! Free up memory.
!---------------------------------------------------------------------------------------------

 call cleanup_input_sfc_data

 if (convert_nst) call cleanup_input_nst_data

 
 call update_landmask

!---------------------------------------------------------------------------------------------
! Write data to file.
!---------------------------------------------------------------------------------------------

 call write_fv3_sfc_data_netcdf(localpet)

!---------------------------------------------------------------------------------------------
! Free up memory.
!---------------------------------------------------------------------------------------------

 if (convert_nst) call cleanup_target_nst_data

 call cleanup_all_target_sfc_data

 call cleanup_static_fields

 return

 end subroutine surface_driver

!> Horizontally interpolate surface fields from input to target FV3
!> grid using esmf routines.
!!
!! @param[in] localpet  ESMF local persistent execution thread
!!
!! @author George Gayno NOAA/EMC
 subroutine interp(localpet)

 use mpi
 use esmf

 use sfc_input_data, only            : canopy_mc_input_grid,  &
                                       f10m_input_grid,  &
                                       ffmm_input_grid,  &
                                       landsea_mask_input_grid, &
                                       q2m_input_grid,  &
                                       seaice_depth_input_grid, &
                                       seaice_fract_input_grid, &
                                       seaice_skin_temp_input_grid, &
                                       skin_temp_input_grid, &
                                       snow_depth_input_grid, &
                                       snow_liq_equiv_input_grid, &
                                       soil_temp_input_grid, &
                                       soil_type_input_grid, &
                                       soilm_tot_input_grid, &
                                       srflag_input_grid, &
                                       t2m_input_grid,  &
                                       tprcp_input_grid,  &
                                       ustar_input_grid,  &
                                       veg_type_input_grid, &
                                       z0_input_grid, &
                                       veg_type_landice_input, &
                                       veg_greenness_input_grid, &
                                       max_veg_greenness_input_grid, &
                                       min_veg_greenness_input_grid, &
                                       lai_input_grid

 use nst_input_data, only           :  c_d_input_grid, &
                                       c_0_input_grid, &
                                       d_conv_input_grid, &
                                       dt_cool_input_grid, &
                                       ifd_input_grid, &
                                       qrain_input_grid, &
                                       tref_input_grid, &
                                       w_d_input_grid, &
                                       w_0_input_grid, &
                                       xs_input_grid, &
                                       xt_input_grid, &
                                       xu_input_grid, &
                                       xv_input_grid, &
                                       xz_input_grid, &
                                       xtts_input_grid, &
                                       xzts_input_grid, &
                                       z_c_input_grid, &
                                       zm_input_grid

 use atm_input_data, only            : terrain_input_grid

 use model_grid, only                : input_grid, target_grid, &
                                       i_target, j_target, &
                                       lsoil_target, &
                                       num_tiles_target_grid, &
                                       landmask_target_grid, &
                                       seamask_target_grid,  &
                                       latitude_target_grid

 use program_setup, only             : convert_nst, &
                                       vgtyp_from_climo, & 
                                       sotyp_from_climo, &
                                       vgfrc_from_climo, &
                                       minmax_vgfrc_from_climo, &
                                       lai_from_climo, &
                                       tg3_from_soil, &
                                       fract_grid
                                       
 use static_data, only               : veg_type_target_grid, &
                                       soil_type_target_grid, &
                                       veg_greenness_target_grid, &
                                       substrate_temp_target_grid,&
                                       min_veg_greenness_target_grid,&
                                       max_veg_greenness_target_grid

 use search_util

 implicit none

 integer, intent(in)                :: localpet

 integer                            :: l(1), u(1)
 integer                            :: i, j, ij, rc, tile
 integer                            :: clb_target(2), cub_target(2)
 integer                            :: isrctermprocessing
 integer                            :: num_fields
 integer                            :: vgfrc_ind, mmvg_ind, lai_ind
 integer, allocatable               :: search_nums(:)
 integer(esmf_kind_i4), pointer     :: unmapped_ptr(:)
 integer(esmf_kind_i4), pointer     :: mask_input_ptr(:,:)
 integer(esmf_kind_i4), pointer     :: mask_target_ptr(:,:)
 integer(esmf_kind_i8), pointer     :: landmask_target_ptr(:,:)
 integer(esmf_kind_i8), allocatable :: mask_target_one_tile(:,:)
 integer(esmf_kind_i8), allocatable :: water_target_one_tile(:,:)
 integer(esmf_kind_i8), allocatable :: land_target_one_tile(:,:)
 integer(esmf_kind_i8), pointer     :: seamask_target_ptr(:,:)

 real(esmf_kind_r8), allocatable    :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable    :: data_one_tile2(:,:)
 real(esmf_kind_r8), allocatable    :: data_one_tile_3d(:,:,:)
 real(esmf_kind_r8), allocatable    :: latitude_one_tile(:,:)
 real(esmf_kind_r8), allocatable    :: fice_target_one_tile(:,:)
 real(esmf_kind_r8), pointer        :: seaice_fract_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: srflag_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: terrain_from_input_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: soil_type_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: landmask_input_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_input_ptr(:,:)
 real(esmf_kind_r8), allocatable    :: veg_type_target_one_tile(:,:)
 real(esmf_kind_r8)                 :: ice_lim

 type(esmf_regridmethod_flag)       :: method
 type(esmf_routehandle)             :: regrid_bl_no_mask
 type(esmf_routehandle)             :: regrid_all_land
 type(esmf_routehandle)             :: regrid_land
 type(esmf_routehandle)             :: regrid_landice
 type(esmf_routehandle)             :: regrid_nonland
 type(esmf_routehandle)             :: regrid_seaice
 type(esmf_routehandle)             :: regrid_water
 
 type(esmf_fieldbundle)             :: bundle_all_target, bundle_all_input
 type(esmf_fieldbundle)             :: bundle_seaice_target, bundle_seaice_input
 type(esmf_fieldbundle)             :: bundle_water_target, bundle_water_input
 type(esmf_fieldbundle)             :: bundle_allland_target, bundle_allland_input
 type(esmf_fieldbundle)             :: bundle_landice_target, bundle_landice_input
 type(esmf_fieldbundle)             :: bundle_nolandice_target, bundle_nolandice_input
 
 logical, allocatable               :: dozero(:)

!-----------------------------------------------------------------------
! Interpolate fieids that do not require 'masked' interpolation.
!-----------------------------------------------------------------------

 method=ESMF_REGRIDMETHOD_BILINEAR

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore FOR NON-MASKED BILINEAR INTERPOLATION."
 call ESMF_FieldRegridStore(t2m_input_grid, &
                            t2m_target_grid, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl_no_mask, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 bundle_all_target = ESMF_FieldBundleCreate(name="all points target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_all_input = ESMF_FieldBundleCreate(name="all points input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
      
 call ESMF_FieldBundleAdd(bundle_all_target, (/t2m_target_grid,q2m_target_grid,tprcp_target_grid, &
                         f10m_target_grid,ffmm_target_grid,ustar_target_grid,srflag_target_grid/), &
                         rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)                        
 call ESMF_FieldBundleAdd(bundle_all_input, (/t2m_input_grid,q2m_input_grid,tprcp_input_grid, &
                         f10m_input_grid,ffmm_input_grid,ustar_input_grid,srflag_input_grid/), &
                         rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)  
                            
 call ESMF_FieldBundleGet(bundle_all_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)  
 
 allocate(dozero(num_fields))
 dozero(:) = .True.  

 call regrid_many(bundle_all_input,bundle_all_target,num_fields,regrid_bl_no_mask,dozero)
 deallocate(dozero) 

 call ESMF_FieldBundleDestroy(bundle_all_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)  
 call ESMF_FieldBundleDestroy(bundle_all_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)  
 
 print*,"- CALL FieldGet FOR SRFLAG."
 call ESMF_FieldGet(srflag_target_grid, &
                    farrayPtr=srflag_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!-----------------------------------------------------------------------
! This is a flag field.  Using neighbor was expensive.  So use
! bilinear and 'nint'.
!-----------------------------------------------------------------------

 srflag_target_ptr = nint(srflag_target_ptr)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl_no_mask, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!-----------------------------------------------------------------------
! First, set the mask on the target and input grids.
!-----------------------------------------------------------------------

 print*,"- CALL GridAddItem FOR TARGET GRID."
 call ESMF_GridAddItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddItem", rc)

 print*,"- CALL GridGetItem FOR TARGET GRID."
 call ESMF_GridGetItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetItem", rc)

 print*,"- CALL FieldGet FOR TARGET GRID SEAMASK."
 call ESMF_FieldGet(seamask_target_grid, &
                    computationalLBound=clb_target, &
                    computationalUBound=cub_target, &
                    farrayPtr=seamask_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 print*,"- CALL FieldGet FOR TARGET GRID LANDMASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    farrayPtr=landmask_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)   

 print*,"- CALL GridAddItem FOR INPUT GRID SEAMASK."
 call ESMF_GridAddItem(input_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddItem", rc)

 print*,"- CALL FieldGet FOR INPUT GRID LANDMASK."
 call ESMF_FieldGet(landsea_mask_input_grid, &
                    farrayPtr=landmask_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL GridGetItem FOR INPUT GRID LANDMASK."
 call ESMF_GridGetItem(input_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetItem", rc)

  if (localpet == 0) then
   allocate(data_one_tile(i_target,j_target))
   allocate(data_one_tile_3d(i_target,j_target,lsoil_target))
   allocate(mask_target_one_tile(i_target,j_target))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_3d(0,0,0))
   allocate(mask_target_one_tile(0,0))
 endif
    
 !-----------------------------------------------------------------------
 ! Interpolate vegetation type to target grid if chosen in namelist and terrain
 ! for use in replacing isolated bad terrain values
 !-----------------------------------------------------------------------
 
 method=ESMF_REGRIDMETHOD_NEAREST_STOD

 isrctermprocessing = 1
 
 mask_input_ptr = 0
 where (nint(landmask_input_ptr) == 1) mask_input_ptr = 1

 mask_target_ptr = 0
 where (landmask_target_ptr == 1) mask_target_ptr = 1
 
 print*,"- CALL FieldCreate FOR TERRAIN FROM INPUT GRID LAND."
 terrain_from_input_grid_land = ESMF_FieldCreate(target_grid, &
                                           typekind=ESMF_TYPEKIND_R8, &
                                           staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)
 
 print*,"- CALL FieldRegridStore for land fields."
 call ESMF_FieldRegridStore(terrain_input_grid, &
                            terrain_from_input_grid_land, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_all_land, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 print*,"- CALL Field_Regrid TERRAIN."
 call ESMF_FieldRegrid(terrain_input_grid, &
                       terrain_from_input_grid_land, &
                       routehandle=regrid_all_land, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ,  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)
    
 print*,"- CALL FieldGet FOR terrain from input grid at land."
 call ESMF_FieldGet(terrain_from_input_grid_land, &
                    farrayPtr=terrain_from_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 l = lbound(unmapped_ptr)
 u = ubound(unmapped_ptr)

 do ij = l(1), u(1)
   call ij_to_i_j(unmapped_ptr(ij), i_target, j_target, i, j)
   terrain_from_input_ptr(i,j) = -9999.9 
 enddo
  nullify(terrain_from_input_ptr) 
  
 do tile = 1, num_tiles_target_grid
 
   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(landmask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
      
   print*,"- CALL FieldGather FOR TERRAIN FROM INPUT GRID: ", tile
   call ESMF_FieldGather(terrain_from_input_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   if (localpet == 0) then
     allocate(land_target_one_tile(i_target,j_target))
     land_target_one_tile = 0
     where(mask_target_one_tile == 1) land_target_one_tile = 1
     call search(data_one_tile, land_target_one_tile, i_target, j_target, tile, 7)
     deallocate(land_target_one_tile)
   endif

   print*,"- CALL FieldScatter FOR TERRAIN FROM INPUT GRID: ", tile
   call ESMF_FieldScatter(terrain_from_input_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo
 
 if(.not. vgtyp_from_climo) then
  
   print*,"- CALL FieldRegrid VEG TYPE."
   call ESMF_FieldRegrid(veg_type_input_grid, &
                         veg_type_target_grid, &
                         routehandle=regrid_all_land, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)
   
   print*,"- CALL FieldGet FOR TARGET grid veg type."
   call ESMF_FieldGet(veg_type_target_grid, &
                      farrayPtr=veg_type_target_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)
        
   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)

   do ij = l(1), u(1)
     call ij_to_i_j(unmapped_ptr(ij), i_target, j_target, i, j)
     veg_type_target_ptr(i,j) = -9999.9 
   enddo

   do tile = 1, num_tiles_target_grid
     print*,"- CALL FieldGather FOR TARGET GRID VEG TYPE TILE: ", tile
     call ESMF_FieldGather(veg_type_target_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", rc)

     print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
     call ESMF_FieldGather(landmask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", rc)

     if (localpet == 0) then
       allocate(land_target_one_tile(i_target,j_target))
       land_target_one_tile = 0
       where(mask_target_one_tile == 1) land_target_one_tile = 1
       call search(data_one_tile, land_target_one_tile, i_target, j_target, tile, 225)
       deallocate(land_target_one_tile)
     endif

     print*,"- CALL FieldScatter FOR TARGET GRID VEG TYPE: ", tile
     call ESMF_FieldScatter(veg_type_target_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)
   enddo
   nullify(veg_type_target_ptr) 
 endif

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_all_land, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)
    
!-----------------------------------------------------------------------
! Next, determine the sea ice fraction on target grid.  
! For fractional grids, the ice fraction is not scaled for the
! fraction of non-land. So if a point is 50% land and non-land,
! an ice frac of 100% means the entire non-land portion is ice covered.
!-----------------------------------------------------------------------

 mask_input_ptr = 1
 where (nint(landmask_input_ptr) == 1) mask_input_ptr = 0
 
! For non-fractional grids, 'seamask_target_ptr' is '0' for land points
! and '1' for non-land points. For fractional grids, 'seamask_target_ptr'
! will be '0' if all land and '1' is at least some non-land.

 mask_target_ptr = int(seamask_target_ptr,kind=esmf_kind_i4)

 method=ESMF_REGRIDMETHOD_CONSERVE

 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for sea ice fraction."
 call ESMF_FieldRegridStore(seaice_fract_input_grid, &
                            seaice_fract_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_nonland, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 print*,"- CALL Field_Regrid for sea ice fraction."
 call ESMF_FieldRegrid(seaice_fract_input_grid, &
                       seaice_fract_target_grid, &
                       routehandle=regrid_nonland, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL FieldGet FOR TARGET grid sea ice fraction."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=seaice_fract_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 l = lbound(unmapped_ptr)
 u = ubound(unmapped_ptr)

 do ij = l(1), u(1)
   call ij_to_i_j(unmapped_ptr(ij), i_target, j_target, i, j)
   seaice_fract_target_ptr(i,j) = -9999.9 ! flag value for missing point
                               ! which will be replaced in routine
                               ! "search".
 enddo

 if (localpet == 0) then
   allocate(latitude_one_tile(i_target,j_target))
 else
   allocate(latitude_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid

   print*,"- CALL FieldGather FOR TARGET GRID SEAICE FRACTION TILE: ", tile
   call ESMF_FieldGather(seaice_fract_target_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
   
   print*,"- CALL FieldGather FOR TARGET GRID MASK TILE: ", tile
   call ESMF_FieldGather(seamask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   print*,"- CALL FieldGather FOR TARGET LATITUDE TILE: ", tile
   call ESMF_FieldGather(latitude_target_grid, latitude_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   if (localpet == 0) then
     call search(data_one_tile, mask_target_one_tile, i_target, j_target, tile, 91, &
                 latitude=latitude_one_tile)
   endif
   
! When running with fractional grids, to reduce the number of points with small amounts of open water, 
! set any point with ice between 95-100% to 100%.

   if (fract_grid) then
     ice_lim = 0.95_esmf_kind_r8
   else
     ice_lim = 1.0_esmf_kind_r8
   endif

   if (localpet == 0) then
     do j = 1, j_target
     do i = 1, i_target
       if (data_one_tile(i,j) > ice_lim) then
         data_one_tile(i,j) = 1.0_esmf_kind_r8
       endif
       if (data_one_tile(i,j) < 0.15_esmf_kind_r8) data_one_tile(i,j) = 0.0_esmf_kind_r8
     enddo
     enddo
   endif

   print*,"- CALL FieldScatter FOR TARGET GRID SEAICE FRACTION TILE: ", tile
   call ESMF_FieldScatter(seaice_fract_target_grid, data_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo

 deallocate(latitude_one_tile)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_nonland, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!---------------------------------------------------------------------------------------------
! Now interpolate other sea ice related fields.  Since we know what points are ice on
! the target grid, reset the target grid mask.
!---------------------------------------------------------------------------------------------

 mask_input_ptr = 0
 where (nint(landmask_input_ptr) == 2) mask_input_ptr = 1

 mask_target_ptr = 0 
 do j = clb_target(2), cub_target(2)
 do i = clb_target(1), cub_target(1)
   if (seaice_fract_target_ptr(i,j) > 0.0) mask_target_ptr(i,j) = 1
 enddo
 enddo

 method=ESMF_REGRIDMETHOD_NEAREST_STOD
 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for 3d seaice fields."
 if(fract_grid)then
 call ESMF_FieldRegridStore(soil_temp_input_grid, &
                            ice_temp_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_seaice, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 else
 call ESMF_FieldRegridStore(soil_temp_input_grid, &
                            soil_temp_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_seaice, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 endif
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 bundle_seaice_target = ESMF_FieldBundleCreate(name="sea ice target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_seaice_input = ESMF_FieldBundleCreate(name="sea ice input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)

 if (fract_grid) then
   call ESMF_FieldBundleAdd(bundle_seaice_target, (/seaice_depth_target_grid, snow_depth_at_ice_target_grid, &
                          snow_liq_equiv_at_ice_target_grid, seaice_skin_temp_target_grid, &
                          ice_temp_target_grid/), rc=rc)
 else
   call ESMF_FieldBundleAdd(bundle_seaice_target, (/seaice_depth_target_grid, snow_depth_target_grid, &
                          snow_liq_equiv_target_grid, seaice_skin_temp_target_grid, &
                          soil_temp_target_grid/), rc=rc)
 endif

  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)

 call ESMF_FieldBundleAdd(bundle_seaice_input, (/seaice_depth_input_grid, snow_depth_input_grid, &
                          snow_liq_equiv_input_grid, seaice_skin_temp_input_grid, &
                          soil_temp_input_grid/), rc=rc)                          
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 call ESMF_FieldBundleGet(bundle_seaice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)
 

 allocate(search_nums(num_fields))
 allocate(dozero(num_fields))

 search_nums = (/92,66,65,21,21/)
 dozero(:) = .True.
 
 l = lbound(unmapped_ptr)
 u = ubound(unmapped_ptr)
 
 call regrid_many(bundle_seaice_input,bundle_seaice_target,num_fields,regrid_seaice,dozero, &
                  unmapped_ptr=unmapped_ptr )
 deallocate(dozero)                 
 call ESMF_FieldBundleDestroy(bundle_seaice_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 if (localpet == 0) then
   allocate(fice_target_one_tile(i_target,j_target))
 else
   allocate(fice_target_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid

   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(landmask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(seaice_fract_target_grid, fice_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

!cfract using ice flag of '2' here. cant do that.
   if (localpet == 0) then   
     where(mask_target_one_tile == 1) mask_target_one_tile = 0
     where(fice_target_one_tile > 0.0) mask_target_one_tile = 1
     call search_many(num_fields,bundle_seaice_target,tile, search_nums,localpet, &
                    mask=mask_target_one_tile)
   else
     call search_many(num_fields,bundle_seaice_target, tile,search_nums,localpet)
   endif

 enddo

 !deallocate(search_nums, fice_target_one_tile)
 deallocate(search_nums)

 call ESMF_FieldBundleDestroy(bundle_seaice_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleDestroy", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_seaice, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!---------------------------------------------------------------------------------------------
! Now interpolate open water fields. 
!---------------------------------------------------------------------------------------------

 mask_input_ptr = 0
 where (nint(landmask_input_ptr) == 0) mask_input_ptr = 1

!cfract include points that are fractional ice.
!cfract use seamask_target, which is 1 if there is some water.
!cfract then remove points where fice is > 0.
!cfract We want points with at least some water, but no ice.

 mask_target_ptr = 0
 where (seamask_target_ptr == 1) mask_target_ptr = 1
 if (fract_grid) then
   where (seaice_fract_target_ptr  == 1.0_esmf_kind_r8) mask_target_ptr = 0
 else
   where (seaice_fract_target_ptr  > 0.0) mask_target_ptr = 0
 endif

 method=ESMF_REGRIDMETHOD_CONSERVE
 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for water fields."
 if(fract_grid)then
 call ESMF_FieldRegridStore(skin_temp_input_grid, &
                            sst_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_water, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 else
 call ESMF_FieldRegridStore(skin_temp_input_grid, &
                            skin_temp_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_water, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 endif
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 bundle_water_target = ESMF_FieldBundleCreate(name="water target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_water_input = ESMF_FieldBundleCreate(name="water input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)

 if(fract_grid)then
 call ESMF_FieldBundleAdd(bundle_water_target, (/sst_target_grid, z0_water_target_grid/), rc=rc)
 else
 call ESMF_FieldBundleAdd(bundle_water_target, (/skin_temp_target_grid, z0_target_grid/), rc=rc)
 endif

  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)

 call ESMF_FieldBundleAdd(bundle_water_input, (/skin_temp_input_grid, z0_input_grid/), rc=rc)  
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)

 if (convert_nst) then

   call ESMF_FieldBundleAdd(bundle_water_target, (/c_d_target_grid,c_0_target_grid,d_conv_target_grid, &
                            dt_cool_target_grid,ifd_target_grid,qrain_target_grid,tref_target_grid, &
                            w_d_target_grid,w_0_target_grid,xs_target_grid,xt_target_grid,xu_target_grid, &
                            xv_target_grid,xz_target_grid,xtts_target_grid,xzts_target_grid, &
                            z_c_target_grid,zm_target_grid/), rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleAdd", rc)
        
   call ESMF_FieldBundleAdd(bundle_water_input, (/c_d_input_grid,c_0_input_grid,d_conv_input_grid, &
                            dt_cool_input_grid,ifd_input_grid,qrain_input_grid,tref_input_grid, &
                            w_d_input_grid,w_0_input_grid,xs_input_grid,xt_input_grid,xu_input_grid, &
                            xv_input_grid,xz_input_grid,xtts_input_grid,xzts_input_grid, &
                            z_c_input_grid,zm_input_grid/), rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleGet(bundle_water_target,fieldCount=num_fields,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleGet", rc)
        
   allocate(search_nums(num_fields))
   allocate(dozero(num_fields))
     
   search_nums(:)=(/11,83,0,0,0,0,1,0,11,0,0,0,0,0,0,30,0,0,0,0/)
   dozero(:) = .True.
             
 else
   call ESMF_FieldBundleGet(bundle_water_target,fieldCount=num_fields,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleGet", rc)
        
   allocate(search_nums(num_fields))
   allocate(dozero(num_fields))
   search_nums(:)=(/11,83/)
   dozero(:) = .True.
 endif

 call regrid_many(bundle_water_input,bundle_water_target,num_fields,regrid_water,dozero, &
                  unmapped_ptr=unmapped_ptr, resetifd=.True.)
 deallocate(dozero)                 
 call ESMF_FieldBundleDestroy(bundle_water_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)


 if (localpet == 0) then
   allocate(latitude_one_tile(i_target,j_target))
 else
   allocate(latitude_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid

   print*,"- CALL FieldGather FOR TARGET SEAMASK TILE: ", tile
   call ESMF_FieldGather(seamask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(seaice_fract_target_grid, fice_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   print*,"- CALL FieldGather FOR TARGET LATITUDE TILE: ", tile
   call ESMF_FieldGather(latitude_target_grid, latitude_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

!cfract - here mask must be points with some water, but no ice.
   if (localpet == 0) then
     allocate(water_target_one_tile(i_target,j_target))
     water_target_one_tile = 0
     where(mask_target_one_tile == 1) water_target_one_tile = 1
     if(fract_grid) then
       where(fice_target_one_tile == 1.0_esmf_kind_r8) water_target_one_tile = 0
     else
       where(fice_target_one_tile > 0.0) water_target_one_tile = 0
     endif
     call search_many(num_fields,bundle_water_target, tile,search_nums,localpet, &
                latitude=latitude_one_tile,mask=water_target_one_tile)
   else
     call search_many(num_fields,bundle_water_target, tile,search_nums,localpet)
   endif

   if (localpet == 0) deallocate(water_target_one_tile)

 enddo

 deallocate(latitude_one_tile,search_nums,fice_target_one_tile)
 
 call ESMF_FieldBundleDestroy(bundle_water_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_water, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!---------------------------------------------------------------------------------------------
! Now interpolate "all land" to "all land".  Here, "all land" means landice and non-land ice.
!---------------------------------------------------------------------------------------------

 mask_input_ptr = 0
 where (nint(landmask_input_ptr) == 1) mask_input_ptr = 1

!cfract this logic should work for fractional grids.
 mask_target_ptr = 0
 where (landmask_target_ptr == 1) mask_target_ptr = 1

 method=ESMF_REGRIDMETHOD_CONSERVE
 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for land fields."
 call ESMF_FieldRegridStore(snow_depth_input_grid, &
                            snow_depth_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_all_land, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 bundle_allland_target = ESMF_FieldBundleCreate(name="all land target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_allland_input = ESMF_FieldBundleCreate(name="all land input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 call ESMF_FieldBundleAdd(bundle_allland_target, (/canopy_mc_target_grid, snow_depth_target_grid, &
                          snow_liq_equiv_target_grid/), rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 call ESMF_FieldBundleAdd(bundle_allland_input, (/canopy_mc_input_grid, snow_depth_input_grid, &
                          snow_liq_equiv_input_grid/), rc=rc)                          
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 call ESMF_FieldBundleGet(bundle_allland_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)
 
 allocate(search_nums(num_fields))
 allocate(dozero(num_fields))

 search_nums = (/223,66,65/)
 dozero=(/.True.,.False.,.False./)
 
 call regrid_many(bundle_allland_input,bundle_allland_target,num_fields,regrid_all_land,dozero, &
                  unmapped_ptr=unmapped_ptr)
 deallocate(dozero)
 call ESMF_FieldBundleDestroy(bundle_allland_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)
  

 do tile = 1, num_tiles_target_grid

   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(landmask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

!cfract 0 - all water; 1 - some land.
   if (localpet == 0) then
     allocate(land_target_one_tile(i_target,j_target))
     land_target_one_tile = 0
     where(mask_target_one_tile == 1) land_target_one_tile = 1
   
     call search_many(num_fields,bundle_allland_target, &
                    tile,search_nums,localpet, mask=land_target_one_tile)
   else
     call search_many(num_fields,bundle_allland_target, tile,search_nums,localpet)
   endif

   if (localpet == 0) deallocate(land_target_one_tile)   
 enddo

 deallocate(search_nums) 
 call ESMF_FieldBundleDestroy(bundle_allland_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)
      
 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_all_land, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!---------------------------------------------------------------------------------------------
! Now interpolate landice points to landice points.
!---------------------------------------------------------------------------------------------

 print*,"- CALL FieldGet FOR INPUT GRID VEG TYPE."
 call ESMF_FieldGet(veg_type_input_grid, &
                    farrayPtr=veg_type_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,'land ice check ',veg_type_landice_input

 mask_input_ptr = 0
 where (nint(veg_type_input_ptr) == veg_type_landice_input) mask_input_ptr = 1

 print*,"- CALL FieldGet FOR TARGET GRID VEG TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract veg_type_target should contain a valid value at points with
!cfract any land. So this logic should work with fractional grids.
 mask_target_ptr = 0
 where (nint(veg_type_target_ptr) == veg_type_landice_target) mask_target_ptr = 1

 method=ESMF_REGRIDMETHOD_NEAREST_STOD
 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for landice fields."
 call ESMF_FieldRegridStore(soil_temp_input_grid, &
                            soil_temp_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_landice, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

  bundle_landice_target = ESMF_FieldBundleCreate(name="landice target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 bundle_landice_input = ESMF_FieldBundleCreate(name="landice input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
 call ESMF_FieldBundleAdd(bundle_landice_target, (/skin_temp_target_grid, terrain_from_input_grid,& 
                          soil_temp_target_grid/), rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 call ESMF_FieldBundleAdd(bundle_landice_input, (/skin_temp_input_grid, terrain_input_grid,&
                          soil_temp_input_grid/), rc=rc)                       
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)
 
 if (.not. sotyp_from_climo) then
   call ESMF_FieldBundleAdd(bundle_landice_input, (/soil_type_input_grid/),rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleAdd(bundle_landice_target,(/soil_type_target_grid/),rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc)
 endif  

 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 call ESMF_FieldBundleGet(bundle_landice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)

 allocate(search_nums(num_fields))
 allocate(dozero(num_fields))
 
 if (sotyp_from_climo) then
   search_nums = (/21,7,21/)
   dozero(:)=.False.
 else
   search_nums = (/21,7,21,231/)
   dozero(:)=(/.False.,.False.,.False.,.True./)
 endif

 call regrid_many(bundle_landice_input,bundle_landice_target,num_fields,regrid_landice,dozero, &
                  unmapped_ptr=unmapped_ptr )  
 deallocate(dozero)                                
 call ESMF_FieldBundleDestroy(bundle_landice_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 if (localpet == 0) then
   allocate (veg_type_target_one_tile(i_target,j_target))
   allocate (land_target_one_tile(i_target,j_target))
   allocate (data_one_tile2(i_target,j_target))
 else
   allocate (veg_type_target_one_tile(0,0))
   allocate (land_target_one_tile(0,0))
   allocate (data_one_tile2(0,0))
 endif

 do tile = 1, num_tiles_target_grid
   print*,"- CALL FieldGather FOR TARGET VEG TYPE TILE: ", tile
   call ESMF_FieldGather(veg_type_target_grid, veg_type_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

!cfract setting land_target_one_tile. this should work for fractional grids.

   if (localpet == 0) then
     land_target_one_tile = 0
     where(nint(veg_type_target_one_tile) == veg_type_landice_target) land_target_one_tile = 1
   endif
   
   print*,"- CALL FieldGather FOR TERRAIN FROM INPUT GRID LAND, TILE: ", tile
    call ESMF_FieldGather(terrain_from_input_grid_land, data_one_tile2, rootPet=0, tile=tile, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
       call error_handler("IN FieldGather", rc)

   if (localpet==0) then
      call search_many(num_fields,bundle_landice_target,tile,search_nums,localpet,&
                terrain_land=data_one_tile2,mask=land_target_one_tile)
   else
      call search_many(num_fields,bundle_landice_target,tile,search_nums,localpet)
   endif
 enddo

 deallocate (veg_type_target_one_tile)
 deallocate (land_target_one_tile)
 deallocate(search_nums)
 
 call ESMF_FieldBundleDestroy(bundle_landice_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_landice, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!---------------------------------------------------------------------------------------------
! Now interpolate land (not including landice pts) to land (not including landice).
!---------------------------------------------------------------------------------------------

 mask_input_ptr = 0
 where (nint(landmask_input_ptr) == 1) mask_input_ptr = 1
 where (nint(veg_type_input_ptr) == veg_type_landice_input) mask_input_ptr = 0

!cfract This should work for fractional grid.
 mask_target_ptr = 0
 where (landmask_target_ptr == 1) mask_target_ptr = 1
 where (nint(veg_type_target_ptr) == veg_type_landice_target) mask_target_ptr = 0

 method=ESMF_REGRIDMETHOD_NEAREST_STOD
 isrctermprocessing = 1

 print*,"- CALL FieldRegridStore for 3d land (but no land ice) fields."
 call ESMF_FieldRegridStore(soilm_tot_input_grid, &
                            soilm_tot_target_grid, &
                            srcmaskvalues=(/0/), &
                            dstmaskvalues=(/0/), &
                            polemethod=ESMF_POLEMETHOD_NONE, &
                            srctermprocessing=isrctermprocessing, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            normtype=ESMF_NORMTYPE_FRACAREA, &
                            routehandle=regrid_land, &
                            regridmethod=method, &
                            unmappedDstList=unmapped_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

  bundle_nolandice_target = ESMF_FieldBundleCreate(name="land no landice target", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
      
 bundle_nolandice_input = ESMF_FieldBundleCreate(name="land no landice input", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleCreate", rc)
      
 call ESMF_FieldBundleAdd(bundle_nolandice_target, (/skin_temp_target_grid, terrain_from_input_grid,& 
                          soil_type_from_input_grid,soilm_tot_target_grid,soil_temp_target_grid/), rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
      
 call ESMF_FieldBundleAdd(bundle_nolandice_input, (/skin_temp_input_grid, terrain_input_grid,&
                          soil_type_input_grid,soilm_tot_input_grid,soil_temp_input_grid/), rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
 
 
 if (.not. sotyp_from_climo) then 
!   call ESMF_FieldBundleAdd(bundle_nolandice_target, (/soil_type_target_grid/), rc=rc)
!   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!      call error_handler("IN FieldBundleAdd", rc)
!   call ESMF_FieldBundleAdd(bundle_nolandice_input, (/soil_type_input_grid/), rc=rc)
!   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!      call error_handler("IN FieldBundleAdd", rc)
   print*,"- CALL Field_Regrid ."
   call ESMF_FieldRegrid(soil_type_input_grid, &
                         soil_type_target_grid, &
                         routehandle=regrid_land, &
                         zeroregion=ESMF_REGION_SELECT, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ,  rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)
    
   call ESMF_FieldGet(soil_type_target_grid, &
                    farrayPtr=soil_type_target_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)
    
   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)

   do ij = l(1), u(1)
     call ij_to_i_j(unmapped_ptr(ij), i_target, j_target, i, j)
     soil_type_target_ptr(i,j) = -9999.9 
   enddo
 !  call ESMF_FieldBundleGet(bundle_nolandice_target,fieldCount=num_fields,rc=rc)
 !  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
 !     call error_handler("IN FieldBundleGet", rc)
 !  sotyp_ind = 3
 endif
 
 if (.not. vgfrc_from_climo) then 
   call ESMF_FieldBundleAdd(bundle_nolandice_target, (/veg_greenness_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleAdd(bundle_nolandice_input, (/veg_greenness_input_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleGet(bundle_nolandice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)
   vgfrc_ind = num_fields
 endif
 
 if (.not. lai_from_climo) then 
   call ESMF_FieldBundleAdd(bundle_nolandice_target, (/lai_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleAdd(bundle_nolandice_input, (/lai_input_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleGet(bundle_nolandice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)
   lai_ind = num_fields
 endif
 
 if (.not. minmax_vgfrc_from_climo) then 
   call ESMF_FieldBundleAdd(bundle_nolandice_target, (/max_veg_greenness_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleAdd(bundle_nolandice_input, (/max_veg_greenness_input_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
      
   call ESMF_FieldBundleAdd(bundle_nolandice_target, (/min_veg_greenness_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
   call ESMF_FieldBundleAdd(bundle_nolandice_input, (/min_veg_greenness_input_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleAdd", rc)
      
   call ESMF_FieldBundleGet(bundle_nolandice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)
      
   mmvg_ind = num_fields-1
 endif
 
 call ESMF_FieldBundleGet(bundle_nolandice_target,fieldCount=num_fields,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleGet", rc)

 allocate(search_nums(num_fields))
 allocate(dozero(num_fields))
 
 search_nums(1:5) = (/85,7,224,85,86/)
 dozero(1:5) = (/.False.,.False.,.True.,.True.,.False./)
 
 !if (.not.sotyp_from_climo) then
 !  search_nums(sotyp_ind) = 226
 !  dozero(sotyp_ind) = .False.
 !endif
 
 if (.not. vgfrc_from_climo) then
   search_nums(vgfrc_ind) = 224
   dozero(vgfrc_ind) = .True.
 endif
 
 if (.not. lai_from_climo) then
   search_nums(lai_ind) = 229
   dozero(lai_ind) = .True.
 endif
 
 if (.not. minmax_vgfrc_from_climo) then
   search_nums(mmvg_ind) = 227
   dozero(mmvg_ind) = .True.
   
   search_nums(mmvg_ind+1) = 228
   dozero(mmvg_ind+1) = .True.
 endif

 call regrid_many(bundle_nolandice_input,bundle_nolandice_target,num_fields,regrid_land,dozero, &
                  unmapped_ptr=unmapped_ptr)
 deallocate(dozero)
 call ESMF_FieldBundleDestroy(bundle_nolandice_input,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 if (localpet == 0) then
   allocate (veg_type_target_one_tile(i_target,j_target))
 else
   allocate (veg_type_target_one_tile(0,0))
 endif

 do tile = 1, num_tiles_target_grid

   print*,"- CALL FieldGather FOR TARGET LANDMASK TILE: ", tile
   call ESMF_FieldGather(landmask_target_grid, mask_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   print*,"- CALL FieldGather FOR TARGET VEG TYPE TILE: ", tile
   call ESMF_FieldGather(veg_type_target_grid, veg_type_target_one_tile, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

!cfract Here mask_target_one_tile is 'landmask_target_grid'. it is
!cfract then modified by veg_type_target_grid. so this should work for 
!cfract grids.
   if (localpet == 0) then
     where(nint(veg_type_target_one_tile) == veg_type_landice_target) mask_target_one_tile = 0
   endif
   
   print*,"- CALL FieldGather FOR SOIL TYPE TARGET GRID, TILE: ", tile
   call ESMF_FieldGather(soil_type_target_grid, data_one_tile2, rootPet=0,tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGather", rc)
   if (localpet==0) then 
      call search_many(num_fields,bundle_nolandice_target,tile,search_nums,localpet, &
                soilt_climo=data_one_tile2, mask=mask_target_one_tile)
   else
      call search_many(num_fields,bundle_nolandice_target, tile,search_nums,localpet)
   endif
   
   print*,"- CALL FieldGather FOR TARGET GRID TOTAL SOIL MOISTURE, TILE: ", tile
   call ESMF_FieldGather(soilm_tot_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)

   if (localpet == 0) then
     do j = 1, lsoil_target
       data_one_tile = data_one_tile_3d(:,:,j)
       call search(data_one_tile, mask_target_one_tile, i_target, j_target, tile, 86)
       data_one_tile_3d(:,:,j) = data_one_tile
     enddo
   endif

   print*,"- CALL FieldGather FOR TARGET GRID SOIL TEMPERATURE, TILE: ", tile
   call ESMF_FieldGather(soil_temp_target_grid, data_one_tile_3d, rootPet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
      
   if (tg3_from_soil) then
     print*,"- CALL FieldScatter FOR TARGET GRID SUBSTRATE TEMPERATURE, TILE: ", tile
     call ESMF_FieldScatter(substrate_temp_target_grid, data_one_tile_3d(:,:,lsoil_target), rootPet=0, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)
   endif
   
   if (.not. sotyp_from_climo) then
     print*,"- CALL FieldGather FOR SOIL TYPE TARGET GRID LAND, TILE: ",tile
     call ESMF_FieldGather(soil_type_target_grid, data_one_tile,rootPet=0,tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)

     if (localpet == 0) then
       call search(data_one_tile, mask_target_one_tile, i_target, j_target,tile,226)
     endif

     print*,"- CALL FieldScatter FOR SOIL TYPE TARGET GRID, TILE: ", tile
     call ESMF_FieldScatter(soil_type_target_grid,data_one_tile,rootPet=0,tile=tile,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)
   endif

 enddo

 deallocate(search_nums) 
 call ESMF_FieldBundleDestroy(bundle_nolandice_target,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldBundleDestroy", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_land, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

 deallocate(veg_type_target_one_tile)

 deallocate(data_one_tile, data_one_tile2)
 deallocate(data_one_tile_3d)
 deallocate(mask_target_one_tile)

 return

 end subroutine interp
 
!> Compute liquid portion of the total soil moisture.
!!
!! @author George Gayno NOAA/EMC
 subroutine calc_liq_soil_moisture

 use esmf

 use model_grid, only                : landmask_target_grid

 use program_setup, only             : maxsmc_target, &
                                       bb_target, &
                                       satpsi_target

 use static_data, only               : soil_type_target_grid, &
                                       veg_type_target_grid

 implicit none
 
 integer                            :: clb(3), cub(3), rc
 integer                            :: i, j, n, soil_type

 integer(esmf_kind_i8), pointer     :: landmask_ptr(:,:)

 real                               :: bx, fk
 real(esmf_kind_r8), pointer        :: soilm_liq_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soilm_tot_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soil_temp_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soil_type_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_ptr(:,:)

 print*,"- COMPUTE LIQUID PORTION OF TOTAL SOIL MOISTURE."

 print*,"- CALL FieldGet FOR TOTAL SOIL MOISTURE."
 call ESMF_FieldGet(soilm_tot_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=soilm_tot_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LIQUID SOIL MOISTURE."
 call ESMF_FieldGet(soilm_liq_target_grid, &
                    farrayPtr=soilm_liq_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SOIL TEMPERATURE."
 call ESMF_FieldGet(soil_temp_target_grid, &
                    farrayPtr=soil_temp_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR VEGETATION TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SOIL TYPE."
 call ESMF_FieldGet(soil_type_target_grid, &
                    farrayPtr=soil_type_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LANDMASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    farrayPtr=landmask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)

!---------------------------------------------------------------------------------------------
! Check land points that are not permanent land ice.  
!---------------------------------------------------------------------------------------------

!cfract - 1 for at least some land.
     if (landmask_ptr(i,j) == 1 .and. nint(veg_type_ptr(i,j)) /= veg_type_landice_target) then

       soil_type = nint(soil_type_ptr(i,j))

       do n = clb(3), cub(3) 

         if (soil_temp_ptr(i,j,n) < (frz_h2o-0.0001)) then

           bx = bb_target(soil_type)

           if (bx .gt. blim) bx = blim

           fk=(((hlice/(grav*(-satpsi_target(soil_type))))*           &
            ((soil_temp_ptr(i,j,n)-frz_h2o)/soil_temp_ptr(i,j,n)))**             &
            (-1/bx))*maxsmc_target(soil_type)

           if (fk .lt. 0.02) fk = 0.02

           soilm_liq_ptr(i,j,n) = min ( fk, soilm_tot_ptr(i,j,n) )

!-----------------------------------------------------------------------
! now use iterative solution for liquid soil water content using
! FUNCTION FRH2O with the initial guess for SH2O from above explicit
! first guess.
!-----------------------------------------------------------------------

           soilm_liq_ptr(i,j,n) = frh2O(soil_temp_ptr(i,j,n),                        &
                           soilm_tot_ptr(i,j,n), soilm_liq_ptr(i,j,n),             &
                           maxsmc_target(soil_type),bb_target(soil_type),    &
                           satpsi_target(soil_type))

         else  ! temp above freezing. all moisture is liquid

           soilm_liq_ptr(i,j,n) = soilm_tot_ptr(i,j,n)

         end if  ! is soil layer below freezing?

       enddo ! soil layer

     end if ! is this point land?

   enddo
 enddo

 end subroutine calc_liq_soil_moisture

!> Calculate supercooled soil moisture
!!
!! Calculate amount of supercooled liquid soil water content if
!! temperature is below 273.15K. Requires Newton-type iteration to
!! solve the nonlinear implicit equation given in eqn 17 of Koren et. al
!! (1999, JGR, VOL 104(D16), 19569-19585).
!!
!! New version (June 2001): Much faster and more accurate Newton
!! iteration achieved by first taking log of eqn cited above -- less than
!! 4 (typically 1 or 2) iterations achieves convergence. Also, explicit
!! 1-step solution option for special case of parameter ck=0, which
!! reduces the original implicit equation to a simpler explicit form,
!! known as the "Flerchinger eqn". Improved handling of solution in the
!! limit of freezing point temperature.
!!
!! @param[in]  tkelv  Temperature (Kelvin)
!! @param[in]  smc    Total soil moisture content (volumetric)
!! @param[in]  sh2O   Liquid soil moisture content (volumetric)
!! @param[in]  smcmax  Saturation soil moisture content
!! @param[in]  bexp    Soil type "b" parameter
!! @param[in]  psis    Saturated soil matric potential
!! @return     frh2O   Supercooled liquid water content
!!
!! @author George Gayno NOAA/EMC @date 2005-05-20
 FUNCTION FRH2O (TKELV,SMC,SH2O,SMCMAX,BEXP,PSIS)

 use esmf

 IMPLICIT NONE

 INTEGER NLOG
 INTEGER KCOUNT

 REAL BEXP
 REAL BX
 REAL DENOM
 REAL DF
 REAL DSWL
 REAL FK
 REAL FRH2O
 REAL PSIS
 REAL(esmf_kind_r8) :: SH2O
 REAL(esmf_kind_r8) :: SMC
 REAL SMCMAX
 REAL SWL
 REAL SWLK
 REAL(esmf_kind_r8) :: TKELV

 REAL, PARAMETER                  :: CK    = 8.0
 REAL, PARAMETER                  :: ERROR = 0.005

! ----------------------------------------------------------------------
! LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
! SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS
! NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES.
! ----------------------------------------------------------------------

 BX = BEXP
 IF (BEXP .GT. BLIM) BX = BLIM

! ----------------------------------------------------------------------
! INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
! ----------------------------------------------------------------------

 NLOG=0
 KCOUNT=0

 IF (CK .NE. 0.0) THEN

! ----------------------------------------------------------------------
! OPTION 1: ITERATED SOLUTION FOR NONZERO CK
! IN KOREN ET AL, JGR, 1999, EQN 17
! ----------------------------------------------------------------------
! INITIAL GUESS FOR SWL (frozen content)
! ----------------------------------------------------------------------

   SWL = SMC-SH2O

! ----------------------------------------------------------------------
! KEEP WITHIN BOUNDS.
! ----------------------------------------------------------------------

   IF (SWL .GT. (SMC-0.02)) SWL = SMC-0.02
   IF (SWL .LT. 0.) SWL = 0.

! ----------------------------------------------------------------------
!  START OF ITERATIONS
! ----------------------------------------------------------------------

   DO WHILE ( (NLOG .LT. 10) .AND. (KCOUNT .EQ. 0) )

     NLOG = NLOG+1
     DF = LOG(( PSIS*GRAV/HLICE ) * ( ( 1.+CK*SWL )**2. ) *      &
        ( SMCMAX/(SMC-SWL) )**BX) - LOG(-(TKELV-frz_h2o)/TKELV)
     DENOM = 2. * CK / ( 1.+CK*SWL ) + BX / ( SMC - SWL )
     SWLK = SWL - DF/DENOM

! ----------------------------------------------------------------------
! BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
! ----------------------------------------------------------------------

     IF (SWLK .GT. (SMC-0.02)) SWLK = SMC - 0.02
     IF (SWLK .LT. 0.) SWLK = 0.

! ----------------------------------------------------------------------
! MATHEMATICAL SOLUTION BOUNDS APPLIED.
! ----------------------------------------------------------------------

     DSWL = ABS(SWLK-SWL)
     SWL = SWLK

! ----------------------------------------------------------------------
! IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
! WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
! ----------------------------------------------------------------------

     IF ( DSWL .LE. ERROR )  THEN
       KCOUNT = KCOUNT+1
     ENDIF

   END DO

! ----------------------------------------------------------------------
!  END OF ITERATIONS
! ----------------------------------------------------------------------
! BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
! ----------------------------------------------------------------------

   FRH2O = SMC - SWL

! ----------------------------------------------------------------------
! END OPTION 1
! ----------------------------------------------------------------------

 ENDIF

!-----------------------------------------------------------------------
! OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
! IN KOREN ET AL., JGR, 1999, EQN 17
! APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
! ----------------------------------------------------------------------

 IF (KCOUNT .EQ. 0) THEN

   FK = (((HLICE/(GRAV*(-PSIS)))*                  &
        ((TKELV-frz_h2o)/TKELV))**(-1/BX))*SMCMAX

   IF (FK .LT. 0.02) FK = 0.02

   FRH2O = MIN (FK, SMC)

 ENDIF

 RETURN

 END function frh2o

!> Adjust soil moisture for changes in soil type between the input and
!! target grids. Works for Noah land model only. Required to preserve
!! latent/sensible heat fluxes.
!!
!! @author George Gayno NOAA/EMC
 subroutine rescale_soil_moisture

 use esmf

 use model_grid, only                : landmask_target_grid

 use program_setup, only             : drysmc_input, drysmc_target, &
                                       maxsmc_input, maxsmc_target, &
                                       refsmc_input, refsmc_target, &
                                       wltsmc_input, wltsmc_target

 use static_data, only               : soil_type_target_grid, &
                                       veg_greenness_target_grid, &
                                       veg_type_target_grid

 implicit none

 integer                            :: clb(3), cub(3), i, j, k, rc
 integer                            :: soilt_input, soilt_target
 integer(esmf_kind_i8), pointer     :: landmask_ptr(:,:)

 real(esmf_kind_r8), pointer        :: soilm_tot_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soil_type_input_ptr(:,:)
 real(esmf_kind_r8), pointer        :: soil_type_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_greenness_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_ptr(:,:)
 real                               :: f1, fn, smcdir, smctra

 print*,"- RESCALE SOIL MOISTURE FOR CHANGES IN SOIL TYPE."

 print*,"- CALL FieldGet FOR TOTAL SOIL MOISTURE."
 call ESMF_FieldGet(soilm_tot_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=soilm_tot_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LAND MASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    farrayPtr=landmask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR VEGETATION TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR VEGETATION GREENNESS."
 call ESMF_FieldGet(veg_greenness_target_grid, &
                    farrayPtr=veg_greenness_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID SOIL TYPE."
 call ESMF_FieldGet(soil_type_target_grid, &
                    farrayPtr=soil_type_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SOIL TYPE FROM INPUT GRID."
 call ESMF_FieldGet(soil_type_from_input_grid, &
                    farrayPtr=soil_type_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
   do i = clb(1), cub(1)

!---------------------------------------------------------------------------------------------
! Check land points that are not permanent land ice.  
!---------------------------------------------------------------------------------------------

!cfract when '1' will contain some land.
     if (landmask_ptr(i,j) == 1 .and. nint(veg_type_ptr(i,j)) /= veg_type_landice_target) then

        soilt_target = nint(soil_type_target_ptr(i,j))
        soilt_input  = nint(soil_type_input_ptr(i,j))

!---------------------------------------------------------------------------------------------
! Rescale soil moisture at points where the soil type between the input and output
! grids is different.  Caution, this logic assumes the input and target grids use the same
! soil type dataset.
!---------------------------------------------------------------------------------------------

        if (soilt_target /= soilt_input) then
!---------------------------------------------------------------------------------------------
! Rescale top layer.  First, determine direct evaporation part:
!---------------------------------------------------------------------------------------------

          f1=(soilm_tot_ptr(i,j,1)-drysmc_input(soilt_input)) /    &
             (maxsmc_input(soilt_input)-drysmc_input(soilt_input))

          smcdir=drysmc_target(soilt_target) + f1 *        &
                (maxsmc_target(soilt_target) - drysmc_target(soilt_target))

!---------------------------------------------------------------------------------------------
! Continue top layer rescale.  Now determine transpiration part:
!---------------------------------------------------------------------------------------------

          if (soilm_tot_ptr(i,j,1) < refsmc_input(soilt_input)) then
            f1=(soilm_tot_ptr(i,j,1) - wltsmc_input(soilt_input)) /       &
               (refsmc_input(soilt_input) - wltsmc_input(soilt_input))
            smctra=wltsmc_target(soilt_target) + f1  *     &
                  (refsmc_target(soilt_target) - wltsmc_target(soilt_target))
          else
            f1=(soilm_tot_ptr(i,j,1) - refsmc_input(soilt_input)) /        &
               (maxsmc_input(soilt_input) - refsmc_input(soilt_input))
            smctra=refsmc_target(soilt_target) + f1 *      &
                  (maxsmc_target(soilt_target) - refsmc_target(soilt_target))
          endif

!---------------------------------------------------------------------------------------------
! Top layer is weighted by green vegetation fraction:
!---------------------------------------------------------------------------------------------

          soilm_tot_ptr(i,j,1) = ((1.0 - veg_greenness_ptr(i,j)) * smcdir)  + &
                                  (veg_greenness_ptr(i,j) * smctra)

!---------------------------------------------------------------------------------------------
! Rescale bottom layers as follows:
!
! - Rescale between wilting point and reference value when wilting < soil m < reference, or
! - Rescale between reference point and maximum value when reference < soil m < max.
!---------------------------------------------------------------------------------------------

          do k = 2, cub(3)
            if (soilm_tot_ptr(i,j,k) < refsmc_input(soilt_input)) then
              fn = (soilm_tot_ptr(i,j,k) - wltsmc_input(soilt_input)) /        &
                (refsmc_input(soilt_input) - wltsmc_input(soilt_input))
              soilm_tot_ptr(i,j,k) = wltsmc_target(soilt_target) + fn *         &
                (refsmc_target(soilt_target) - wltsmc_target(soilt_target))
            else
              fn = (soilm_tot_ptr(i,j,k) - refsmc_input(soilt_input)) /         &
                (maxsmc_input(soilt_input) - refsmc_input(soilt_input))
              soilm_tot_ptr(i,j,k) = refsmc_target(soilt_target) + fn *         &
                (maxsmc_target(soilt_target) - refsmc_target(soilt_target))
            endif
          enddo

        endif ! is soil type different?

!---------------------------------------------------------------------------------------------
! Range check all layers.
!---------------------------------------------------------------------------------------------

        soilm_tot_ptr(i,j,1)=min(soilm_tot_ptr(i,j,1),maxsmc_target(soilt_target))
        soilm_tot_ptr(i,j,1)=max(drysmc_target(soilt_target),soilm_tot_ptr(i,j,1))

        do k = 2, cub(3)
          soilm_tot_ptr(i,j,k)=min(soilm_tot_ptr(i,j,k),maxsmc_target(soilt_target))
          soilm_tot_ptr(i,j,k)=max(wltsmc_target(soilt_target),soilm_tot_ptr(i,j,k))
        enddo

     endif ! is this a land point?

   enddo
 enddo

 return

 end subroutine rescale_soil_moisture

!> Adjust soil temperature for changes in terrain height between the input and
!! target grids.
!!
!! @author George Gayno NOAA/EMC
 subroutine adjust_soilt_for_terrain

 use model_grid, only                : landmask_target_grid,  &
                                       terrain_target_grid

 use static_data, only               : veg_type_target_grid

 implicit none

 integer                            :: clb(3), cub(3), i, j, k, rc
 integer(esmf_kind_i8), pointer     :: landmask_ptr(:,:)

 real, parameter                    :: lapse_rate  = 6.5e-03
 real                               :: terrain_diff
 real(esmf_kind_r8), pointer        :: terrain_input_ptr(:,:)
 real(esmf_kind_r8), pointer        :: terrain_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_target_ptr(:,:)
 real(esmf_kind_r8), pointer        :: soil_temp_target_ptr(:,:,:)

 print*,"- CALL FieldGet FOR TARGET GRID LAND-SEA MASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    farrayPtr=landmask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID VEGETATION TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID TERRAIN."
 call ESMF_FieldGet(terrain_target_grid, &
                    farrayPtr=terrain_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TERRAIN INTERP TO TARGET GRID."
 call ESMF_FieldGet(terrain_from_input_grid, &
                    farrayPtr=terrain_input_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SOIL TEMP TARGET GRID."
 call ESMF_FieldGet(soil_temp_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=soil_temp_target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
 
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
!cfract at least some land when equal to '1'.
   if (landmask_ptr(i,j) == 1) then
     terrain_diff = abs(terrain_input_ptr(i,j) - terrain_target_ptr(i,j))
     if (terrain_diff > 100.0) then
       do k = clb(3), cub(3)
         soil_temp_target_ptr(i,j,k) = soil_temp_target_ptr(i,j,k) + &
              ((terrain_input_ptr(i,j) - terrain_target_ptr(i,j)) * lapse_rate)
         if (nint(veg_type_target_ptr(i,j)) == veg_type_landice_target) then
           soil_temp_target_ptr(i,j,k) = min(soil_temp_target_ptr(i,j,k), 273.16)
         endif
       enddo
     endif
   endif
 enddo
 enddo

 end subroutine adjust_soilt_for_terrain

!> Adjust soil levels of the input grid if there is a mismatch between input and
!! target grids. Presently can only convert from 9 to 4 levels. 
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames
!! @author Jeff Beck
 subroutine adjust_soil_levels(localpet)
 use model_grid, only       : lsoil_target, i_input, j_input, input_grid
 use sfc_input_data, only   : lsoil_input, soil_temp_input_grid, &
                              soilm_liq_input_grid, soilm_tot_input_grid
 implicit none
 integer, intent(in)                   :: localpet
 character(len=500)       :: msg
 character(len=2)         :: lsoil_input_ch, lsoil_target_ch
 integer                  :: rc
 real(esmf_kind_r8)          :: tmp(i_input,j_input), &
                                data_one_tile(i_input,j_input,lsoil_input), &
                                tmp3d(i_input,j_input,lsoil_target)
 if (lsoil_input == 9 .and. lsoil_target == 4) then
   print*, "CONVERTING FROM 9 INPUT SOIL LEVELS TO 4 TARGET SOIL LEVELS"
   call ESMF_FieldGather(soil_temp_input_grid, data_one_tile, rootPet=0, tile=1, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
      
   call ESMF_FieldDestroy(soil_temp_input_grid,rc=rc)
   soil_temp_input_grid = ESMF_FieldCreate(input_grid, &
                         typekind=ESMF_TYPEKIND_R8, &
                         staggerloc=ESMF_STAGGERLOC_CENTER, &
                         ungriddedLBound=(/1/), &
                         ungriddedUBound=(/lsoil_target/), rc=rc)
                                         
   if(localpet==0)then
      tmp3d(:,:,1)= (data_one_tile(:,:,1) + data_one_tile(:,:,2))/2.0 * 0.1 + &
                                      (data_one_tile(:,:,2) + data_one_tile(:,:,3))/2.0 * 0.3 + &
                                      (data_one_tile(:,:,3) + data_one_tile(:,:,4))/2.0 * 0.6
      tmp = (data_one_tile(:,:,6) - data_one_tile(:,:,5)) / 30.0 * 10.0 + data_one_tile(:,:,5) !Linear approx. of 40 cm obs
      tmp3d(:,:,2)= (data_one_tile(:,:,4) + data_one_tile(:,:,5)) / 2.0 * 0.75 + &
                                      (data_one_tile(:,:,5) + tmp) / 2.0 * 0.25
      tmp3d(:,:,3)= (tmp + data_one_tile(:,:,6)) /2.0 * (1.0/3.0) + &
                                      (data_one_tile(:,:,6) + data_one_tile(:,:,7)) / 2.0 * (2.0/3.0)
      tmp = (data_one_tile(:,:,9) - data_one_tile(:,:,9)) / 140.0 * 40.0 + data_one_tile(:,:,8) !Linear approx of 200 cm obs
      tmp3d(:,:,4)= (data_one_tile(:,:,7) + data_one_tile(:,:,8)) / 2.0 * 0.6 + &
                                      (data_one_tile(:,:,8) + tmp) / 2.0 * 0.4
   endif
  
   call ESMF_FieldScatter(soil_temp_input_grid, tmp3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)   
                                                                              
   call ESMF_FieldGather(soilm_tot_input_grid, data_one_tile, rootPet=0, tile=1, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
      
   call ESMF_FieldDestroy(soilm_tot_input_grid,rc=rc)
   soilm_tot_input_grid = ESMF_FieldCreate(input_grid, &
                         typekind=ESMF_TYPEKIND_R8, &
                         staggerloc=ESMF_STAGGERLOC_CENTER, &
                         ungriddedLBound=(/1/), &
                         ungriddedUBound=(/lsoil_target/), rc=rc)
                                         
  if(localpet==0) then
      tmp3d(:,:,1)= (data_one_tile(:,:,1) + data_one_tile(:,:,2))/2.0 * 0.1 + &
                                      (data_one_tile(:,:,2) + data_one_tile(:,:,3))/2.0 * 0.3 + &
                                      (data_one_tile(:,:,3) + data_one_tile(:,:,4))/2.0 * 0.6
      tmp = (data_one_tile(:,:,6) - data_one_tile(:,:,5)) / 30.0 * 10.0 + data_one_tile(:,:,5) !Linear approx. of 40 cm obs
      tmp3d(:,:,2)= (data_one_tile(:,:,4) + data_one_tile(:,:,5)) / 2.0 * 0.75 + &
                                      (data_one_tile(:,:,5) + tmp) / 2.0 * 0.25
      tmp3d(:,:,3)= (tmp + data_one_tile(:,:,6)) /2.0 * (1.0/3.0) + &
                                      (data_one_tile(:,:,6) + data_one_tile(:,:,7)) / 2.0 * (2.0/3.0)
      tmp = (data_one_tile(:,:,9) - data_one_tile(:,:,9)) / 140.0 * 40.0 + data_one_tile(:,:,8) !Linear approx of 200 cm obs
      tmp3d(:,:,4)= (data_one_tile(:,:,7) + data_one_tile(:,:,8)) / 2.0 * 0.6 + &
                                      (data_one_tile(:,:,8) + tmp) / 2.0 * 0.4
   endif
  
   call ESMF_FieldScatter(soilm_tot_input_grid, tmp3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)   
  
   call ESMF_FieldGather(soilm_liq_input_grid, data_one_tile, rootPet=0, tile=1, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", rc)
      
   call ESMF_FieldDestroy(soilm_liq_input_grid,rc=rc)
   soilm_liq_input_grid = ESMF_FieldCreate(input_grid, &
                         typekind=ESMF_TYPEKIND_R8, &
                         staggerloc=ESMF_STAGGERLOC_CENTER, &
                         ungriddedLBound=(/1/), &
                         ungriddedUBound=(/lsoil_target/), rc=rc)
  if(localpet==0) then
      tmp3d(:,:,1)= (data_one_tile(:,:,1) + data_one_tile(:,:,2))/2.0 * 0.1 + &
                                      (data_one_tile(:,:,2) + data_one_tile(:,:,3))/2.0 * 0.3 + &
                                      (data_one_tile(:,:,3) + data_one_tile(:,:,4))/2.0 * 0.6
      tmp = (data_one_tile(:,:,6) - data_one_tile(:,:,5)) / 30.0 * 10.0 + data_one_tile(:,:,5) !Linear approx. of 40 cm obs
      tmp3d(:,:,2)= (data_one_tile(:,:,4) + data_one_tile(:,:,5)) / 2.0 * 0.75 + &
                                      (data_one_tile(:,:,5) + tmp) / 2.0 * 0.25
      tmp3d(:,:,3)= (tmp + data_one_tile(:,:,6)) /2.0 * (1.0/3.0) + &
                                      (data_one_tile(:,:,6) + data_one_tile(:,:,7)) / 2.0 * (2.0/3.0)
      tmp = (data_one_tile(:,:,9) - data_one_tile(:,:,9)) / 140.0 * 40.0 + data_one_tile(:,:,8) !Linear approx of 200 cm obs
      tmp3d(:,:,4)= (data_one_tile(:,:,7) + data_one_tile(:,:,8)) / 2.0 * 0.6 + &
                                      (data_one_tile(:,:,8) + tmp) / 2.0 * 0.4
   endif
  
   call ESMF_FieldScatter(soilm_liq_input_grid, tmp3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)   
 
 elseif (lsoil_input /= lsoil_target) then
  rc = -1
  write(lsoil_input_ch, '(i2)') lsoil_input
  write(lsoil_target_ch, '(i2)') lsoil_target
  msg="NUMBER OF SOIL LEVELS IN INPUT " // lsoil_input_ch // " AND OUTPUT " &
      // lsoil_target_ch // " MUST EITHER BE EQUAL OR 9 AND 4 RESPECTIVELY."
  call error_handler(msg, rc)
 endif
 
 end subroutine adjust_soil_levels

!> Set roughness length at land and sea ice. At land, roughness is
!! set from a lookup table based on the vegetation type. At sea ice,
!! roughness is set to 1 cm.
!!
!! @author George Gayno NOAA/EMC
 subroutine roughness

 use model_grid, only                : landmask_target_grid, &
                                       seamask_target_grid
 use static_data, only               : veg_type_target_grid
 use program_setup, only             : fract_grid

 implicit none

 integer                            :: clb(2), cub(2), i, j, rc
 integer(esmf_kind_i8), pointer     :: landmask_ptr(:,:)
 integer(esmf_kind_i8), pointer     :: seamask_ptr(:,:)

 real                               :: z0_igbp(20)
 real(esmf_kind_r8), pointer        :: data_ptr(:,:)
 real(esmf_kind_r8), pointer        :: data_ptr2(:,:)
 real(esmf_kind_r8), pointer        :: data_ptr3(:,:)
 real(esmf_kind_r8), pointer        :: fice_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_ptr(:,:)

 data z0_igbp /1.089, 2.653, 0.854, 0.826, 0.800, 0.050,  &
               0.030, 0.856, 0.856, 0.150, 0.040, 0.130,  &
               1.000, 0.250, 0.011, 0.011, 0.001, 0.076,  &
               0.050, 0.030/

 print*,"- CALL FieldGet FOR TARGET GRID LAND-SEA MASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=landmask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID SEA ICE."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=fice_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID VEGETATION TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID Z0."
 call ESMF_FieldGet(z0_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if(fract_grid)then

   print*,"- CALL FieldGet FOR TARGET GRID Z0 ICE."
   call ESMF_FieldGet(z0_ice_target_grid, &
                      farrayPtr=data_ptr2, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) > 0.0) then
       data_ptr2(i,j) = 1.0
     else
       data_ptr2(i,j) = -1.e20
     endif
   enddo
   enddo

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 1) then
       data_ptr(i,j) = z0_igbp(nint(veg_type_ptr(i,j))) * 100.0
     else
       data_ptr(i,j) = -1.e20
     endif
   enddo
   enddo

   print*,"- CALL FieldGet FOR TARGET GRID Z0 WATER."
   call ESMF_FieldGet(z0_water_target_grid, &
                      farrayPtr=data_ptr3, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

   print*,"- CALL FieldGet FOR TARGET SEA MASK."
   call ESMF_FieldGet(seamask_target_grid, &
                      farrayPtr=seamask_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) == 1.0_esmf_kind_r8 .or. seamask_ptr(i,j) == 0) then
       data_ptr3(i,j) = -1.e20
     endif
   enddo
   enddo

 else ! non-fractional grid

!cfract should check for fice instead? under fractional
!cfract grids need to preserve original landmask_target.
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) > 0.0) then
       data_ptr(i,j) = 1.0
     elseif (landmask_ptr(i,j) == 1) then
       data_ptr(i,j) = z0_igbp(nint(veg_type_ptr(i,j))) * 100.0
     endif
   enddo
   enddo

 endif

 end subroutine roughness

!> Perform some quality control checks before output.
!!
!! @author George Gayno NOAA/EMC
 subroutine qc_check

 use program_setup, only             : fract_grid

 use model_grid, only                : landmask_target_grid, &
                                       seamask_target_grid

 use static_data, only               : alvsf_target_grid, &
                                       alvwf_target_grid, &
                                       alnsf_target_grid, &
                                       alnwf_target_grid, &
                                       facsf_target_grid, &
                                       facwf_target_grid, &
                                       mxsno_albedo_target_grid, &
                                       max_veg_greenness_target_grid, &
                                       min_veg_greenness_target_grid, &
                                       slope_type_target_grid, &
                                       soil_type_target_grid, &
                                       substrate_temp_target_grid, &
                                       veg_greenness_target_grid, &
                                       veg_type_target_grid

 implicit none

 integer                            :: clb(2), cub(2), i, j, rc
 integer(esmf_kind_i8), pointer     :: landmask_ptr(:,:)
 integer(esmf_kind_i8), pointer     :: seamask_ptr(:,:)

 real(esmf_kind_r8), pointer        :: data_ptr(:,:)
 real(esmf_kind_r8), pointer        :: data3d_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: ice_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soilmt_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: soilml_ptr(:,:,:)
 real(esmf_kind_r8), pointer        :: veg_greenness_ptr(:,:)
 real(esmf_kind_r8), pointer        :: veg_type_ptr(:,:)
 real(esmf_kind_r8), pointer        :: seaice_skint_ptr(:,:)
 real(esmf_kind_r8), pointer        :: skint_ptr(:,:)
 real(esmf_kind_r8), pointer        :: fice_ptr(:,:)
 real(esmf_kind_r8), pointer        :: hice_ptr(:,:)
 real(esmf_kind_r8), pointer        :: tg3_ptr(:,:)
 real(esmf_kind_r8), pointer        :: tg3_ice_ptr(:,:)
 real(esmf_kind_r8), pointer        :: snod_ptr(:,:)
 real(esmf_kind_r8), pointer        :: snol_ptr(:,:)

 print*,"- CALL FieldGet FOR TARGET GRID LAND-SEA MASK."
 call ESMF_FieldGet(landmask_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=landmask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (fract_grid) then
   print*,"- CALL FieldGet FOR TARGET GRID SEA MASK."
   call ESMF_FieldGet(seamask_target_grid, &
                      farrayPtr=seamask_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
 endif

 print*,"- CALL FieldGet FOR TARGET GRID SEA ICE FRACTION."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=fice_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- SET NON-LAND FLAG FOR TARGET GRID SLOPE TYPE."
 call ESMF_FieldGet(slope_type_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract Setting to flag value at water. With fractional grid,
!cfract restrict this to points that are all water. i.e.,
!cfract where landmask_ptr = 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID SOIL TYPE."
 call ESMF_FieldGet(soil_type_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract Setting to flag value at water. With fractional grid,
!cfract restrict this to points that are all water. i.e.,
!cfract where landmask_ptr = 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID VEGETATION TYPE."
 call ESMF_FieldGet(veg_type_target_grid, &
                    farrayPtr=veg_type_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract Setting to flag value at water. With fractional grid,
!cfract restrict this to points that are all water. i.e.,
!cfract where landmask_ptr = 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) veg_type_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET TARGET GRID ALVSF AT NON-LAND."
 call ESMF_FieldGet(alvsf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract This array contains albedo at points with at least
!cfract some land. At all other points set to flag value.
!cfract That means points where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.06 ! gfs physics flag value
   if (fract_grid .and. landmask_ptr(i,j) /= 1) data_ptr(i,j) = -1.e20 ! gfs physics flag value
 enddo
 enddo

 if(fract_grid)then
   print*,"- SET TARGET GRID ALVSF_NL AT NON-LAND."
   call ESMF_FieldGet(alvsf_nl_target_grid, &
                      farrayPtr=data_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (seamask_ptr(i,j) == 1) data_ptr(i,j) = 0.06 ! gfs flag value at any non-land
   enddo
   enddo
 endif

 print*,"- SET TARGET GRID ALVWF AT NON-LAND."
 call ESMF_FieldGet(alvwf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.06 ! gfs physics flag value
   if (fract_grid .and. landmask_ptr(i,j) /= 1) data_ptr(i,j) = -1.e20 ! gfs physics flag value
 enddo
 enddo

 if(fract_grid)then
   print*,"- SET TARGET GRID ALVWF_NL AT NON-LAND."
   call ESMF_FieldGet(alvwf_nl_target_grid, &
                      farrayPtr=data_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (seamask_ptr(i,j) == 1) data_ptr(i,j) = 0.06 ! gfs flag value at any non-land
   enddo
   enddo
 endif

 print*,"- SET TARGET GRID ALNSF AT NON-LAND."
 call ESMF_FieldGet(alnsf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.06 ! gfs physics flag value
   if (fract_grid .and. landmask_ptr(i,j) /= 1) data_ptr(i,j) = -1.e20 ! gfs physics flag value
 enddo
 enddo

 if(fract_grid)then
   print*,"- SET TARGET GRID ALNSF_NL AT NON-LAND."
   call ESMF_FieldGet(alnsf_nl_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (seamask_ptr(i,j) == 1) data_ptr(i,j) = 0.06 ! gfs flag value at any non-land
   enddo
   enddo
 endif

 print*,"- SET TARGET GRID ALNWF AT NON-LAND."
 call ESMF_FieldGet(alnwf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.06 ! gfs physics flag value
   if (fract_grid .and. landmask_ptr(i,j) /= 1) data_ptr(i,j) = -1.e20 ! gfs physics flag value
 enddo
 enddo

 if(fract_grid)then
   print*,"- SET TARGET GRID ALNWF_NL AT NON-LAND."
   call ESMF_FieldGet(alnwf_nl_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", rc)
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (seamask_ptr(i,j) == 1) data_ptr(i,j) = 0.06  ! gfs flag value at any non-land.
   enddo
   enddo
 endif

 print*,"- SET NON-LAND FLAG FOR TARGET GRID FACSF."
 call ESMF_FieldGet(facsf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID FACWF."
 call ESMF_FieldGet(facwf_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID MAXIMUM GREENNESS."
 call ESMF_FieldGet(max_veg_greenness_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID MINIMUM GREENNESS."
 call ESMF_FieldGet(min_veg_greenness_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID VEGETATION GREENNESS."
 call ESMF_FieldGet(veg_greenness_target_grid, &
                    farrayPtr=veg_greenness_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) veg_greenness_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- SET NON-LAND FLAG FOR TARGET GRID MAX SNOW ALBEDO."
 call ESMF_FieldGet(mxsno_albedo_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!cfract for fract grid, this is where landmask_ptr is 0.
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) /= 1) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- ZERO OUT TARGET GRID CANOPY MOISTURE CONTENT WHERE NO PLANTS."
 call ESMF_FieldGet(canopy_mc_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (veg_greenness_ptr(i,j) <= 0.01) data_ptr(i,j) = 0.0
 enddo
 enddo

 print*,"- CALL FieldGet FOR TARGET GRID ICE SKIN TEMP."
 call ESMF_FieldGet(seaice_skin_temp_target_grid, &
                    farrayPtr=seaice_skint_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- SET TARGET GRID SKIN TEMP AT ICE POINTS."
 call ESMF_FieldGet(skin_temp_target_grid, &
                    farrayPtr=skint_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- SET TARGET GRID SEA ICE DEPTH TO ZERO AT NON-ICE POINTS."
 call ESMF_FieldGet(seaice_depth_target_grid, &
                    farrayPtr=hice_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if(fract_grid)then
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (fice_ptr(i,j) > 0.0) then
!    skint_ptr(i,j) = (fice_ptr(i,j) * seaice_skint_ptr(i,j)) +  &
!                     ( (1.0 - fice_ptr(i,j)) * frz_ice )
   else
     seaice_skint_ptr(i,j) = -1.e20
     hice_ptr(i,j) = 0.0
   endif
 enddo
 enddo
 else
 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (fice_ptr(i,j) > 0.0) then
     skint_ptr(i,j) = (fice_ptr(i,j) * seaice_skint_ptr(i,j)) +  &
                      ( (1.0 - fice_ptr(i,j)) * frz_ice )
   else
     seaice_skint_ptr(i,j) = skint_ptr(i,j)
     hice_ptr(i,j) = 0.0
   endif
 enddo
 enddo
 endif

 if (fract_grid) then

   print*,"- SET TARGET GRID SST FLAG VALUE."

   call ESMF_FieldGet(sst_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) == 1.0_esmf_kind_r8 .or. seamask_ptr(i,j) == 0.0) then
       data_ptr(i,j) = -1.e20
     endif
   enddo
   enddo

 endif

 print*,"- SET TARGET GRID SUBSTRATE TEMP AT ICE."

 if (fract_grid) then

   call ESMF_FieldGet(seaice_substrate_temp_target_grid, &
                    farrayPtr=tg3_ice_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) > 0.0) then  ! sea ice
       tg3_ice_ptr(i,j) = frz_ice
     else
       tg3_ice_ptr(i,j) = -1.e20
     endif
   enddo
   enddo

   call ESMF_FieldGet(substrate_temp_target_grid, &
                    farrayPtr=tg3_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 0.0) then  ! sea ice
       tg3_ptr(i,j) = -1.e20
     endif
   enddo
   enddo

 else

 call ESMF_FieldGet(substrate_temp_target_grid, &
                    farrayPtr=tg3_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (fice_ptr(i,j) > 0.0) then  ! sea ice
     tg3_ptr(i,j) = frz_ice
   elseif (landmask_ptr(i,j) == 0) then  ! open water flag value.
     tg3_ptr(i,j) = skint_ptr(i,j)
   endif
 enddo
 enddo

 endif

 if (fract_grid) then

   print*,"- SET MISSING FLAG AT TARGET GRID SNOW DEPTH AT ICE."
   call ESMF_FieldGet(snow_depth_at_ice_target_grid, &
                    farrayPtr=snod_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", rc)

   print*,"- SET MISSING FLAG AT TARGET GRID SNOW LIQ EQUIV AT ICE."
   call ESMF_FieldGet(snow_liq_equiv_at_ice_target_grid, &
                    farrayPtr=snol_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) == 0.0) then
       snol_ptr(i,j) = -1.e20
       snod_ptr(i,j) = -1.e20
     end if
   enddo
   enddo

 endif

 print*,"- ZERO OUT TARGET GRID SNOW DEPTH AT OPEN WATER."
 call ESMF_FieldGet(snow_depth_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (fract_grid) then
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) /= 1) then  ! not land
       data_ptr(i,j) = -1.e20
     end if
   enddo
   enddo
 else
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 0 .and. fice_ptr(i,j) == 0.0) then
       data_ptr(i,j) = 0.0
     end if
   enddo
   enddo
 endif

 print*,"- ZERO OUT TARGET GRID SNOW LIQ AT OPEN WATER."
 call ESMF_FieldGet(snow_liq_equiv_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (fract_grid) then
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) /= 1) then  ! not land
       data_ptr(i,j) = -1.e20
     end if
   enddo
   enddo
 else
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 0 .and. fice_ptr(i,j) == 0.0) then
       data_ptr(i,j) = 0.0
     end if
   enddo
   enddo
 endif

 print*,"- SET NON-LAND FLAG VALUE FOR TARGET GRID TOTAL SOIL MOISTURE."
 call ESMF_FieldGet(soilm_tot_target_grid, &
                    farrayPtr=soilmt_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- SET NON-LAND FLAG VALUE FOR  TARGET GRID LIQUID SOIL MOISTURE."
 call ESMF_FieldGet(soilm_liq_target_grid, &
                    farrayPtr=soilml_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (fract_grid) then
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 0 .or. &
         nint(veg_type_ptr(i,j)) == veg_type_landice_target) then
       soilmt_ptr(i,j,:) = 1.0
       soilml_ptr(i,j,:) = 1.0
     endif
   enddo
   enddo
 else
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (fice_ptr(i,j) > 0.0 .or. landmask_ptr(i,j) == 0 .or. &
         nint(veg_type_ptr(i,j)) == veg_type_landice_target) then
       soilmt_ptr(i,j,:) = 1.0
       soilml_ptr(i,j,:) = 1.0
     endif
   enddo
   enddo
 endif

 print*,"- SET OPEN WATER FLAG FOR TARGET GRID SOIL TEMPERATURE."
 call ESMF_FieldGet(soil_temp_target_grid, &
                    farrayPtr=data3d_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldGet", rc)

 if (fract_grid) then

  do j = clb(2), cub(2)
  do i = clb(1), cub(1)
    if (landmask_ptr(i,j) == 0) then
      data3d_ptr(i,j,:) = -1.e20
    endif
  enddo
  enddo

  print*,"- SET FLAG FOR TARGET GRID ICE TEMPERATURE."
  call ESMF_FieldGet(ice_temp_target_grid, &
                     farrayPtr=ice_ptr, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

  do j = clb(2), cub(2)
  do i = clb(1), cub(1)
    if (fice_ptr(i,j) == 0.0) then
      ice_ptr(i,j,:) = -1.e20
    endif
  enddo
  enddo

 else

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)
   if (landmask_ptr(i,j) == 0 .and. fice_ptr(i,j) == 0.0) then
     data3d_ptr(i,j,:) = skint_ptr(i,j)  ! open water flag value.
   endif
 enddo
 enddo

 endif

 if (fract_grid) then ! set flag value at non-land
   do j = clb(2), cub(2)
   do i = clb(1), cub(1)
     if (landmask_ptr(i,j) == 0) then
       skint_ptr(i,j) = -1.e20
     endif
   enddo
   enddo
 endif

 return

 end subroutine qc_check

!> nst is not active at land or sea ice points.  Set nst fields to flag values at these
!! points.
!!
!! @author George Gayno NOAA/EMC
 subroutine nst_land_fill

 use model_grid, only         : seamask_target_grid

 use program_setup, only      : fract_grid

 implicit none

 integer(esmf_kind_i8), pointer     :: mask_ptr(:,:)
 integer                            :: rc,i
 integer, PARAMETER                 :: num_nst_fields_minus2 = 16
 integer, PARAMETER                 :: xz_fill = 30.0
 integer, PARAMETER                 :: nst_fill = 0.0

 real(esmf_kind_r8), pointer        :: data_ptr(:,:)
 real(esmf_kind_r8), pointer        :: fice_ptr(:,:)
 real(esmf_kind_r8), pointer        :: skint_ptr(:,:)

 type(esmf_field)                   :: temp_field
 type(esmf_fieldbundle)             :: nst_bundle

 print*,"- CALL FieldGet FOR TARGET GRID SEAMASK."
 call ESMF_FieldGet(seamask_target_grid, &
                    farrayPtr=mask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TARGET GRID SEAICE FRACT."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=fice_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", rc)
    
 nst_bundle = ESMF_FieldBundleCreate(name="nst_bundle", rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleCreate", rc)

 call ESMF_FieldBundleAdd(nst_bundle, (/c_d_target_grid,c_0_target_grid,d_conv_target_grid, &
                          dt_cool_target_grid,ifd_target_grid,qrain_target_grid,&
                          w_d_target_grid,w_0_target_grid,xs_target_grid,xt_target_grid,&
                          xu_target_grid,xv_target_grid,xtts_target_grid,xzts_target_grid, &
                          z_c_target_grid, zm_target_grid/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleAdd", rc) 

 print*,"- CALL FieldGet FOR TREF."
 call ESMF_FieldGet(tref_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SKIN T."
 call ESMF_FieldGet(skin_temp_target_grid, &
                    farrayPtr=skint_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", rc)

!cfract Setting filler value for tref at ice and land points. 
!cfract Under fractional grids skin t is not currently defined.
!cfract So, set to ice temp.

 where(mask_ptr == 0) data_ptr = skint_ptr
 if (fract_grid) then
   where(fice_ptr > 0.0) data_ptr = frz_ice
 else
   where(fice_ptr > 0.0) data_ptr = skint_ptr
 endif

! xz

 print*,"- CALL FieldGet FOR XZ."
 call ESMF_FieldGet(xz_target_grid, &
                    farrayPtr=data_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", rc)

!cfract same as above.
 where(mask_ptr == 0) data_ptr = xz_fill
 where(fice_ptr > 0.0) data_ptr = xz_fill

 do i = 1,num_nst_fields_minus2
   
   call ESMF_FieldBundleGet(nst_bundle,i,temp_field,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN FieldBundleGet", rc)
       
   call ESMF_FieldGet(temp_field,farrayPtr=data_ptr,rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
     call error_handler("IN FieldGet", rc)
     
!cfract same as above.
   where(mask_ptr == 0) data_ptr = nst_fill
   where(fice_ptr > 0.0) data_ptr = nst_fill

 enddo

 call ESMF_FieldBundleDestroy(nst_bundle,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldBundleDestroy", rc)  
      
 end subroutine nst_land_fill

!> Create ESMF fields for the target grid surface variables
!!
!! @author George Gayno NOAA/EMC
 subroutine create_surface_esmf_fields

 use model_grid, only         : target_grid, lsoil_target

 use program_setup, only      : fract_grid

 implicit none

 integer                        :: rc

 real(esmf_kind_r8), pointer    :: target_ptr(:,:), target_ptr_3d(:,:,:)
 real                           :: init_val = -999.9

 print*,"- CALL FieldCreate FOR TARGET GRID T2M."
 t2m_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="t2m_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid t2m."
 call ESMF_FieldGet(t2m_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID Q2M."
 q2m_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="q2m_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid q2m."
 call ESMF_FieldGet(q2m_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID TPRCP."
 tprcp_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="tprcp_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid tprcp."
 call ESMF_FieldGet(tprcp_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID F10M."
 f10m_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="f10m_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid f10m."
 call ESMF_FieldGet(f10m_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID FFMM."
 ffmm_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="ffmm_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid ffmm."
 call ESMF_FieldGet(ffmm_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID USTAR."
 ustar_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name="ustar_target_grid", &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid ustar."
 call ESMF_FieldGet(ustar_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID SNOW LIQ EQUIV."
 snow_liq_equiv_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="snow_liq_equiv_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid snow liq equiv."
 call ESMF_FieldGet(snow_liq_equiv_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if (fract_grid) then
 print*,"- CALL FieldCreate FOR TARGET GRID SNOW LIQ EQUIV AT SEA ICE."
 snow_liq_equiv_at_ice_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="snow_liq_equiv_at_ice_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid snow liq equiv at sea ice."
 call ESMF_FieldGet(snow_liq_equiv_at_ice_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 endif ! fractional grid

 print*,"- CALL FieldCreate FOR TARGET GRID SNOW DEPTH."
 snow_depth_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="snow_depth_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid snow depth."
 call ESMF_FieldGet(snow_depth_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if (fract_grid) then
 print*,"- CALL FieldCreate FOR TARGET GRID SNOW DEPTH AT SEA ICE."
 snow_depth_at_ice_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="snow_depth_at_ice_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid snow depth at sea ice."
 call ESMF_FieldGet(snow_depth_at_ice_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val
 endif

 print*,"- CALL FieldCreate FOR TARGET GRID SEA ICE FRACTION."
 seaice_fract_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="seaice_fract_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid sea ice fraction."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID SEA ICE DEPTH."
 seaice_depth_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="seaice_depth_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET sea ice depth."
 call ESMF_FieldGet(seaice_depth_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if(fract_grid)then
 print*,"- CALL FieldCreate FOR TARGET GRID sst."
 sst_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="sst_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET sst."
 call ESMF_FieldGet(sst_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val
 endif

 print*,"- CALL FieldCreate FOR TARGET GRID SEA ICE SKIN TEMP."
 seaice_skin_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="seaice_skin_temp_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET sea ice skin temp."
 call ESMF_FieldGet(seaice_skin_temp_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if (fract_grid) then

   print*,"- CALL FieldCreate FOR TARGET GRID SEA ICE SUBSTRATE TEMP."
   seaice_substrate_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="seaice_substrate_temp_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldCreate", rc)

   print*,"- INITIALIZE TARGET sea ice substrate temp."
   call ESMF_FieldGet(seaice_substrate_temp_target_grid, &
                      farrayPtr=target_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

   target_ptr = init_val

 endif

 print*,"- CALL FieldCreate FOR TARGET GRID SRFLAG."
 srflag_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="srflag_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET srflag."
 call ESMF_FieldGet(srflag_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID SKIN TEMPERATURE."
 skin_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="skin_temp_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid skin temp."
 call ESMF_FieldGet(skin_temp_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if(fract_grid)then
 print*,"- CALL FieldCreate FOR TARGET ALVSF AT NON-LAND."
 alvsf_nl_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="alvsf_nl_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET ALVSF AT NON-LAND."
 call ESMF_FieldGet(alvsf_nl_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = -1.e20

 print*,"- CALL FieldCreate FOR TARGET ALVWF AT NON-LAND."
 alvwf_nl_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="alvwf_nl_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET ALVWF AT NON-LAND."
 call ESMF_FieldGet(alvwf_nl_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = -1.e20

 print*,"- CALL FieldCreate FOR TARGET ALNSF AT NON-LAND."
 alnsf_nl_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="alnsf_nl_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET ALNSF AT NON-LAND."
 call ESMF_FieldGet(alnsf_nl_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = -1.e20

 print*,"- CALL FieldCreate FOR TARGET ALNWF AT NON-LAND."
 alnwf_nl_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="alnwf_nl_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET ALNWF AT NON-LAND."
 call ESMF_FieldGet(alnwf_nl_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = -1.e20
 endif ! fract_grid

 print*,"- CALL FieldCreate FOR TARGET GRID CANOPY MOISTURE CONTENT."
 canopy_mc_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="canopy_mc_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid canopy moisture."
 call ESMF_FieldGet(canopy_mc_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val
 
 print*,"- CALL FieldCreate FOR TARGET GRID LEAF AREA INDEX."
 lai_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="lai_target_grid",&
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET leaf area index."
 call ESMF_FieldGet(lai_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID Z0."
 z0_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="z0_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid z0."
 call ESMF_FieldGet(z0_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if (fract_grid)then
 print*,"- CALL FieldCreate FOR TARGET GRID Z0_ICE."
 z0_ice_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="z0_ice_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid z0_ice."
 call ESMF_FieldGet(z0_ice_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID Z0_WATER."
 z0_water_target_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="z0_water_target_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid z0_water."
 call ESMF_FieldGet(z0_water_target_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val
 endif

 print*,"- CALL FieldCreate FOR INTERPOLATED TARGET GRID TERRAIN."
 terrain_from_input_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     name="terrain_from_input_grid", &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid interpolated terrain."
 call ESMF_FieldGet(terrain_from_input_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 print*,"- CALL FieldCreate FOR INTERPOLATED TARGET GRID SOIL TYPE."
 soil_type_from_input_grid = ESMF_FieldCreate(target_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     name="soil_type_from_input_grid", rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid soil type"
 call ESMF_FieldGet(soil_type_from_input_grid, &
                    farrayPtr=target_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr = init_val

 if(fract_grid)then
 print*,"- CALL FieldCreate FOR TARGET GRID sea ice column TEMPERATURE."
 ice_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="ice_temp_target_grid", &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid ice temp"
 call ESMF_FieldGet(ice_temp_target_grid, &
                    farrayPtr=target_ptr_3d, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr_3d = init_val
 endif

 print*,"- CALL FieldCreate FOR TARGET GRID SOIL TEMPERATURE."
 soil_temp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="soil_temp_target_grid", &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid soil temp"
 call ESMF_FieldGet(soil_temp_target_grid, &
                    farrayPtr=target_ptr_3d, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr_3d = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID TOTAL SOIL MOISTURE."
 soilm_tot_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="soilm_tot_target_grid", &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid soil moist"
 call ESMF_FieldGet(soilm_tot_target_grid, &
                    farrayPtr=target_ptr_3d, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr_3d = init_val

 print*,"- CALL FieldCreate FOR TARGET GRID LIQUID SOIL MOISTURE."
 soilm_liq_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="soilm_liq_target_grid", &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- INITIALIZE TARGET grid soil liq"
 call ESMF_FieldGet(soilm_liq_target_grid, &
                    farrayPtr=target_ptr_3d, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 target_ptr_3d = init_val

 end subroutine create_surface_esmf_fields

!> Create ESMF fields for the target grid nst variables
!!
!! @author George Gayno
 subroutine create_nst_esmf_fields

 use model_grid, only               : target_grid

 implicit none

 integer                           :: rc

 print*,"- CALL FieldCreate FOR TARGET GRID C_D."
 c_d_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='c_d', &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID C_0."
 c_0_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='c_0', &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID D_CONV."
 d_conv_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='d_conv',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID DT_COOL."
 dt_cool_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='dt_cool',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID IFD."
 ifd_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='ifd',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID QRAIN."
 qrain_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='qrain',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID TREF."
 tref_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='tref',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID W_D."
 w_d_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='w_d',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID W_0."
 w_0_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='w_0',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XS."
 xs_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xs',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XT."
 xt_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xt',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XU."
 xu_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xu',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XV."
 xv_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xv',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XZ."
 xz_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xz',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XTTS."
 xtts_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xtts',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID XZTS."
 xzts_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='xzts',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID Z_C."
 z_c_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='z_c',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID ZM."
 zm_target_grid = ESMF_FieldCreate(target_grid, &
                                    typekind=ESMF_TYPEKIND_R8, &
                                     name='zm',&
                                    staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 end subroutine create_nst_esmf_fields

!> Update landmask for sea ice.
!!
!! @author George Gayno
 subroutine update_landmask

 use model_grid, only : landmask_target_grid, land_frac_target_grid

 use program_setup, only : fract_grid

 implicit none

 integer                        :: i, j, rc, clb(2), cub(2)
 integer(esmf_kind_i8), pointer :: mask_ptr(:,:)

 real(esmf_kind_r8), pointer    :: ice_ptr(:,:)
 real(esmf_kind_r8), pointer    :: land_frac_ptr(:,:)

 print*,"- UPDATE TARGET LANDMASK WITH ICE RECORD."

 print*,"- GET TARGET grid sea ice fraction."
 call ESMF_FieldGet(seaice_fract_target_grid, &
                    farrayPtr=ice_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- GET TARGET landmask."
 call ESMF_FieldGet(landmask_target_grid, &
                    farrayPtr=mask_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if (.not.fract_grid) then

   where(ice_ptr > 0.0) mask_ptr = 2

 else

   print*,"- GET TARGET land fraction."
   call ESMF_FieldGet(land_frac_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=land_frac_ptr, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

   do j = clb(2), cub(2)
   do i = clb(1), cub(1)

     if(ice_ptr(i,j) > 0.0) then
       mask_ptr(i,j) = 2
     else
       mask_ptr(i,j) = int(land_frac_ptr(i,j))
     endif
  
   enddo
   enddo

 endif

 end subroutine update_landmask

!> Convert 1d index to 2d indices.
!!
!! @param[in] ij  the 1d index
!! @param[in] itile  i-dimension of the tile
!! @param[in] jtile  j-dimension of the tile
!! @param[out] i  the "i" index
!! @param[out] j  the "j" index
!! @author George Gayno NOAA/EMC
 subroutine ij_to_i_j(ij, itile, jtile, i, j)

 implicit none

 integer(esmf_kind_i4), intent(in)  :: ij
 integer              , intent(in)  :: itile, jtile

 integer              , intent(out) :: i, j

 integer                            :: tile_num
 integer                            :: pt_loc_this_tile

 tile_num = ((ij-1) / (itile*jtile)) ! tile number minus 1
 pt_loc_this_tile = ij - (tile_num * itile * jtile)
                                     ! "ij" location of point within tile.

 j = (pt_loc_this_tile - 1) / itile + 1
 i = mod(pt_loc_this_tile, itile)

 if (i==0) i = itile

 return

 end subroutine ij_to_i_j

!> Regrid multiple ESMF fields from input to target grid
!!
!! @param[in] bundle_pre  ESMF fieldBundle on input grid
!! @param[in] bundle_post  ESMF fieldBundle on target grid
!! @param[in] num_field  Number of fields in target field pointer
!! @param[inout] route  Route handle to saved ESMF regridding instructions
!! @param[in]  dozero  Logical length num_field for whether field should be zeroed out before regridding
!! @param[inout]  unmapped_ptr (optional) Pointer to unmapped points from FieldRegrid
!! @param[in]  resetifd (optional) Logical for whether to reset ifd (only for water where nst data is used)
!! @author Larissa Reames, OU CIMMS/NOAA/NSSL
 subroutine regrid_many(bundle_pre,bundle_post, num_field,route,dozero, &
                        unmapped_ptr,resetifd)
 
 use esmf  
 use program_setup, only                : convert_nst
 use model_grid, only                   : i_target, j_target

 implicit none
 
 integer, intent(in)                    :: num_field
 type(esmf_routehandle), intent(inout)  :: route
 type(esmf_fieldbundle), intent(in)     :: bundle_pre, bundle_post
 logical, intent(in)                    :: dozero(num_field)
 logical, intent(in), optional       :: resetifd
 integer(esmf_kind_i4), intent(inout), optional  :: unmapped_ptr(:)
 
 type(esmf_field)                       :: field_pre,field_post
 real(esmf_kind_r8), pointer            :: tmp_ptr(:,:)
 type(realptr_2d),allocatable           :: ptr_2d(:)
 type(realptr_3d),allocatable           :: ptr_3d(:)
 logical                                :: is2d(num_field)
 character(len=50)                      :: fname
 integer :: i, j, k, ij, ind_2d, ind_3d, rc, ndims,n2d, n3d,localpet, l(1), u(1)
 type(esmf_vm) :: vm

 ind_2d = 0
 ind_3d = 0

 if(present(unmapped_ptr)) then
   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)
 endif
 
 do i = 1, num_field
   call ESMF_FieldBundleGet(bundle_pre,i,field_pre,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldBundleGet", rc)

   call ESMF_FieldBundleGet(bundle_post,i,field_post,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldBundleGet", rc)

   call ESMF_FieldGet(field_post,dimCount=ndims,name=fname,rc=rc)   
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
     call error_handler("IN FieldGet", rc)
   
   call ESMF_VMGetGlobal(vm, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN VMGetGlobal", rc)
   call ESMF_VMGet(vm, localPet=localpet, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN VMGet", rc)
   if(localpet==0) print*, "in regrid_many fname = ", fname, ndims
   if (ndims == 2) is2d(i) = .True.
   if (ndims == 3) is2d(i) = .False.
   
   if (dozero(i)) then
     call ESMF_FieldRegrid(field_pre, &
                           field_post, &
                           routehandle=route, &
                           termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)
   else
     call ESMF_FieldRegrid(field_pre, &
                           field_post, &
                           routehandle=route, &
                           zeroregion=ESMF_REGION_SELECT, &
                           termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)
   endif
 enddo
 
 if (present(resetifd)) then
   if( resetifd .and. convert_nst) then 
     call ESMF_FieldGet(ifd_target_grid,farrayPtr=tmp_ptr,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
     tmp_ptr = float(nint(tmp_ptr))
   endif
 endif
 
 n2d = count(is2d(:))
 n3d = count(.not.is2d(:))
 if(localpet==0) print*, is2d(:) 
 if (present(unmapped_ptr)) then
   allocate(ptr_2d(n2d))
   if (n3d .ne. 0) allocate(ptr_3d(n3d))
   do i=1, num_field
     if (is2d(i)) then 
       ind_2d = ind_2d + 1
       call ESMF_FieldBundleGet(bundle_post,i,field_post,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldBundleGet", rc)
       call ESMF_FieldGet(field_post, farrayPtr=ptr_2d(ind_2d)%p, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc) 
       call ESMF_FieldGet(field_post,name=fname,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGet", rc)
       if (localpet==0) print*, "in doreplace loop, 2d field = ", trim(fname)
     else
       ind_3d = ind_3d + 1
       call ESMF_FieldBundleGet(bundle_post,i,field_post,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldBundleGet", rc)
       call ESMF_FieldGet(field_post,name=fname,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGet", rc)
       if (localpet==0) print*, "in doreplace loop, 3d field = ", trim(fname)
       call ESMF_FieldGet(field_post, farrayPtr=ptr_3d(ind_3d)%p, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", rc) 
     endif
   end do
 
   do ij = l(1), u(1)
     call ij_to_i_j(unmapped_ptr(ij), i_target, j_target, i, j)
     do k = 1,n2d
       ptr_2d(k)%p(i,j) = -9999.9
     enddo 
     do k = 1,n3d
       ptr_3d(k)%p(i,j,:) = -9999.9
     enddo
   enddo
   deallocate(ptr_2d)
   if(n3d .ne. 0) deallocate(ptr_3d)
 endif
 end subroutine regrid_many

!> Execute the search function for multple fields
!!
!! @param[in] num_field  Number of fields to process.
!! @param[inout] bundle_target  ESMF FieldBundle holding target fields to search
!! @param[in] tile  Current cubed sphere tile.
!! @param[inout]  search_nums  Array length num_field holding search field numbers corresponding to each field provided for searching.
!! @param[in]  localpet  ESMF local persistent execution thread.
!! @param[in]  latitude  (optional) A real array size i_target,j_target of latitude on the target grid 
!! @param[in]  terrain_land  (optional) A real array size i_target,j_target of terrain height (m) on the target grid 
!! @param[in]  soilt_climo  (optional) A real array size i_target,j_target of climatological soil type on the target grid 
!! @param[inout] mask  (optional) An integer array of size i_target,j_target that holds masked (0) and unmasked (1)
!!                     values indicating where to execute search (only at
!unmasked points).
!! @author Larissa Reames, OU CIMMS/NOAA/NSSL
 subroutine search_many(num_field,bundle_target,tile,search_nums,localpet,latitude, &
                        terrain_land,soilt_climo, mask)

 use model_grid, only                  : i_target,j_target, lsoil_target
 use program_setup, only               : external_model, input_type
 use search_util

 implicit none

 integer, intent(in)                         :: num_field
 type(esmf_fieldbundle), intent(inout)       :: bundle_target

 real(esmf_kind_r8), intent(inout), optional :: latitude(i_target,j_target)
 real(esmf_kind_r8), intent(inout), optional :: terrain_land(i_target,j_target)
 real(esmf_kind_r8), intent(inout), optional :: soilt_climo(i_target,j_target)
 integer(esmf_kind_i8), intent(inout), optional  :: mask(i_target,j_target)
 
 real(esmf_kind_r8), allocatable :: field_data_2d(:,:)   
 real(esmf_kind_r8), allocatable :: field_data_3d(:,:,:)  
 integer, intent(in)             :: tile,localpet
 integer, intent(inout)          :: search_nums(num_field)
 
 type(esmf_field)                :: temp_field
 character(len=50)               :: fname
 integer, parameter              :: SOTYP_LAND_FIELD_NUM = 224
 integer, parameter              :: SST_FIELD_NUM = 11
 integer, parameter              :: TERRAIN_FIELD_NUM= 7
 integer :: j,k, rc, ndims

 
 do k = 1,num_field
   call ESMF_FieldBundleGet(bundle_target,k,temp_field, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", rc)
   call ESMF_FieldGet(temp_field, name=fname, dimcount=ndims,rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldGet", rc)
   if (localpet==0) then
       allocate(field_data_2d(i_target,j_target))
   else
       allocate(field_data_2d(0,0))
   endif
   if (ndims .eq. 2) then
       call ESMF_FieldGather(temp_field,field_data_2d,rootPet=0,tile=tile, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)
     if (localpet == 0) then 
       if (present(latitude) .and. search_nums(k).eq.SST_FIELD_NUM) then
         ! Sea surface temperatures; pass latitude field to search
         call search(field_data_2d, mask, i_target, j_target, tile,search_nums(k),latitude=latitude)
       elseif (present(terrain_land) .and. search_nums(k) .eq. TERRAIN_FIELD_NUM) then
         ! Terrain height; pass optional climo terrain array to search
         call search(field_data_2d, mask, i_target, j_target, tile,search_nums(k),terrain_land=terrain_land)
       elseif (search_nums(k) .eq. SOTYP_LAND_FIELD_NUM) then
         ! Soil type over land    
         if (fname .eq. "soil_type_target_grid") then
           ! Soil type over land when interpolating input data to target grid
           ! *with* the intention of retaining interpolated data in output
           call search(field_data_2d, mask, i_target, j_target, tile,search_nums(k),soilt_climo=soilt_climo)
         elseif (present(soilt_climo)) then
           if (maxval(field_data_2d) > 0 .and. (trim(external_model) .ne. "GFS" .or. trim(input_type) .ne. "grib2")) then
             ! Soil type over land when interpolating input data to target grid
             ! *without* the intention of retaining data in output file
             call search(field_data_2d, mask, i_target, j_target, tile, search_nums(k))
           else 
             ! If no soil type field exists in input data (e.g., GFS grib2) then don't search
             ! but simply set data to the climo field. This may result in
             ! somewhat inaccurate soil moistures as no scaling will occur 
             field_data_2d = soilt_climo
           endif !check field value   
         endif !sotype from target grid
       else
         ! Any field that doesn't require any of the special treatments or
         ! passing of additional variables as in those above
         call search(field_data_2d, mask, i_target, j_target, tile,search_nums(k))
       endif !if present  
     endif !localpet
     call ESMF_FieldScatter(temp_field, field_data_2d, rootPet=0, tile=tile,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)
   else
     if (localpet==0) then
         allocate(field_data_3d(i_target,j_target,lsoil_target))
     else
         allocate(field_data_3d(0,0,0))
     endif
 
     ! Process 3d fields soil temperature, moisture, and liquid
     call ESMF_FieldGather(temp_field,field_data_3d,rootPet=0,tile=tile,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldGather", rc)
        
     if (localpet==0) then 
       do j = 1, lsoil_target
         field_data_2d = field_data_3d(:,:,j)
         call search(field_data_2d, mask, i_target, j_target, tile, 21)
         field_data_3d(:,:,j) = field_data_2d
       enddo
     endif
     call ESMF_FieldScatter(temp_field, field_data_3d, rootPet=0, tile=tile,rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)
     deallocate(field_data_3d)
   endif !ndims
   deallocate(field_data_2d)
 end do !fields

 end subroutine search_many

!> Free up memory once the target grid surface fields are
!! no longer needed.
!!
!! @author George Gayno NOAA/EMC
 subroutine cleanup_all_target_sfc_data

 use surface_target_data, only : cleanup_target_sfc_data

 implicit none

 integer                     :: rc

 print*,"- DESTROY LOCAL TARGET GRID SURFACE FIELDS."

 call ESMF_FieldDestroy(terrain_from_input_grid, rc=rc)
 call ESMF_FieldDestroy(terrain_from_input_grid_land, rc=rc)
 call ESMF_FieldDestroy(soil_type_from_input_grid, rc=rc)

 call cleanup_target_sfc_data

 end subroutine cleanup_all_target_sfc_data

 end module surface
