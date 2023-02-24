!> @file
!! @brief Process atmospheric fields.
!! @author George Gayno NCEP/EMC

!> Process atmospheric fields. Horizontally interpolate from input to
!! target FV3 grid using ESMF regridding. Adjust surface pressure
!! according to terrain differences between input and target
!! grid. Vertically interpolate to target FV3 grid vertical
!! levels. Processing based on the spectral GFS version of CHGRES.
!! 
!! For variables "b4adj" indicates fields on the target grid before
!! vertical adjustment. "target" indicates data on target grid.
!! "input" indicates data on input grid. "_s" indicates fields on the
!! 'south' edge of the grid box.  "_w" indicate fields on the 'west'
!! edge of the grid box.  Otherwise, fields are at the center of the
!! grid box.
!!
!! @author George Gayno NCEP/EMC
 module atmosphere

 use esmf

 use atmosphere_target_data, only    : lev_target, levp1_target, nvcoord_target, &
                                       vcoord_target, delp_target_grid, &
                                       dzdt_target_grid, ps_target_grid, &
                                       temp_target_grid, tracers_target_grid, &
                                       u_s_target_grid, v_s_target_grid, &
                                       u_w_target_grid, v_w_target_grid, &
                                       zh_target_grid, qnwfa_climo_target_grid, &
                                       qnifa_climo_target_grid

 use atm_input_data, only            : lev_input, &
                                       levp1_input, &
                                       tracers_input_grid, &
                                       dzdt_input_grid, &
                                       ps_input_grid, &
                                       xwind_input_grid,   &
                                       ywind_input_grid,   &
                                       zwind_input_grid,   &
                                       temp_input_grid,   &
                                       pres_input_grid,   &
                                       terrain_input_grid, &
                                       read_input_atm_data, &
                                       cleanup_input_atm_data

 use model_grid, only                : target_grid,  &
                                       latitude_s_target_grid,  &
                                       longitude_s_target_grid, &
                                       latitude_w_target_grid,  &
                                       longitude_w_target_grid, &
                                       terrain_target_grid

 use program_setup, only             : vcoord_file_target_grid, &
                                       wam_cold_start, wam_parm_file, & 
                                       cycle_year, cycle_mon,     &
                                       cycle_day, cycle_hour,     &
                                       regional, &
                                       tracers, num_tracers,      &
                                       num_tracers_input,         & 
                                       atm_weight_file, &
                                       use_thomp_mp_climo

 use thompson_mp_climo_data, only    : read_thomp_mp_climo_data,  &
                                       cleanup_thomp_mp_climo_input_data, &
                                       qnifa_climo_input_grid, &
                                       qnwfa_climo_input_grid, &
                                       thomp_pres_climo_input_grid, &
                                       lev_thomp_mp_climo

 use write_data, only                : write_fv3_atm_header_netcdf, &
                                       write_fv3_atm_bndy_data_netcdf, &
                                       write_fv3_atm_data_netcdf

 use utilities, only                 : error_handler

 implicit none

 private

 type(esmf_field)                       :: dzdt_b4adj_target_grid !< vertical vel before vert adj
 type(esmf_field), allocatable          :: tracers_b4adj_target_grid(:) !< tracers before vert adj
 type(esmf_field)                       :: ps_b4adj_target_grid !< sfc pres before terrain adj
 type(esmf_field)                       :: pres_target_grid !< 3-d pressure
 type(esmf_field)                       :: pres_b4adj_target_grid !< 3-d pres before terrain adj
 type(esmf_field)                       :: temp_b4adj_target_grid !< temp before vert adj
 type(esmf_field)                       :: terrain_interp_to_target_grid !< Input grid terrain interpolated to target grid.   
 type(esmf_field)                       :: xwind_target_grid !< x-component wind, grid box center
 type(esmf_field)                       :: ywind_target_grid !< y-component wind, grid box center
 type(esmf_field)                       :: zwind_target_grid !< z-component wind, grid box center
 type(esmf_field)                       :: xwind_b4adj_target_grid !< x-component wind, before vert adj
 type(esmf_field)                       :: ywind_b4adj_target_grid !< y-component wind, before vert adj
 type(esmf_field)                       :: zwind_b4adj_target_grid !< z-component wind, before vert adj
 type(esmf_field)                       :: xwind_s_target_grid !< x-component wind, 'south' edge
 type(esmf_field)                       :: ywind_s_target_grid !< y-component wind, 'south' edge
 type(esmf_field)                       :: zwind_s_target_grid !< z-component wind, 'south' edge
 type(esmf_field)                       :: xwind_w_target_grid !< x-component wind, 'west' edge
 type(esmf_field)                       :: ywind_w_target_grid !< y-component wind, 'west' edge
 type(esmf_field)                       :: zwind_w_target_grid !< z-component wind, 'west' edge

! Fields associated with thompson microphysics climatological tracers.

 type(esmf_field)                       :: qnifa_climo_b4adj_target_grid !< number concentration of ice
                                           !! friendly aerosols before vert adj
 type(esmf_field)                       :: qnwfa_climo_b4adj_target_grid !< number concentration of water
                                           !! friendly aerosols before vert adj
 type(esmf_field)                       :: thomp_pres_climo_b4adj_target_grid !< pressure of each level on
                                           !! target grid

 public :: atmosphere_driver
 public :: read_vcoord_info

 contains

!> Driver routine to process for atmospheric fields.
!!
!! @param[in] localpet ESMF local persistent execution thread 
!! @author George Gayno
 subroutine atmosphere_driver(localpet)

 use mpi

 implicit none

 integer, intent(in)                :: localpet

 integer                            :: isrctermprocessing
 integer                            :: rc, n

 type(esmf_regridmethod_flag)       :: method
 type(esmf_routehandle)             :: regrid_bl

 real(esmf_kind_r8), parameter      :: p0=101325.0
 real(esmf_kind_r8), parameter      :: rd = 287.058
 real(esmf_kind_r8), parameter      :: grav = 9.81
 real(esmf_kind_r8), parameter      :: lapse = -6.5e-03

 real(esmf_kind_r8), parameter      :: exponent = rd*lapse/grav
 real(esmf_kind_r8), parameter      :: one_over_exponent = 1.0 / exponent

 real(esmf_kind_r8), pointer        :: psptr(:,:), tempptr(:,:,:)

!-----------------------------------------------------------------------------------
! Read atmospheric fields on the input grid.
!-----------------------------------------------------------------------------------

 call read_input_atm_data(localpet)

!-----------------------------------------------------------------------------------
! Read vertical coordinate info for target grid.
!-----------------------------------------------------------------------------------

 call read_vcoord_info

!-----------------------------------------------------------------------------------
! Create target grid field objects to hold data before vertical adjustment.
!-----------------------------------------------------------------------------------

 call create_atm_b4adj_esmf_fields

!-----------------------------------------------------------------------------------
! Horizontally interpolate.  If specified, use weights from file.
!-----------------------------------------------------------------------------------

 isrctermprocessing = 1

 if (trim(atm_weight_file) /= "NULL") then

   print*,"- CALL FieldSMMStore FOR ATMOSPHERIC FIELDS."

   call ESMF_FieldSMMStore(temp_input_grid, &
                           temp_b4adj_target_grid, &
                           atm_weight_file, &
                           routehandle=regrid_bl, &
                           srctermprocessing=isrctermprocessing, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldSMMStore", rc)

 else

   print*,"- CALL FieldRegridStore FOR ATMOSPHERIC FIELDS."

   method=ESMF_REGRIDMETHOD_BILINEAR

   call ESMF_FieldRegridStore(temp_input_grid, &
                              temp_b4adj_target_grid, &
                              polemethod=ESMF_POLEMETHOD_ALLAVG, &
                              srctermprocessing=isrctermprocessing, &
                              routehandle=regrid_bl, &
                              regridmethod=method, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridStore", rc)

 endif

 print*,"- CALL Field_Regrid FOR TEMPERATURE."
 call ESMF_FieldRegrid(temp_input_grid, &
                       temp_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR PRESSURE."
 call ESMF_FieldRegrid(pres_input_grid, &
                       pres_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 do n = 1, num_tracers_input
   print*,"- CALL Field_Regrid FOR TRACER ", trim(tracers(n))
   call ESMF_FieldRegrid(tracers_input_grid(n), &
                         tracers_b4adj_target_grid(n), &
                         routehandle=regrid_bl, &
                         termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)
      
 enddo

 print*,"- CALL Field_Regrid FOR VERTICAL VELOCITY."
 call ESMF_FieldRegrid(dzdt_input_grid, &
                       dzdt_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)
    
 nullify(tempptr)
 print*,"- CALL FieldGet FOR INPUT GRID VERTICAL VEL."
 call ESMF_FieldGet(dzdt_input_grid, &
                    farrayPtr=tempptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 print*, "MIN MAX W INPUT = ", minval(tempptr), maxval(tempptr)

 nullify(tempptr)
 print*,"- CALL FieldGet FOR VERTICAL VEL B4ADJ."
 call ESMF_FieldGet(dzdt_b4adj_target_grid, &
                    farrayPtr=tempptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 print*, "MIN MAX W B4ADJ = ", minval(tempptr), maxval(tempptr)
 
 nullify(psptr)
 print*,"- CALL FieldGet FOR INPUT SURFACE PRESSURE."
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!------------------------------------------------------------------------------------
! Assume standard lapse rate when interpolating pressure (per Phil Pegion).
!------------------------------------------------------------------------------------

 psptr = (psptr/p0)**exponent

 print*,"- CALL Field_Regrid FOR SURFACE PRESSURE."
 call ESMF_FieldRegrid(ps_input_grid, &
                       ps_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 nullify(psptr)
 print*,"- CALL FieldGet FOR INPUT SURFACE PRESSURE B4ADJ."
 call ESMF_FieldGet(ps_b4adj_target_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 psptr = p0 * psptr**one_over_exponent

 print*,"- CALL Field_Regrid FOR TERRAIN."
 call ESMF_FieldRegrid(terrain_input_grid, &
                       terrain_interp_to_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, &
                       rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR x WIND."
 call ESMF_FieldRegrid(xwind_input_grid, &
                       xwind_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR y WIND."
 call ESMF_FieldRegrid(ywind_input_grid, &
                       ywind_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR z WIND."
 call ESMF_FieldRegrid(zwind_input_grid, &
                       zwind_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", rc)



 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)

!-----------------------------------------------------------------------------------
! Deallocate input fields.
!-----------------------------------------------------------------------------------

 call cleanup_input_atm_data

!-----------------------------------------------------------------------------------
! Create target grid field objects to hold data after vertical interpolation.
!-----------------------------------------------------------------------------------

 call create_atm_esmf_fields

!-----------------------------------------------------------------------------------
! Adjust surface pressure for terrain differences.
!-----------------------------------------------------------------------------------

 call newps(localpet)

!-----------------------------------------------------------------------------------
! Compute 3-d pressure based on adjusted surface pressure.
!-----------------------------------------------------------------------------------

 call newpr1(localpet)

!-----------------------------------------------------------------------------------
! Vertically interpolate.
!-----------------------------------------------------------------------------------

 call vintg

 if( wam_cold_start ) then 
   call vintg_wam (cycle_year,cycle_mon,cycle_day,cycle_hour,wam_parm_file)
 endif

!-----------------------------------------------------------------------------------
! Compute height.
!-----------------------------------------------------------------------------------

 call compute_zh

!-----------------------------------------------------------------------------------
! Free up memory.
!-----------------------------------------------------------------------------------

 call cleanup_target_atm_b4adj_data

!-----------------------------------------------------------------------------------
! Interpolate winds to 'd' grid.
!-----------------------------------------------------------------------------------

 isrctermprocessing = 1
 method=ESMF_REGRIDMETHOD_BILINEAR

 print*,"- CALL FieldRegridStore FOR X-WIND WEST EDGE."
 call ESMF_FieldRegridStore(xwind_target_grid, &
                            xwind_w_target_grid, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 print*,"- CALL Field_Regrid FOR X-WIND WEST EDGE."
 call ESMF_FieldRegrid(xwind_target_grid, &
                       xwind_w_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR Y-WIND WEST EDGE."
 call ESMF_FieldRegrid(ywind_target_grid, &
                       ywind_w_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR Z-WIND WEST EDGE."
 call ESMF_FieldRegrid(zwind_target_grid, &
                       zwind_w_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

 isrctermprocessing = 1
 method=ESMF_REGRIDMETHOD_BILINEAR

 print*,"- CALL FieldRegridStore FOR X-WIND SOUTH EDGE."
 call ESMF_FieldRegridStore(xwind_target_grid, &
                            xwind_s_target_grid, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", rc)

 print*,"- CALL Field_Regrid FOR X-WIND SOUTH EDGE."
 call ESMF_FieldRegrid(xwind_target_grid, &
                       xwind_s_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR Y-WIND SOUTH EDGE."
 call ESMF_FieldRegrid(ywind_target_grid, &
                       ywind_s_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR Z-WIND SOUTH EDGE."
 call ESMF_FieldRegrid(zwind_target_grid, &
                       zwind_s_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", rc)

!-----------------------------------------------------------------------------------
! Convert from 3-d to 2-d cartesian winds.
!-----------------------------------------------------------------------------------

 call convert_winds_to_uv
 
!-----------------------------------------------------------------------------------
! If selected, process thompson microphysics climatological fields.
!-----------------------------------------------------------------------------------

 if (use_thomp_mp_climo) then
   call read_thomp_mp_climo_data
   call horiz_interp_thomp_mp_climo
   call vintg_thomp_mp_climo
 endif 

!-----------------------------------------------------------------------------------
! Write target data to file.
!-----------------------------------------------------------------------------------

 call write_fv3_atm_header_netcdf(localpet)
 if (regional <= 1) call write_fv3_atm_data_netcdf(localpet)
 if (regional >= 1) call write_fv3_atm_bndy_data_netcdf(localpet)

!-----------------------------------------------------------------------------------
! Free up memory.
!-----------------------------------------------------------------------------------

 call cleanup_all_target_atm_data

 end subroutine atmosphere_driver

!> Create target grid field objects to hold data before vertical
!! interpolation. These will be defined with the same number of
!! vertical levels as the input grid.
!!
!! @author George Gayno
 subroutine create_atm_b4adj_esmf_fields

 implicit none

 integer                          :: rc, n

 allocate(tracers_b4adj_target_grid(num_tracers_input))

 do n = 1, num_tracers_input
   print*,"- CALL FieldCreate FOR TARGET GRID TRACER BEFORE ADJUSTMENT ", trim(tracers(n))
   tracers_b4adj_target_grid(n) = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldCreate", rc)
 enddo

 print*,"- CALL FieldCreate FOR TARGET GRID TEMPERATURE BEFORE ADJUSTMENT."
 temp_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID PRESSURE BEFORE ADJUSTMENT."
 pres_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID VERTICAL VELOCITY BEFORE ADJUSTMENT."
 dzdt_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID xwind."
 xwind_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID ywind."
 ywind_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID zwind."
 zwind_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET TERRAIN."
 terrain_interp_to_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET SURFACE PRESSURE BEFORE ADJUSTMENT."
 ps_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 end subroutine create_atm_b4adj_esmf_fields

!> Create target grid field objects.
!!
!! @author George Gayno
 subroutine create_atm_esmf_fields

 implicit none

 integer                          :: rc, n

 allocate(tracers_target_grid(num_tracers))

 do n = 1, num_tracers
    print*,"- CALL FieldCreate FOR TARGET GRID TRACERS ", trim(tracers(n))    
    tracers_target_grid(n) = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldCreate", rc)
 enddo

 print*,"- CALL FieldCreate FOR TARGET GRID TEMPERATURE."
 temp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID PRESSURE."
 pres_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID VERTICAL VELOCITY."
 dzdt_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID DELP."
 delp_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID xwind."
 xwind_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID ywind."
 ywind_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID zwind."
 zwind_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET HEIGHT."
 zh_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/levp1_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET U_S."
 u_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET V_S."
 v_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET xwind_S."
 xwind_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET ywind_S."
 ywind_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET zwind_S."
 zwind_s_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET U_W."
 u_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET V_W."
 v_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET xwind_W."
 xwind_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET ywind_W."
 ywind_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET zwind_W."
 zwind_w_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET SURFACE PRESSURE."
 ps_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 end subroutine create_atm_esmf_fields

!> Convert 3-d component winds to u and v.
!!
!! @author George Gayno
 subroutine convert_winds_to_uv
 
 implicit none

 integer                         :: clb(3), cub(3)
 integer                         :: i, j, k, rc

 real(esmf_kind_r8), pointer     :: latptr(:,:)
 real(esmf_kind_r8), pointer     :: lonptr(:,:)
 real(esmf_kind_r8), pointer     :: uptr(:,:,:)
 real(esmf_kind_r8), pointer     :: vptr(:,:,:)
 real(esmf_kind_r8), pointer     :: xwindptr(:,:,:)
 real(esmf_kind_r8), pointer     :: ywindptr(:,:,:)
 real(esmf_kind_r8), pointer     :: zwindptr(:,:,:)
 real(esmf_kind_r8)              :: latrad, lonrad

!-----------------------------------------------------------------------------------
! Convert from 3-d cartesian to 2-cartesian winds
!-----------------------------------------------------------------------------------

 print*,'- CONVERT WINDS.'

 print*,"- CALL FieldGet FOR xwind_S."
 call ESMF_FieldGet(xwind_s_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=xwindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ywind_S."
 call ESMF_FieldGet(ywind_s_target_grid, &
                    farrayPtr=ywindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR zwind_S."
 call ESMF_FieldGet(zwind_s_target_grid, &
                    farrayPtr=zwindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR U_S."
 call ESMF_FieldGet(u_s_target_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR V_S."
 call ESMF_FieldGet(v_s_target_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LATITUDE_S."
 call ESMF_FieldGet(latitude_s_target_grid, &
                    farrayPtr=latptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LONGITUDE_S."
 call ESMF_FieldGet(longitude_s_target_grid, &
                    farrayPtr=lonptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     latrad = latptr(i,j) * acos(-1.) / 180.0
     lonrad = lonptr(i,j) * acos(-1.) / 180.0
     do k = clb(3), cub(3)
       uptr(i,j,k) = xwindptr(i,j,k) * cos(lonrad) + ywindptr(i,j,k) * sin(lonrad)
       vptr(i,j,k) = -xwindptr(i,j,k) * sin(latrad) * sin(lonrad) + &
                      ywindptr(i,j,k) * sin(latrad) * cos(lonrad) + &
                      zwindptr(i,j,k) * cos(latrad)
     enddo
   enddo
 enddo

 print*,"- CALL FieldGet FOR xwind_w."
 call ESMF_FieldGet(xwind_w_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=xwindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ywind_w."
 call ESMF_FieldGet(ywind_w_target_grid, &
                    farrayPtr=ywindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR zwind_w."
 call ESMF_FieldGet(zwind_w_target_grid, &
                    farrayPtr=zwindptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR U_W."
 call ESMF_FieldGet(u_w_target_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR V_W."
 call ESMF_FieldGet(v_w_target_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LATITUDE_W."
 call ESMF_FieldGet(latitude_w_target_grid, &
                    farrayPtr=latptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LONGITUDE_W."
 call ESMF_FieldGet(longitude_w_target_grid, &
                    farrayPtr=lonptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     latrad = latptr(i,j) * acos(-1.) / 180.0
     lonrad = lonptr(i,j) * acos(-1.) / 180.0
     do k = clb(3), cub(3)
       uptr(i,j,k) = xwindptr(i,j,k) * cos(lonrad) + ywindptr(i,j,k) * sin(lonrad)
       vptr(i,j,k) = -xwindptr(i,j,k) * sin(latrad) * sin(lonrad) + &
                      ywindptr(i,j,k) * sin(latrad) * cos(lonrad) + &
                      zwindptr(i,j,k) * cos(latrad)
     enddo
   enddo
 enddo

 end subroutine convert_winds_to_uv

!> Computes 3-D pressure given an adjusted surface pressure.
!!                                                                       
!! program history log:                                                  
!! 2005-04-11  Hann-Ming Henry Juang    hybrid sigma, sigma-p, and sigma-
!! - PRGMMR: Henry Juang    ORG: W/NMC23     DATE: 2005-04-11            
!! - PRGMMR: Fanglin Yang   ORG: W/NMC23     DATE: 2006-11-28            
!! - PRGMMR: S. Moorthi     ORG: NCEP/EMC    DATE: 2006-12-12            
!! - PRGMMR: S. Moorthi     ORG: NCEP/EMC    DATE: 2007-01-02            
!!                                                                       
!!   INPUT ARGUMENT LIST:                                                
!!     IM           INTEGER NUMBER OF POINTS TO COMPUTE                  
!!     KM           INTEGER NUMBER OF LEVELS                             
!!     IDVC         INTEGER VERTICAL COORDINATE ID                       
!!                  (1 FOR SIGMA AND 2 FOR HYBRID)                       
!!     IDSL         INTEGER TYPE OF SIGMA STRUCTURE                      
!!                  (1 FOR PHILLIPS OR 2 FOR MEAN)                       
!!     NVCOORD      INTEGER NUMBER OF VERTICAL COORDINATES               
!!     VCOORD       REAL (KM+1,NVCOORD) VERTICAL COORDINATE VALUES       
!!                  FOR IDVC=1, NVCOORD=1: SIGMA INTERFACE               
!!                  FOR IDVC=2, NVCOORD=2: HYBRID INTERFACE A AND B      
!!                  FOR IDVC=3, NVCOORD=3: JUANG GENERAL HYBRID INTERFACE
!!                     AK  REAL (KM+1) HYBRID INTERFACE A                
!!                     BK  REAL (KM+1) HYBRID INTERFACE B                
!!     PS           REAL (IX) SURFACE PRESSURE (PA)                      
!!   OUTPUT ARGUMENT LIST:                                               
!!     PM           REAL (IX,KM) MID-LAYER PRESSURE (PA)                 
!!     DP           REAL (IX,KM) LAYER DELTA PRESSURE (PA)
!!
!! @param[in] localpet ESMF local persistent execution thread  
!!
!! @author Hann Ming Henry Juang, Juang, Fanglin Yang, S. Moorthi
 subroutine newpr1(localpet)
 implicit none 

 integer, intent(in) :: localpet

 integer                         :: idsl, idvc, rc
 integer                         :: i, j, k, clb(3), cub(3)

 real(esmf_kind_r8), parameter   :: rd=287.05
 real(esmf_kind_r8), parameter   :: cp=1004.6
 real(esmf_kind_r8), parameter   :: rocp=rd/cp
 real(esmf_kind_r8), parameter   :: rocp1=rocp+1
 real(esmf_kind_r8), parameter   :: rocpr=1/rocp

 real(esmf_kind_r8), pointer     :: delp_ptr(:,:,:)
 real(esmf_kind_r8), pointer     :: pptr(:,:,:)    ! adjusted 3-d p.
 real(esmf_kind_r8), pointer     :: psptr(:,:)  ! adjusted surface p.
 real(esmf_kind_r8)              :: ak, bk
 real(esmf_kind_r8), allocatable :: pi(:,:,:)
                
 print*,"COMPUTE 3-D PRESSURE FROM ADJUSTED SURFACE PRESSURE."

 idvc = 2 ! hard wire for now.
 idsl = 2 ! hard wire for now.

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR DELP."
 call ESMF_FieldGet(delp_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=delp_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE AFTER ADJUSTMENT"
 call ESMF_FieldGet(ps_target_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(pi(clb(1):cub(1),clb(2):cub(2),1:levp1_target))
 
 if(idvc.eq.2) then
   do k=1,levp1_target
     ak = vcoord_target(k,1) 
     bk = vcoord_target(k,2) 
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pi(i,j,k) = ak + bk*psptr(i,j)
       enddo
     enddo
   enddo 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         delp_ptr(i,j,k) = pi(i,j,k) - pi(i,j,k+1)
       enddo
     enddo
   enddo
 else 
   call error_handler("PROGRAM ONLY WORKS WITH IDVC 2", 1)
 endif

 if(idsl.eq.2) then 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pptr(i,j,k) = (pi(i,j,k)+pi(i,j,k+1))/2.0
       enddo
     enddo
   enddo 
 else 
   do k=1,lev_target
     do i= clb(1), cub(1)
       do j= clb(2), cub(2)
         pptr(i,j,k) = ((pi(i,j,k)**rocp1-pi(i,j,k+1)**rocp1)/        &
                        (rocp1*(pi(i,j,k)-pi(i,j,k+1))))**rocpr        
       enddo
     enddo
   enddo 
 endif 

 deallocate(pi)

 if (localpet == 0) then
    print*,'new pres ',pptr(clb(1),clb(2),:)
    print*,'delp     ',delp_ptr(clb(1),clb(2),:)
 endif

 end subroutine newpr1 

!> Computes adjusted surface pressure given a new terrain height.
!!
!! Computes a new surface pressure given a new orography. The new
!! pressure is computed assuming a hydrostatic balance and a constant
!! temperature lapse rate. Below ground, the lapse rate is assumed to
!! be -6.5 k/km.
!!
!! program history log:
!! -  91-10-31  mark iredell
!! -  2018-apr  adapt for fv3. george gayno
!!
!! @param[in] localpet ESMF local persistent execution thread 
!! @author Mark Iredell, George Gayno @date 92-10-31
 subroutine newps(localpet)

 implicit none

 integer, intent(in) :: localpet
 integer                         :: i, j, k, ii
 integer                         :: clb(3), cub(3), ls, rc

 real(esmf_kind_r8), pointer     :: pptr(:,:,:)
 real(esmf_kind_r8), pointer     :: psptr(:,:)
 real(esmf_kind_r8), pointer     :: psnewptr(:,:)  ! adjusted surface p.
 real(esmf_kind_r8), pointer     :: tptr(:,:,:)
 real(esmf_kind_r8), pointer     :: qptr(:,:,:)
 real(esmf_kind_r8), pointer     :: zsptr(:,:)
 real(esmf_kind_r8), pointer     :: zsnewptr(:,:)
 real(esmf_kind_r8), allocatable :: zu(:,:)
 real(esmf_kind_r8), parameter   :: beta=-6.5E-3
 real(esmf_kind_r8), parameter   :: epsilon=1.E-9
 real(esmf_kind_r8), parameter   :: g=9.80665
 real(esmf_kind_r8), parameter   :: rd=287.05
 real(esmf_kind_r8), parameter   :: rv=461.50
 real(esmf_kind_r8), parameter   :: gor=g/rd
 real(esmf_kind_r8), parameter   :: fv=rv/rd-1.
 real(esmf_kind_r8)              :: ftv, fgam, apu, fz0
 real(esmf_kind_r8)              :: atvu, atv, fz1, fp0
 real(esmf_kind_r8)              :: apd, azd, agam, azu
 real(esmf_kind_r8)              :: atvd, fp1, gamma, pu
 real(esmf_kind_r8)              :: tvu, pd, tvd
 real(esmf_kind_r8)              :: at, aq, ap, az

 ftv(at,aq)=at*(1+fv*aq)
 fgam(apu,atvu,apd,atvd)=-gor*log(atvd/atvu)/log(apd/apu)
 fz0(ap,atv,azd,apd)=azd+atv/gor*log(apd/ap)
 fz1(ap,atv,azd,apd,agam)=azd-atv/agam*((apd/ap)**(-agam/gor)-1)
 fp0(az,azu,apu,atvu)=apu*exp(-gor/atvu*(az-azu))
 fp1(az,azu,apu,atvu,agam)=apu*(1+agam/atvu*(az-azu))**(-gor/agam)

 print*,"- ADJUST SURFACE PRESSURE FOR NEW TERRAIN."

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 if(localpet==0) then
   print*,'old pres ',pptr(clb(1),clb(2),:)
 endif

 print*,"- CALL FieldGet FOR TEMPERATURE"
 call ESMF_FieldGet(temp_b4adj_target_grid, &
                    farrayPtr=tptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
! Find specific humidity in the array of tracer fields.

 do ii = 1, num_tracers
   if (trim(tracers(ii)) == "sphum") exit
 enddo

 print*,"- CALL FieldGet FOR SPECIFIC HUMIDITY"
 call ESMF_FieldGet(tracers_b4adj_target_grid(ii), &
                    farrayPtr=qptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 print*,"- CALL FieldGet FOR SURFACE PRESSURE BEFORE ADJUSTMENT"
 call ESMF_FieldGet(ps_b4adj_target_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE AFTER ADJUSTMENT"
 call ESMF_FieldGet(ps_target_grid, &
                    farrayPtr=psnewptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR OLD TERRAIN"
 call ESMF_FieldGet(terrain_interp_to_target_grid, &
                    farrayPtr=zsptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR NEW TERRAIN"
 call ESMF_FieldGet(terrain_target_grid, &
                    farrayPtr=zsnewptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(zu(clb(1):cub(1),clb(2):cub(2)))

!-----------------------------------------------------------------------------------
! Note, this routine was adapted from the spectral GFS which labeled the lowest
! model layer as '1'.  
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! Compute surface pressure below the original ground.
!-----------------------------------------------------------------------------------

 ls=0
 k=1
 gamma=beta
 do i=clb(1), cub(1)
 do j=clb(2), cub(2)
   pu=pptr(i,j,k)
   tvu=ftv(tptr(i,j,k),qptr(i,j,k))
   zu(i,j)=fz1(pu,tvu,zsptr(i,j),psptr(i,j),gamma)
   if(zsnewptr(i,j).le.zu(i,j)) then
     pu=pptr(i,j,k)
     tvu=ftv(tptr(i,j,k),qptr(i,j,k))
     if(abs(gamma).gt.epsilon) then
       psnewptr(i,j)=fp1(zsnewptr(i,j),zu(i,j),pu,tvu,gamma)
     else
       psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
     endif
   else
     psnewptr(i,j)=0
     ls=ls+1
   endif
 enddo
 enddo

!-----------------------------------------------------------------------------------
! Compute surface pressure above the original ground.
!-----------------------------------------------------------------------------------

 do k=2,cub(3)
   if(ls.gt.0) then
     do i=clb(1),cub(1)
     do j=clb(2),cub(2)
       if(psnewptr(i,j).eq.0) then
         pu=pptr(i,j,k)
         tvu=ftv(tptr(i,j,k),qptr(i,j,k))
         pd=pptr(i,j,k-1)
         tvd=ftv(tptr(i,j,k-1),qptr(i,j,k-1))
         gamma=fgam(pu,tvu,pd,tvd)
         if(abs(gamma).gt.epsilon) then
           zu(i,j)=fz1(pu,tvu,zu(i,j),pd,gamma)
         else
           zu(i,j)=fz0(pu,tvu,zu(i,j),pd)
         endif
         if(zsnewptr(i,j).le.zu(i,j)) then
           if(abs(gamma).gt.epsilon) then
             psnewptr(i,j)=fp1(zsnewptr(i,j),zu(i,j),pu,tvu,gamma)
           else
             psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
           endif
           ls=ls-1
         endif
       endif
     enddo
     enddo
   endif
 enddo

!-----------------------------------------------------------------------------------
! Compute surface pressure over the top.
!-----------------------------------------------------------------------------------


 if(ls.gt.0) then
   k=cub(3)
   gamma=0
   do i=clb(1),cub(1)
   do j=clb(2),cub(2)
     if(psnewptr(i,j).eq.0) then
       pu=pptr(i,j,k)
       tvu=ftv(tptr(i,j,k),qptr(i,j,k))
       psnewptr(i,j)=fp0(zsnewptr(i,j),zu(i,j),pu,tvu)
     endif
   enddo
   enddo
 endif

 deallocate(zu)

 if (localpet == 0) then
!  do i=clb(1),cub(1)
!  do j=clb(2),cub(2)
   do i=clb(1),clb(1)
   do j=clb(2),clb(2)
     print*,'sfcp adjust ',(zsnewptr(i,j)-zsptr(i,j)), psptr(i,j),psnewptr(i,j)
   enddo
   enddo
 endif

 end subroutine newps

!> Reads model vertical coordinate definition file (as specified by
!! namelist variable vcoord_file_target_grid).
!!
!! @author George Gayno
 subroutine read_vcoord_info
 implicit none

 integer                    :: istat, n, k

 print*
 print*,"OPEN VERTICAL COORD FILE: ", trim(vcoord_file_target_grid)
 open(14, file=trim(vcoord_file_target_grid), form='formatted', iostat=istat)
 if (istat /= 0) then
   call error_handler("OPENING VERTICAL COORD FILE", istat)
 endif

 read(14, *, iostat=istat) nvcoord_target, lev_target
 if (istat /= 0) then
   call error_handler("READING VERTICAL COORD FILE", istat)
 endif

 levp1_target = lev_target + 1

 allocate(vcoord_target(levp1_target, nvcoord_target))
 read(14, *, iostat=istat) ((vcoord_target(n,k), k=1,nvcoord_target), n=1,levp1_target)
 if (istat /= 0) then
   call error_handler("READING VERTICAL COORD FILE", istat)
 endif

 print*
 
 close(14)

 end subroutine read_vcoord_info

!> Horizontally interpolate thompson microphysics data to the target
!! model grid.
!!
!! @author George Gayno
 subroutine horiz_interp_thomp_mp_climo

 implicit none

 integer  :: isrctermprocessing, rc

 type(esmf_regridmethod_flag)       :: method
 type(esmf_routehandle)             :: regrid_bl

 isrctermprocessing=1

 print*,"- CALL FieldCreate FOR TARGET GRID THOMP CLIMO QNIFA BEFORE ADJUSTMENT."
 qnifa_climo_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID THOMP CLIMO QNWFA BEFORE ADJUSTMENT."
 qnwfa_climo_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID THOMP CLIMO PRESSURE BEFORE ADJUSTMENT."
 thomp_pres_climo_b4adj_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID THOMP CLIMO QNIFA."
 qnifa_climo_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR TARGET GRID THOMP CLIMO QNWFA."
 qnwfa_climo_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_target/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldRegridStore FOR THOMPSON CLIMO FIELDS."

 method=ESMF_REGRIDMETHOD_BILINEAR

 call ESMF_FieldRegridStore(qnifa_climo_input_grid, &
                            qnifa_climo_b4adj_target_grid, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            srctermprocessing=isrctermprocessing, &
                            routehandle=regrid_bl, &
                            regridmethod=method, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridStore", rc)

 print*,"- CALL Field_Regrid FOR THOMP CLIMO QNIFA."
 call ESMF_FieldRegrid(qnifa_climo_input_grid, &
                       qnifa_climo_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR THOMP CLIMO QNWFA."
 call ESMF_FieldRegrid(qnwfa_climo_input_grid, &
                       qnwfa_climo_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL Field_Regrid FOR THOMP PRESSURE."
 call ESMF_FieldRegrid(thomp_pres_climo_input_grid, &
                       thomp_pres_climo_b4adj_target_grid, &
                       routehandle=regrid_bl, &
                       termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegrid", rc)

 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=regrid_bl, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)

!-----------------------------------------------------------------------------------
! Free up input data memory.
!-----------------------------------------------------------------------------------

 call cleanup_thomp_mp_climo_input_data

 end subroutine horiz_interp_thomp_mp_climo

!> Vertically interpolate atmospheric fields to target FV3 grid.
!!
!! Vertically interpolate thompson microphysics climo tracers to the
!! target model levels.
!!
!! @author George Gayno
 SUBROUTINE VINTG_THOMP_MP_CLIMO

 implicit none

 INTEGER                         :: CLB(3), CUB(3), RC
 INTEGER                         :: IM, KM1, KM2, NT
 INTEGER                         :: I, J, K

 REAL(ESMF_KIND_R8), ALLOCATABLE :: Z1(:,:,:), Z2(:,:,:)
 REAL(ESMF_KIND_R8), ALLOCATABLE :: C1(:,:,:,:),C2(:,:,:,:)

 REAL(ESMF_KIND_R8), POINTER     :: QNIFA1PTR(:,:,:)       ! input
 REAL(ESMF_KIND_R8), POINTER     :: QNIFA2PTR(:,:,:)       ! target
 REAL(ESMF_KIND_R8), POINTER     :: QNWFA1PTR(:,:,:)       ! input
 REAL(ESMF_KIND_R8), POINTER     :: QNWFA2PTR(:,:,:)       ! target
 REAL(ESMF_KIND_R8), POINTER     :: P1PTR(:,:,:)       ! input pressure
 REAL(ESMF_KIND_R8), POINTER     :: P2PTR(:,:,:)       ! target pressure

 print*,"- VERTICALY INTERPOLATE THOMP MP CLIMO TRACERS."

 print*,"- CALL FieldGet FOR 3-D THOMP PRES."
 call ESMF_FieldGet(thomp_pres_climo_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=p1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The '1'/'2' arrays hold fields before/after interpolation.  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NT=  2  ! number of thomp tracers

 ALLOCATE(Z1(CLB(1):CUB(1),CLB(2):CUB(2),lev_thomp_mp_climo))
 ALLOCATE(Z2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET))
 ALLOCATE(C1(CLB(1):CUB(1),CLB(2):CUB(2),lev_thomp_mp_climo,NT))
 ALLOCATE(C2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET,NT))

 Z1 = -LOG(P1PTR)

 print*,"- CALL FieldGet FOR 3-D ADJUSTED PRESS"
 call ESMF_FieldGet(pres_target_grid, &
                    farrayPtr=P2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 Z2 = -LOG(P2PTR)

!print*,'pres check 1 ', p1ptr(clb(1),clb(2),:)
!print*,'pres check 2 ', p2ptr(clb(1),clb(2),:)

 print*,"- CALL FieldGet FOR qnifa before vertical adjustment."
 call ESMF_FieldGet(qnifa_climo_b4adj_target_grid, &
                    farrayPtr=QNIFA1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,1) =  QNIFA1PTR(:,:,:)

 print*,"- CALL FieldGet FOR qnwfa before vertical adjustment."
 call ESMF_FieldGet(qnwfa_climo_b4adj_target_grid, &
                    farrayPtr=QNWFA1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,2) =  QNWFA1PTR(:,:,:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  PERFORM LAGRANGIAN ONE-DIMENSIONAL INTERPOLATION
!  THAT IS 4TH-ORDER IN INTERIOR, 2ND-ORDER IN OUTSIDE INTERVALS
!  AND 1ST-ORDER FOR EXTRAPOLATION.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 IM = (CUB(1)-CLB(1)+1) * (CUB(2)-CLB(2)+1)
 KM1= LEV_THOMP_MP_CLIMO
 KM2= LEV_TARGET

 CALL TERP3(IM,1,1,1,1,NT,(IM*KM1),(IM*KM2), &
            KM1,IM,IM,Z1,C1,KM2,IM,IM,Z2,C2)

 print*,"- CALL FieldGet FOR ADJUSTED climo qnifa."
 call ESMF_FieldGet(qnifa_climo_target_grid, &
                    farrayPtr=QNIFA2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ADJUSTED climo qnwfa."
 call ESMF_FieldGet(qnwfa_climo_target_grid, &
                    farrayPtr=QNWFA2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 DO K=1,LEV_TARGET
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     QNIFA2PTR(I,J,K) = C2(I,J,K,1)
     QNWFA2PTR(I,J,K) = C2(I,J,K,2)
   ENDDO
   ENDDO
 ENDDO

 DEALLOCATE (Z1, Z2, C1, C2)

 call ESMF_FieldDestroy(qnifa_climo_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(qnwfa_climo_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(thomp_pres_climo_b4adj_target_grid, rc=rc)

 END SUBROUTINE VINTG_THOMP_MP_CLIMO


!> Vertically extend model top into thermosphere for whole atmosphere model.
!!
!! Use climatological data to extent model top into thermosphere for
!! temperature and consoder primary compositions of neutral atmosphere
!! in term of specific values of oxygen, single oxygen, and ozone.
!!
!! @param [in] year  initial year
!! @param [in] month  initial month
!! @param [in] day  initial day
!! @param [in] hour  initial hour
!! @param [in] pf    path to MSIS2.1 parm file
!!
!! @author Hann-Ming Henry Juang NCEP/EMC
 SUBROUTINE VINTG_WAM (YEAR,MONTH,DAY,HOUR,PF)

 IMPLICIT NONE

 include 'mpif.h'

 INTEGER, INTENT(IN)             :: YEAR,MONTH,DAY,HOUR
 CHARACTER(*), INTENT(IN)        :: PF

 REAL(ESMF_KIND_R8), PARAMETER   :: AMO  = 15.9994  ! molecular weight of o
 REAL(ESMF_KIND_R8), PARAMETER   :: AMO2 = 31.999   !molecular weight of o2
 REAL(ESMF_KIND_R8), PARAMETER   :: AMN2 = 28.013   !molecular weight of n2

 REAL(ESMF_KIND_R8)              :: COE,WFUN(10),DEGLAT,HOLD
 REAL(ESMF_KIND_R8)              :: SUMMASS,QVMASS,O3MASS
 INTEGER                         :: I, J, K, II, CLB(3), CUB(3), RC, KREF
 INTEGER                         :: IDAT(8),JDOW,JDAY,ICDAY

 REAL(ESMF_KIND_R8), ALLOCATABLE :: TEMP(:),ON(:),O2N(:),N2N(:),PRMB(:)
        
 REAL(ESMF_KIND_R8), POINTER     :: LATPTR(:,:)        ! output latitude
 REAL(ESMF_KIND_R8), POINTER     :: P1PTR(:,:,:)       ! input pressure
 REAL(ESMF_KIND_R8), POINTER     :: P2PTR(:,:,:)       ! output pressure
 REAL(ESMF_KIND_R8), POINTER     :: DZDT2PTR(:,:,:)    ! output vvel
 REAL(ESMF_KIND_R8), POINTER     :: T2PTR(:,:,:)       ! output temperature
 REAL(ESMF_KIND_R8), POINTER     :: Q2PTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: QVPTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: QOPTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: O2PTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: O3PTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: XWIND2PTR(:,:,:)  ! output wind (x component)
 REAL(ESMF_KIND_R8), POINTER     :: YWIND2PTR(:,:,:)  ! output wind (y component)
 REAL(ESMF_KIND_R8), POINTER     :: ZWIND2PTR(:,:,:)  ! output wind (z component)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 print*,"VINTG_WAM:- VERTICALY EXTEND FIELDS FOR WAM COLD START."

! prepare date
 IDAT = 0
 JDOW = 0
 JDAY = 0
 ICDAY = 0
 IDAT(1)=year
 IDAT(2)=month
 IDAT(3)=day
 IDAT(5)=hour
 CALL W3DOXDAT(IDAT,JDOW,ICDAY,JDAY)
 print *,"VINTG_WAM: WAM START DATE FOR ICDAY=",ICDAY

! prepare weighting function
 DO K=1,10
   WFUN(K) = (K-1.0) / 9.0
 ENDDO

 ALLOCATE(TEMP(LEV_TARGET))
 ALLOCATE(PRMB(LEV_TARGET))
 ALLOCATE(  ON(LEV_TARGET))
 ALLOCATE( O2N(LEV_TARGET))
 ALLOCATE( N2N(LEV_TARGET))

! p1 (pascal)
 print*,"VINTG_WAM:- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=p1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
!print*,"VINTG_WAM: p1ptr ",(p1ptr(1,1,k),k=1,LEV_INPUT)

! p2 (pascal)
 print*,"VINTG_WAM:- CALL FieldGet FOR 3-D ADJUSTED PRESS"
 call ESMF_FieldGet(pres_target_grid, &
                    farrayPtr=P2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
!print*,"VINTG_WAM: p2ptr ",(p2ptr(1,1,k),k=1,LEV_TARGET)

! latitude in degree
 print*,"VINTG_WAM - CALL FieldGet FOR LATITUDE_S."
 call ESMF_FieldGet(latitude_s_target_grid, &
                    farrayPtr=LATPTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
!print*,"VINTG_WAM: latptr ",(latptr(1,j),j=clb(2),cub(2))

! temp
 print*,"VINTG_WAM:- CALL FieldGet FOR 3-D ADJUSTED TEMP."
 call ESMF_FieldGet(temp_target_grid, &
                    farrayPtr=T2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

! dzdt
 print*,"VINTG_WAM:- CALL FieldGet FOR ADJUSTED VERTICAL VELOCITY."
 call ESMF_FieldGet(dzdt_target_grid, &
                    farrayPtr=DZDT2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

! wind
 print*,"VINTG_WAM:- CALL FieldGet FOR ADJUSTED WIND COMPONENTS."

 call ESMF_FieldGet(xwind_target_grid, &
                    farrayPtr=XWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 call ESMF_FieldGet(ywind_target_grid, &
                    farrayPtr=YWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 call ESMF_FieldGet(zwind_target_grid, &
                    farrayPtr=ZWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

!
! determine vertical blending point and modified extrapolation values
!
 DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)

     DO K=1,LEV_TARGET
       IF(P2PTR(I,J,K).le.P1PTR(I,J,LEV_INPUT)) THEN
         KREF     =K-1
         EXIT
       ENDIF
     ENDDO
!
     DO K=KREF,LEV_TARGET
       COE = P2PTR(I,J,K) / P2PTR(I,J,KREF)
       XWIND2PTR(I,J,K) = COE*XWIND2PTR(I,J,K)
       YWIND2PTR(I,J,K) = COE*YWIND2PTR(I,J,K)
       ZWIND2PTR(I,J,K) = COE*ZWIND2PTR(I,J,K)
       DZDT2PTR(I,J,K)   = COE*DZDT2PTR(I,J,K)
     ENDDO

   ENDDO
 ENDDO

! 
! point necessary tracers
!
 DO II = 1, NUM_TRACERS

   print*,"VINTG_WAM:- CALL FieldGet FOR 3-D TRACER ", trim(tracers(ii))
   call ESMF_FieldGet(tracers_target_grid(ii), &
                      farrayPtr=Q2PTR, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", rc)

   DO J=CLB(2),CUB(2)
     DO I=CLB(1),CUB(1)
       DO K=1,LEV_TARGET
         IF(P2PTR(I,J,K).le.P1PTR(I,J,LEV_INPUT)) THEN
           KREF     =K-1
           EXIT
         ENDIF
       ENDDO
!
       DO K=KREF,LEV_TARGET
         COE = MIN(1.0, P2PTR(I,J,K) / P2PTR(I,J,KREF) )
         Q2PTR(I,J,K) = COE * Q2PTR(I,J,K)
       ENDDO
     ENDDO
   ENDDO

   IF (TRIM(TRACERS(II)) == "sphum") QVPTR => Q2PTR
   IF (TRIM(TRACERS(II)) == "spo"  ) QOPTR => Q2PTR
   IF (TRIM(TRACERS(II)) == "spo2" ) O2PTR => Q2PTR
   IF (TRIM(TRACERS(II)) == "spo3" ) O3PTR => Q2PTR

 ENDDO

!
! obtained wam gases distribution and temperature profile
!
 DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
!
     DEGLAT = LATPTR(I,J)
     DO K=1,LEV_TARGET
       PRMB(K) = P2PTR(I,J,K) * 0.01
     ENDDO
     CALL GETTEMP(ICDAY,1,DEGLAT,1,PRMB,LEV_TARGET,PF,TEMP,ON,O2N,N2N)
!
     DO K=1,LEV_TARGET
       SUMMASS = ON(K)*AMO+O2N(K)*AMO2+N2N(K)*AMN2
       QVMASS  = SUMMASS*QVPTR(I,J,K)/(1.-QVPTR(I,J,K))
       SUMMASS = SUMMASS+QVMASS
       O3MASS  = SUMMASS*O3PTR(I,J,K)
       SUMMASS = SUMMASS+O3MASS
       HOLD    = 1.0 / SUMMASS
       QOPTR(I,J,K) = ON (K)*AMO *HOLD
       O2PTR(I,J,K) = O2N(K)*AMO2*HOLD
       O3PTR(I,J,K) = O3MASS * HOLD
       QVPTR(I,J,K) = QVMASS * HOLD
     ENDDO
!
     DO K=1,LEV_TARGET
       IF(P2PTR(I,J,K).le.P1PTR(I,J,LEV_INPUT)) THEN
         KREF     =K-1
         EXIT
       ENDIF
     ENDDO
!
     DO K=KREF,LEV_TARGET
       T2PTR(I,J,K) = TEMP(K)
     ENDDO
     DO K=KREF-10,KREF-1
       T2PTR(I,J,K) = WFUN(K-KREF+11)  * TEMP(K) + &
                 (1.- WFUN(K-KREF+11)) * T2PTR(I,J,K)
     ENDDO
   ENDDO
 ENDDO

 DEALLOCATE (TEMP, PRMB, ON, O2N, N2N)

 END SUBROUTINE VINTG_WAM

!> Vertically interpolate upper-air fields.
!!
!! Vertically interpolate upper-air fields. Wind, temperature,
!! humidity and other tracers are interpolated. The interpolation is
!! cubic lagrangian in log pressure with a monotonic constraint in the
!! center of the domain. In the outer intervals it is linear in log
!! pressure. Outside the domain, fields are generally held constant,
!! except for temperature and humidity below the input domain, where
!! the temperature lapse rate is held fixed at -6.5 k/km and the
!! relative humidity is held constant. This routine expects fields
!! ordered from bottom to top of atmosphere.
!!
!! @author Mark Iredell @date 92-10-31
 SUBROUTINE VINTG
 use mpi

 IMPLICIT NONE

 REAL(ESMF_KIND_R8), PARAMETER   :: DLTDZ=-6.5E-3*287.05/9.80665
 REAL(ESMF_KIND_R8), PARAMETER   :: DLPVDRT=-2.5E6/461.50
 REAL(ESMF_KIND_R8), PARAMETER   :: ONE = 1.0_ESMF_KIND_R8

 INTEGER                         :: I, J, K, CLB(3), CUB(3), RC
 INTEGER                         :: IM, KM1, KM2, NT, II

 REAL(ESMF_KIND_R8)              :: DZ
 REAL(ESMF_KIND_R8), ALLOCATABLE :: Z1(:,:,:), Z2(:,:,:)
 REAL(ESMF_KIND_R8), ALLOCATABLE :: C1(:,:,:,:),C2(:,:,:,:)
        
 REAL(ESMF_KIND_R8), POINTER     :: P1PTR(:,:,:)       ! input pressure
 REAL(ESMF_KIND_R8), POINTER     :: P2PTR(:,:,:)       ! output pressure
 REAL(ESMF_KIND_R8), POINTER     :: DZDT1PTR(:,:,:)    ! input vvel
 REAL(ESMF_KIND_R8), POINTER     :: DZDT2PTR(:,:,:)    ! output vvel
 REAL(ESMF_KIND_R8), POINTER     :: T1PTR(:,:,:)       ! input temperature
 REAL(ESMF_KIND_R8), POINTER     :: T2PTR(:,:,:)       ! output temperature
 REAL(ESMF_KIND_R8), POINTER     :: Q1PTR(:,:,:)       ! input tracer
 REAL(ESMF_KIND_R8), POINTER     :: Q2PTR(:,:,:)       ! output tracer
 REAL(ESMF_KIND_R8), POINTER     :: XWIND1PTR(:,:,:)  ! input wind (x component)
 REAL(ESMF_KIND_R8), POINTER     :: YWIND1PTR(:,:,:)  ! input wind (y component)
 REAL(ESMF_KIND_R8), POINTER     :: ZWIND1PTR(:,:,:)  ! input wind (z component)
 REAL(ESMF_KIND_R8), POINTER     :: XWIND2PTR(:,:,:)  ! output wind (x component)
 REAL(ESMF_KIND_R8), POINTER     :: YWIND2PTR(:,:,:)  ! output wind (y component)
 REAL(ESMF_KIND_R8), POINTER     :: ZWIND2PTR(:,:,:)  ! output wind (z component)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE LOG PRESSURE INTERPOLATING COORDINATE
!  AND COPY INPUT WIND, TEMPERATURE, HUMIDITY AND OTHER TRACERS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 print*,"- VERTICALY INTERPOLATE FIELDS."

 print*,"- CALL FieldGet FOR 3-D PRES."
 call ESMF_FieldGet(pres_b4adj_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=p1ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The '1'/'2' arrays hold fields before/after interpolation.  
! Note the 'z' component of the horizontal wind will be treated as a
! tracer.  So add one extra third dimension to these 3-d arrays.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ALLOCATE(Z1(CLB(1):CUB(1),CLB(2):CUB(2),LEV_INPUT))
 ALLOCATE(Z2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET))
 ALLOCATE(C1(CLB(1):CUB(1),CLB(2):CUB(2),LEV_INPUT,NUM_TRACERS_INPUT+5))
 ALLOCATE(C2(CLB(1):CUB(1),CLB(2):CUB(2),LEV_TARGET,NUM_TRACERS_INPUT+5))

 Z1 = -LOG(P1PTR)

 print*,"- CALL FieldGet FOR 3-D ADJUSTED PRESS"
 call ESMF_FieldGet(pres_target_grid, &
                    farrayPtr=P2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 Z2 = -LOG(P2PTR)
 
 print*,"- CALL FieldGet FOR x WIND."
 call ESMF_FieldGet(xwind_b4adj_target_grid, &
                    farrayPtr=XWIND1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,1) =  XWIND1PTR(:,:,:)

 print*,"- CALL FieldGet FOR y WIND."
 call ESMF_FieldGet(ywind_b4adj_target_grid, &
                    farrayPtr=YWIND1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,2) =  YWIND1PTR(:,:,:)

 print*,"- CALL FieldGet FOR z WIND."
 call ESMF_FieldGet(zwind_b4adj_target_grid, &
                    farrayPtr=ZWIND1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,3) =  ZWIND1PTR(:,:,:)

 print*,"- CALL FieldGet FOR VERTICAL VELOCITY."
 call ESMF_FieldGet(dzdt_b4adj_target_grid, &
                    farrayPtr=DZDT1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,4) =  DZDT1PTR(:,:,:)
 print*,"MIN MAX W TARGETB4 IN VINTG = ", minval(DZDT1PTR(:,:,:)), maxval(DZDT1PTR(:,:,:))

 print*,"- CALL FieldGet FOR 3-D TEMP."
 call ESMF_FieldGet(temp_b4adj_target_grid, &
                    farrayPtr=T1PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 C1(:,:,:,5) =  T1PTR(:,:,:)

 DO I = 1, NUM_TRACERS_INPUT

   print*,"- CALL FieldGet FOR 3-D TRACERS ", trim(tracers(i))
   call ESMF_FieldGet(tracers_b4adj_target_grid(i), &
                      farrayPtr=Q1PTR, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", rc)

   C1(:,:,:,5+I) =  Q1PTR(:,:,:)

 ENDDO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  PERFORM LAGRANGIAN ONE-DIMENSIONAL INTERPOLATION
!  THAT IS 4TH-ORDER IN INTERIOR, 2ND-ORDER IN OUTSIDE INTERVALS
!  AND 1ST-ORDER FOR EXTRAPOLATION.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 IM = (CUB(1)-CLB(1)+1) * (CUB(2)-CLB(2)+1)
 KM1= LEV_INPUT
 KM2= LEV_TARGET
 NT=  NUM_TRACERS_INPUT + 1 ! treat 'z' wind as tracer.

 CALL TERP3(IM,1,1,1,1,4+NT,(IM*KM1),(IM*KM2), &
            KM1,IM,IM,Z1,C1,KM2,IM,IM,Z2,C2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COPY OUTPUT WIND, TEMPERATURE, HUMIDITY AND OTHER TRACERS
!  EXCEPT BELOW THE INPUT DOMAIN, LET TEMPERATURE INCREASE WITH A FIXED
!  LAPSE RATE AND LET THE RELATIVE HUMIDITY REMAIN CONSTANT.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 print*,"- CALL FieldGet FOR 3-D ADJUSTED TEMP."
 call ESMF_FieldGet(temp_target_grid, &
                    farrayPtr=T2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ADJUSTED VERTICAL VELOCITY."
 call ESMF_FieldGet(dzdt_target_grid, &
                    farrayPtr=DZDT2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ADJUSTED xwind."
 call ESMF_FieldGet(xwind_target_grid, &
                    farrayPtr=XWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ADJUSTED ywind."
 call ESMF_FieldGet(ywind_target_grid, &
                    farrayPtr=YWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR ADJUSTED zwind."
 call ESMF_FieldGet(zwind_target_grid, &
                    farrayPtr=ZWIND2PTR, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 DO K=1,LEV_TARGET
   DO I=CLB(1),CUB(1)
   DO J=CLB(2),CUB(2)
     XWIND2PTR(I,J,K)=C2(I,J,K,1)
     YWIND2PTR(I,J,K)=C2(I,J,K,2)
     ZWIND2PTR(I,J,K)=C2(I,J,K,3)
     DZDT2PTR(I,J,K)=C2(I,J,K,4)
     DZ=Z2(I,J,K)-Z1(I,J,1)
     IF(DZ.GE.0) THEN
       T2PTR(I,J,K)=C2(I,J,K,5)
     ELSE
       T2PTR(I,J,K)=C1(I,J,1,5)*EXP(DLTDZ*DZ)
     ENDIF
   ENDDO
   ENDDO
 ENDDO

 DO II = 1, NUM_TRACERS_INPUT

   print*,"- CALL FieldGet FOR 3-D TRACER ", trim(tracers(ii))
   call ESMF_FieldGet(tracers_target_grid(ii), &
                      farrayPtr=Q2PTR, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", rc)

   IF (TRIM(TRACERS(II)) == "sphum") THEN  ! specific humidity

     DO K=1,LEV_TARGET
       DO I=CLB(1),CUB(1)
       DO J=CLB(2),CUB(2)
         DZ=Z2(I,J,K)-Z1(I,J,1)
         IF(DZ.GE.0) THEN
           Q2PTR(I,J,K) = C2(I,J,K,5+II)
         ELSE
           Q2PTR(I,J,K) = C1(I,J,1,5+II)*EXP(DLPVDRT*(ONE/T2PTR(I,J,K)-ONE/T1PTR(I,J,1))-DZ)
         ENDIF
       ENDDO
       ENDDO
     ENDDO

   ELSE ! all other tracers

     DO K=1,LEV_TARGET
       DO I=CLB(1),CUB(1)
       DO J=CLB(2),CUB(2)
         Q2PTR(I,J,K) = C2(I,J,K,5+II)
       ENDDO
       ENDDO
     ENDDO

   ENDIF

 ENDDO

 DEALLOCATE (Z1, Z2, C1, C2)

 END SUBROUTINE VINTG

!> Cubically interpolate in one dimension.
!!                                                                       
!! Interpolate field(s) in one dimension along the column(s). The
!! interpolation is cubic lagrangian with a monotonic constraint in
!! the center of the domain. In the outer intervals it is linear.
!! Outside the domain, fields are held constant.
!!                                                                       
!! PROGRAM HISTORY LOG:                                                  
!! -  98-05-01  MARK IREDELL                                              
!! - 1999-01-04  IREDELL  USE ESSL SEARCH                                  
!!                                                                       
!! @param[in] im integer number of columns                            
!! @param[in] ixz1 integer column skip number for z1                    
!! @param[in] ixq1 integer column skip number for q1                    
!! @param[in] ixz2 integer column skip number for z2                    
!! @param[in] ixq2 integer column skip number for q2                    
!! @param[in] nm integer number of fields per column                  
!! @param[in] nxq1 integer field skip number for q1                     
!! @param[in] nxq2 integer field skip number for q2                     
!! @param[in] km1 integer number of input points                       
!! @param[in] kxz1 integer point skip number for z1                     
!! @param[in] kxq1 integer point skip number for q1                     
!! @param[in] z1 real (1+(im-1)*ixz1+(km1-1)*kxz1)                    
!!                  input coordinate values in which to interpolate      
!!                  (z1 must be strictly monotonic in either direction)  
!! @param[in] q1 real (1+(im-1)*ixq1+(km1-1)*kxq1+(nm-1)*nxq1)        
!!                  input fields to interpolate                          
!! @param[in] km2 integer number of output points                      
!! @param[in] kxz2 integer point skip number for z2                     
!! @param[in] kxq2 integer point skip number for q2                     
!! @param[in] z2 real (1+(im-1)*ixz2+(km2-1)*kxz2)                    
!!                  output coordinate values to which to interpolate     
!!                  (z2 need not be monotonic)                           
!! @param[out] q2 real (1+(im-1)*ixq2+(km2-1)*kxq2+(nm-1)*nxq2)        
!!                  output interpolated fields                           
!! @author Mark Iredell @date 98-05-01            
 SUBROUTINE TERP3(IM,IXZ1,IXQ1,IXZ2,IXQ2,NM,NXQ1,NXQ2,             &
                  KM1,KXZ1,KXQ1,Z1,Q1,KM2,KXZ2,KXQ2,Z2,Q2)      
      IMPLICIT NONE 
      INTEGER IM,IXZ1,IXQ1,IXZ2,IXQ2,NM,NXQ1,NXQ2 
      INTEGER KM1,KXZ1,KXQ1,KM2,KXZ2,KXQ2 
      INTEGER I,K1,K2,N 
      INTEGER K1S(IM,KM2) 
      REAL(ESMF_KIND_R8), PARAMETER :: ONE = 1.0_ESMF_KIND_R8
      REAL(ESMF_KIND_R8) :: Z1(1+(IM-1)*IXZ1+(KM1-1)*KXZ1) 
      REAL(ESMF_KIND_R8) :: Q1(1+(IM-1)*IXQ1+(KM1-1)*KXQ1+(NM-1)*NXQ1) 
      REAL(ESMF_KIND_R8) :: Z2(1+(IM-1)*IXZ2+(KM2-1)*KXZ2) 
      REAL(ESMF_KIND_R8) :: Q2(1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2) 
!     REAL(ESMF_KIND_R8) :: J2(1+(IM-1)*IXQ2+(KM2-1)*KXQ2+(NM-1)*NXQ2) 
      REAL(ESMF_KIND_R8) :: FFA(IM),FFB(IM),FFC(IM),FFD(IM) 
      REAL(ESMF_KIND_R8) :: GGA(IM),GGB(IM),GGC(IM),GGD(IM) 
      REAL(ESMF_KIND_R8) :: Z1A,Z1B,Z1C,Z1D,Q1A,Q1B,Q1C,Q1D,Z2S,Q2S
!     REAL(ESMF_KIND_R8) :: J2S 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  FIND THE SURROUNDING INPUT INTERVAL FOR EACH OUTPUT POINT.
      CALL RSEARCH(IM,KM1,IXZ1,KXZ1,Z1,KM2,IXZ2,KXZ2,Z2,1,IM,K1S) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  GENERALLY INTERPOLATE CUBICALLY WITH MONOTONIC CONSTRAINT            
!  FROM TWO NEAREST INPUT POINTS ON EITHER SIDE OF THE OUTPUT POINT,    
!  BUT WITHIN THE TWO EDGE INTERVALS INTERPOLATE LINEARLY.              
!  KEEP THE OUTPUT FIELDS CONSTANT OUTSIDE THE INPUT DOMAIN.            
                                                                        
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(IM,IXZ1,IXQ1,IXZ2), &
!$OMP& SHARED(IXQ2,NM,NXQ1,NXQ2,KM1,KXZ1,KXQ1,Z1,Q1,KM2,KXZ2), &
!$OMP& SHARED(KXQ2,Z2,Q2,K1S)
      DO K2=1,KM2 
        DO I=1,IM 
          K1=K1S(I,K2) 
          IF(K1.EQ.1.OR.K1.EQ.KM1-1) THEN 
            Z2S=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            Z1A=Z1(1+(I-1)*IXZ1+(K1-1)*KXZ1) 
            Z1B=Z1(1+(I-1)*IXZ1+(K1+0)*KXZ1) 
            FFA(I)=(Z2S-Z1B)/(Z1A-Z1B) 
            FFB(I)=(Z2S-Z1A)/(Z1B-Z1A) 
            GGA(I)=ONE/(Z1A-Z1B) 
            GGB(I)=ONE/(Z1B-Z1A) 
          ELSEIF(K1.GT.1.AND.K1.LT.KM1-1) THEN 
            Z2S=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
            Z1A=Z1(1+(I-1)*IXZ1+(K1-2)*KXZ1) 
            Z1B=Z1(1+(I-1)*IXZ1+(K1-1)*KXZ1) 
            Z1C=Z1(1+(I-1)*IXZ1+(K1+0)*KXZ1) 
            Z1D=Z1(1+(I-1)*IXZ1+(K1+1)*KXZ1) 
            FFA(I)=(Z2S-Z1B)/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)                                  
            FFB(I)=(Z2S-Z1A)/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)                                  
            FFC(I)=(Z2S-Z1A)/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)                                  
            FFD(I)=(Z2S-Z1A)/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)                                  
            GGA(I)=      ONE/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)+                                 &
                   (Z2S-Z1B)/(Z1A-Z1B)*                                 &
                         ONE/(Z1A-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1A-Z1D)+                                 &
                   (Z2S-Z1B)/(Z1A-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1A-Z1C)*                                 &
                         ONE/(Z1A-Z1D)                                  
            GGB(I)=      ONE/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1B-Z1A)*                                 &
                         ONE/(Z1B-Z1C)*                                 &
                   (Z2S-Z1D)/(Z1B-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1B-Z1A)*                                 &
                   (Z2S-Z1C)/(Z1B-Z1C)*                                 &
                         ONE/(Z1B-Z1D)                                  
            GGC(I)=      ONE/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1C-Z1A)*                                 &
                         ONE/(Z1C-Z1B)*                                 &
                   (Z2S-Z1D)/(Z1C-Z1D)+                                 &
                   (Z2S-Z1A)/(Z1C-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1C-Z1B)*                                 &
                         ONE/(Z1C-Z1D)                                  
            GGD(I)=      ONE/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)+                                 &
                   (Z2S-Z1A)/(Z1D-Z1A)*                                 &
                         ONE/(Z1D-Z1B)*                                 &
                   (Z2S-Z1C)/(Z1D-Z1C)+                                 &
                   (Z2S-Z1A)/(Z1D-Z1A)*                                 &
                   (Z2S-Z1B)/(Z1D-Z1B)*                                 &
                         ONE/(Z1D-Z1C)                                  
          ENDIF 
        ENDDO 

!  INTERPOLATE.                                                         
        DO N=1,NM 
          DO I=1,IM 
            K1=K1S(I,K2) 
            IF(K1.EQ.0) THEN 
              Q2S=Q1(1+(I-1)*IXQ1+(N-1)*NXQ1) 
!             J2S=0 
            ELSEIF(K1.EQ.KM1) THEN 
              Q2S=Q1(1+(I-1)*IXQ1+(KM1-1)*KXQ1+(N-1)*NXQ1) 
!             J2S=0 
            ELSEIF(K1.EQ.1.OR.K1.EQ.KM1-1) THEN 
              Q1A=Q1(1+(I-1)*IXQ1+(K1-1)*KXQ1+(N-1)*NXQ1) 
              Q1B=Q1(1+(I-1)*IXQ1+(K1+0)*KXQ1+(N-1)*NXQ1) 
              Q2S=FFA(I)*Q1A+FFB(I)*Q1B 
!             J2S=GGA(I)*Q1A+GGB(I)*Q1B 
            ELSE 
              Q1A=Q1(1+(I-1)*IXQ1+(K1-2)*KXQ1+(N-1)*NXQ1) 
              Q1B=Q1(1+(I-1)*IXQ1+(K1-1)*KXQ1+(N-1)*NXQ1) 
              Q1C=Q1(1+(I-1)*IXQ1+(K1+0)*KXQ1+(N-1)*NXQ1) 
              Q1D=Q1(1+(I-1)*IXQ1+(K1+1)*KXQ1+(N-1)*NXQ1) 
              Q2S=FFA(I)*Q1A+FFB(I)*Q1B+FFC(I)*Q1C+FFD(I)*Q1D 
!             J2S=GGA(I)*Q1A+GGB(I)*Q1B+GGC(I)*Q1C+GGD(I)*Q1D 
              IF(Q2S.LT.MIN(Q1B,Q1C)) THEN 
                Q2S=MIN(Q1B,Q1C) 
!               J2S=0 
              ELSEIF(Q2S.GT.MAX(Q1B,Q1C)) THEN 
                Q2S=MAX(Q1B,Q1C) 
!               J2S=0 
              ENDIF 
            ENDIF 
            Q2(1+(I-1)*IXQ2+(K2-1)*KXQ2+(N-1)*NXQ2)=Q2S 
!           J2(1+(I-1)*IXQ2+(K2-1)*KXQ2+(N-1)*NXQ2)=J2S 
          ENDDO 
        ENDDO 
      ENDDO 
!$OMP END PARALLEL DO                                                   

 END SUBROUTINE TERP3 

!> Search for a surrounding real interval.
!!                                                                       
!! This subprogram searches monotonic sequences of real numbers for
!! intervals that surround a given search set of real numbers. The
!! sequences may be monotonic in either direction; the real numbers
!! may be single or double precision; the input sequences and sets and
!! the output locations may be arbitrarily dimensioned.
!!                                                                       
!! If the array z1 is dimensioned (im,km1), then the skip numbers are
!! ixz1=1 and kxz1=im; if it is dimensioned (km1,im), then the skip
!! numbers are ixz1=km1 and kxz1=1; if it is dimensioned (im,jm,km1),
!! then the skip numbers are ixz1=1 and kxz1=im*jm; etcetera. Similar
!! examples apply to the skip numbers for z2 and l2.
!!                                                                       
!! Returned values of 0 or km1 indicate that the given search value    
!! is outside the range of the sequence. 
!!                                                                       
!! If a search value is identical to one of the sequence values then
!! the location returned points to the identical value. If the
!! sequence is not strictly monotonic and a search value is identical
!! to more than one of the sequence values, then the location returned
!! may point to any of the identical values.
!!                                                                       
!! to be exact, for each i from 1 to im and for each k from 1 to km2,
!! z=z2(1+(i-1)*ixz2+(k-1)*kxz2) is the search value and
!! l=l2(1+(i-1)*ixl2+(k-1)*kxl2) is the location returned.  if l=0,
!! then z is less than the start point z1(1+(i-1)*ixz1) for ascending
!! sequences (or greater than for descending sequences).  if l=km1,
!! then z is greater than or equal to the end point
!! z1(1+(i-1)*ixz1+(km1-1)*kxz1) for ascending sequences (or less than
!! or equal to for descending sequences).  otherwise z is between the
!! values z1(1+(i-1)*ixz1+(l-1)*kxz1) and z1(1+(i-1)*ixz1+(l-0)*kxz1)
!! and may equal the former.
!!                                                                       
!! @param[in] im integer number of sequences to search                
!! @param[in] km1 integer number of points in each sequence            
!! @param[in] ixz1 integer sequence skip number for z1                  
!! @param[in] kxz1 integer point skip number for z1                     
!! @param[in] z1 real (1+(im-1)*ixz1+(km1-1)*kxz1)                    
!!                  sequence values to search                            
!!                  (z1 must be monotonic in either direction)           
!! @param[in] km2 integer number of points to search for               
!!                  in each respective sequence                          
!! @param[in] ixz2 integer sequence skip number for z2                  
!! @param[in] kxz2 integer point skip number for z2                     
!! @param[in] z2 real (1+(im-1)*ixz2+(km2-1)*kxz2)                    
!!                  set of values to search for                          
!!                  (z2 need not be monotonic)                           
!! @param[in] ixl2 integer sequence skip number for l2                  
!! @param[in] kxl2 integer point skip number for l2                     
!!                                                                       
!! @param[out] l2 integer (1+(im-1)*ixl2+(km2-1)*kxl2)                 
!!                  interval locations having values from 0 to km1       
!!                  (z2 will be between z1(l2) and z1(l2+1))             
!!                                                                       
!! @author Mark Iredell @date 98-05-01            
 SUBROUTINE RSEARCH(IM,KM1,IXZ1,KXZ1,Z1,KM2,IXZ2,KXZ2,Z2,IXL2,KXL2,L2)
 IMPLICIT NONE 

 INTEGER,INTENT(IN)    :: IM,KM1,IXZ1,KXZ1,KM2,IXZ2,KXZ2,IXL2,KXL2 
 INTEGER,INTENT(OUT)   :: L2(1+(IM-1)*IXL2+(KM2-1)*KXL2) 

 REAL(ESMF_KIND_R8),INTENT(IN) :: Z1(1+(IM-1)*IXZ1+(KM1-1)*KXZ1) 
 REAL(ESMF_KIND_R8),INTENT(IN) :: Z2(1+(IM-1)*IXZ2+(KM2-1)*KXZ2) 

 INTEGER                       :: I,K2,L

 REAL(ESMF_KIND_R8)            :: Z 

  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  FIND THE SURROUNDING INPUT INTERVAL FOR EACH OUTPUT POINT.          
 DO I=1,IM 
   IF (Z1(1+(I-1)*IXZ1).LE.Z1(1+(I-1)*IXZ1+(KM1-1)*KXZ1)) THEN 
!  INPUT COORDINATE IS MONOTONICALLY ASCENDING.                        
     DO K2=1,KM2
       Z=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2)
       L=0 
       DO 
         IF(Z.LT.Z1(1+(I-1)*IXZ1+L*KXZ1)) EXIT 
         L=L+1 
         IF(L.EQ.KM1) EXIT 
       ENDDO
       L2(1+(I-1)*IXL2+(K2-1)*KXL2)=L 
     ENDDO 
   ELSE 
!   INPUT COORDINATE IS MONOTONICALLY DESCENDING.                       
     DO K2=1,KM2 
       Z=Z2(1+(I-1)*IXZ2+(K2-1)*KXZ2) 
       L=0 
       DO 
         IF(Z.GT.Z1(1+(I-1)*IXZ1+L*KXZ1)) EXIT 
         L=L+1 
         IF(L.EQ.KM1) EXIT 
       ENDDO
       L2(1+(I-1)*IXL2+(K2-1)*KXL2)=L 
     ENDDO 
   ENDIF 
 ENDDO 
                                                                        
 END SUBROUTINE RSEARCH 

!> Compute vertical level height
!! @author George Gayno
 subroutine compute_zh

 implicit none 

 integer                          :: i,ii, j,k, rc, clb(2), cub(2)

 real(esmf_kind_r8), allocatable  :: pe0(:), pn0(:)
 real(esmf_kind_r8), pointer      :: psptr(:,:)
 real(esmf_kind_r8), pointer      :: zhsfcptr(:,:)
 real(esmf_kind_r8), pointer      :: zhptr(:,:,:)
 real(esmf_kind_r8), pointer      :: tptr(:,:,:)
 real(esmf_kind_r8), pointer      :: qptr(:,:,:)
 real(esmf_kind_r8)               :: ak, bk, zvir, grd
 real(esmf_kind_r8), parameter    :: grav  = 9.80665 
 real(esmf_kind_r8), parameter    :: rdgas = 287.05 
 real(esmf_kind_r8), parameter    :: rvgas = 461.50 

 print*,"- COMPUTE HEIGHT"

 print*,"- CALL FieldGet FOR SURFACE PRESSURE"
 call ESMF_FieldGet(ps_target_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TERRAIN HEIGHT"
 call ESMF_FieldGet(terrain_target_grid, &
                    farrayPtr=zhsfcptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR HEIGHT"
 call ESMF_FieldGet(zh_target_grid, &
                    farrayPtr=zhptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR TEMPERATURE"
 call ESMF_FieldGet(temp_target_grid, &
                    farrayPtr=tptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 do ii = 1, num_tracers
   if (trim(tracers(ii)) == "sphum") exit
 enddo

 print*,"- CALL FieldGet FOR SPECIFIC HUMIDITY"
 call ESMF_FieldGet(tracers_target_grid(ii), &
                    farrayPtr=qptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

 grd = grav/rdgas 
 zvir = rvgas/rdgas - 1.0_esmf_kind_r8
                                                                        
 allocate(pe0(levp1_target))
 allocate(pn0(levp1_target))

 do j = clb(2), cub(2)
 do i = clb(1), cub(1)

   do k = 1, levp1_target
     ak = vcoord_target(k,1)
     ak = max(ak, 1.e-9)
     bk = vcoord_target(k,2)

     pe0(k) = ak + bk*psptr(i,j)
     pn0(k) = log(pe0(k))
   enddo

   zhptr(i,j,1) = zhsfcptr(i,j)

   do k = 2, levp1_target
     zhptr(i,j,k) = zhptr(i,j,k-1)+tptr(i,j,k-1)*(1.+zvir*qptr(i,j,k-1))*     &
              (pn0(k-1)-pn0(k))/grd
   enddo

 enddo
 enddo

 deallocate(pe0, pn0)

 end subroutine compute_zh 
 
!> Cleanup atmospheric field (before adjustment) objects.
!!
!! @author George Gayno
 subroutine cleanup_target_atm_b4adj_data

 implicit none

 integer                     :: i, rc

 print*,"- DESTROY TARGET GRID ATMOSPHERIC BEFORE ADJUSTMENT FIELDS."

 call ESMF_FieldDestroy(xwind_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(ywind_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(zwind_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(dzdt_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(ps_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(pres_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(temp_b4adj_target_grid, rc=rc)
 call ESMF_FieldDestroy(terrain_interp_to_target_grid, rc=rc)

 do i = 1, num_tracers_input
   call ESMF_FieldDestroy(tracers_b4adj_target_grid(i), rc=rc)
 enddo

 deallocate(tracers_b4adj_target_grid)

 end subroutine cleanup_target_atm_b4adj_data

!> Cleanup target grid atmospheric field objects.
!! @author George Gayno
 subroutine cleanup_all_target_atm_data

 use atmosphere_target_data, only : cleanup_atmosphere_target_data

 implicit none

 integer                     :: rc

 print*,"- DESTROY LOCAL TARGET GRID ATMOSPHERIC FIELDS."

 call ESMF_FieldDestroy(pres_target_grid, rc=rc)
 call ESMF_FieldDestroy(xwind_target_grid, rc=rc)
 call ESMF_FieldDestroy(ywind_target_grid, rc=rc)
 call ESMF_FieldDestroy(zwind_target_grid, rc=rc)
 call ESMF_FieldDestroy(xwind_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(ywind_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(zwind_s_target_grid, rc=rc)
 call ESMF_FieldDestroy(xwind_w_target_grid, rc=rc)
 call ESMF_FieldDestroy(ywind_w_target_grid, rc=rc)
 call ESMF_FieldDestroy(zwind_w_target_grid, rc=rc)

 call cleanup_atmosphere_target_data

 end subroutine cleanup_all_target_atm_data

 end module atmosphere
