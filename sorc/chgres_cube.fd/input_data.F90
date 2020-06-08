 module input_data

!--------------------------------------------------------------------------
! Module input_data
!
! Abstract: Read atmospheric, surface and nst data on the input grid.
!    Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!    'history' files, fv3 gaussian history files, spectral gfs
!    gaussian nemsio files, and spectral gfs sigio/sfcio files.
!
! Public Subroutines:
! -----------------
! read_input_atm_data         Driver routine to read atmospheric data
! cleanup_input_atm_data      Free up memory associated with atm data
! read_input_sfc_data         Driver routine to read surface data
! cleanup_input_sfc_data      Free up memory associated with sfc data
! read_input_nst_data         Driver routine to read nst data
! cleanup_input_nst_data      Free up memory associated with nst data
!
! Public variables:
! -----------------
! Defined below.  "input" indicates field associated with the input grid.
!
!--------------------------------------------------------------------------

 use esmf
 use netcdf
 use nemsio_module

 use program_setup, only          : data_dir_input_grid, &
                                    nst_files_input_grid, &
                                    sfc_files_input_grid, &
                                    atm_files_input_grid, &
                                    grib2_file_input_grid, &
                                    atm_core_files_input_grid, &
                                    atm_tracer_files_input_grid, &
                                    convert_nst, &
                                    orog_dir_input_grid, &
                                    orog_files_input_grid, &
                                    tracers_input, num_tracers, &
                                    input_type, tracers, &
                                    get_var_cond, read_from_input

 use model_grid, only             : input_grid,        &
                                    i_input, j_input,  &
                                    ip1_input, jp1_input,  &
                                    num_tiles_input_grid, &
                                    latitude_input_grid, &
                                    longitude_input_grid, &
                                    inv_file

 implicit none

 private

! Fields associated with the atmospheric model.

 type(esmf_field), public              :: dzdt_input_grid       ! vert velocity
 type(esmf_field)                      :: dpres_input_grid      ! pressure thickness
 type(esmf_field), public              :: pres_input_grid       ! 3-d pressure
 type(esmf_field), public              :: ps_input_grid         ! surface pressure
 type(esmf_field), public              :: terrain_input_grid    ! terrain height
 type(esmf_field), public              :: temp_input_grid       ! temperature
 type(esmf_field)                      :: u_input_grid          ! u/v wind at grid
 type(esmf_field)                      :: v_input_grid          ! box center
 type(esmf_field), public              :: wind_input_grid       ! 3-component wind
 type(esmf_field), allocatable, public :: tracers_input_grid(:) ! tracers

 integer, public                 :: lev_input      ! # of atmospheric layers
 integer, public                 :: levp1_input    ! # of atmos layer interfaces

! Fields associated with the land-surface model.

 integer, public                 :: veg_type_landice_input = 15 ! NOAH land ice option
                                                                ! defined at this veg type.
                                                                ! Default is igbp.

 type(esmf_field), public        :: canopy_mc_input_grid    ! canopy moist content
 type(esmf_field), public        :: f10m_input_grid         ! log((z0+10)*1/z0)
 type(esmf_field), public        :: ffmm_input_grid         ! log((z0+z1)*1/z0)
                                                            ! See sfc_diff.f for details.
 type(esmf_field), public        :: landsea_mask_input_grid ! land sea mask;
                                                            ! 0-water, 1-land, 2-ice
 type(esmf_field), public        :: q2m_input_grid          ! 2-m spec hum
 type(esmf_field), public        :: seaice_depth_input_grid ! sea ice depth
 type(esmf_field), public        :: seaice_fract_input_grid ! sea ice fraction
 type(esmf_field), public        :: seaice_skin_temp_input_grid  ! sea ice skin temp
 type(esmf_field), public        :: skin_temp_input_grid    ! skin temp/sst
 type(esmf_field), public        :: snow_depth_input_grid   ! snow dpeth
 type(esmf_field), public        :: snow_liq_equiv_input_grid ! snow liq equiv depth
 type(esmf_field), public        :: soil_temp_input_grid    ! 3-d soil temp
 type(esmf_field), public        :: soil_type_input_grid    ! soil type
 type(esmf_field), public        :: soilm_liq_input_grid    ! 3-d liquid soil moisture
 type(esmf_field), public        :: soilm_tot_input_grid    ! 3-d total soil moisture
 type(esmf_field), public        :: srflag_input_grid       ! snow/rain flag
 type(esmf_field), public        :: t2m_input_grid          ! 2-m temperature
 type(esmf_field), public        :: tprcp_input_grid        ! precip
 type(esmf_field), public        :: ustar_input_grid        ! fric velocity
 type(esmf_field), public        :: veg_type_input_grid     ! vegetation type
 type(esmf_field), public        :: z0_input_grid           ! roughness length

 integer, public      :: lsoil_input=4  ! # of soil layers,
                                                   ! # hardwire for now
 
 character(len=50), private, allocatable :: slevs(:)                           

! Fields associated with the nst model.

 type(esmf_field), public        :: c_d_input_grid
 type(esmf_field), public        :: c_0_input_grid
 type(esmf_field), public        :: d_conv_input_grid
 type(esmf_field), public        :: dt_cool_input_grid
 type(esmf_field), public        :: ifd_input_grid
 type(esmf_field), public        :: qrain_input_grid
 type(esmf_field), public        :: tref_input_grid  ! reference temperature
 type(esmf_field), public        :: w_d_input_grid
 type(esmf_field), public        :: w_0_input_grid
 type(esmf_field), public        :: xs_input_grid
 type(esmf_field), public        :: xt_input_grid
 type(esmf_field), public        :: xu_input_grid
 type(esmf_field), public        :: xv_input_grid
 type(esmf_field), public        :: xz_input_grid
 type(esmf_field), public        :: xtts_input_grid
 type(esmf_field), public        :: xzts_input_grid
 type(esmf_field), public        :: z_c_input_grid
 type(esmf_field), public        :: zm_input_grid

 public :: read_input_atm_data
 public :: cleanup_input_atm_data
 public :: read_input_sfc_data
 public :: cleanup_input_sfc_data
 public :: read_input_nst_data
 public :: cleanup_input_nst_data
 
 contains

!---------------------------------------------------------------------------
! Read input grid atmospheric data driver
!---------------------------------------------------------------------------

 subroutine read_input_atm_data(localpet)

 implicit none

 integer, intent(in)             :: localpet

!-------------------------------------------------------------------------------
! Read the tiled 'warm' restart files.
!-------------------------------------------------------------------------------

 if (trim(input_type) == "restart") then

   call read_input_atm_restart_file(localpet)

!-------------------------------------------------------------------------------
! Read the gaussian history files in netcdf format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gaussian_netcdf") then

   call read_input_atm_gaussian_netcdf_file(localpet)

!-------------------------------------------------------------------------------
! Read the tiled history files in netcdf format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "history") then

   call read_input_atm_tiled_history_file(localpet)

!-------------------------------------------------------------------------------
! Read the gaussian history files in nemsio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gaussian_nemsio") then  ! fv3gfs gaussian nemsio

   call read_input_atm_gaussian_nemsio_file(localpet)

!-------------------------------------------------------------------------------
! Read the spectral gfs gaussian history files in nemsio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gfs_gaussian_nemsio") then ! spectral gfs gaussian 
                                                         ! nemsio.
   call read_input_atm_gfs_gaussian_nemsio_file(localpet)

!-------------------------------------------------------------------------------
! Read the spectral gfs gaussian history files in sigio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gfs_sigio") then ! spectral gfs sigio format.
 
   call read_input_atm_gfs_sigio_file(localpet)

!-------------------------------------------------------------------------------
! Read fv3gfs data in grib2 format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "grib2") then

   call read_input_atm_grib2_file(localpet)

 endif

 end subroutine read_input_atm_data

!---------------------------------------------------------------------------
! Read input grid nst data driver
!---------------------------------------------------------------------------

 subroutine read_input_nst_data(localpet)

 implicit none

 integer, intent(in)             :: localpet

 integer                         :: rc

 print*,"- READ INPUT GRID NST DATA."

 print*,"- CALL FieldCreate FOR INPUT GRID C_D."
 c_d_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID C_0."
 c_0_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID D_CONV."
 d_conv_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID DT_COOL."
 dt_cool_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID IFD."
 ifd_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID QRAIN."
 qrain_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID TREF."
 tref_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID W_D."
 w_d_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID W_0."
 w_0_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XS."
 xs_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XT."
 xt_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XU."
 xu_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XV."
 xv_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XZ."
 xz_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XTTS."
 xtts_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID XZTS."
 xzts_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID Z_C."
 z_c_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID ZM."
 zm_input_grid = ESMF_FieldCreate(input_grid, &
                                  typekind=ESMF_TYPEKIND_R8, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN FieldCreate", rc)
 
!--------------------------------------------------------------------------
! Read input grid nst data from a fv3 gaussian nemsio history file or
! spectral GFS nemsio file.
!--------------------------------------------------------------------------

 if (trim(input_type) == "gaussian_nemsio" .or. trim(input_type) == "gfs_gaussian_nemsio") then

   call read_input_nst_nemsio_file(localpet)

!---------------------------------------------------------------------------
! Read nst data from these netcdf formatted fv3 files: tiled history,
! tiled warm restart, and gaussian history.
!---------------------------------------------------------------------------

 else

   call read_input_nst_netcdf_file(localpet)

 endif

 end subroutine read_input_nst_data

!---------------------------------------------------------------------------
! Read input grid surface data driver.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_data(localpet)

 implicit none

 integer, intent(in)             :: localpet

 integer                         :: rc

 print*,"- CALL FieldCreate FOR INPUT GRID LANDSEA MASK."
 landsea_mask_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID Z0."
 z0_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID VEGETATION TYPE."
 veg_type_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID CANOPY MOISTURE CONTENT."
 canopy_mc_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SEAICE FRACTION."
 seaice_fract_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SEAICE DEPTH."
 seaice_depth_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SEAICE SKIN TEMPERATURE."
 seaice_skin_temp_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SNOW DEPTH."
 snow_depth_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SNOW LIQUID EQUIVALENT."
 snow_liq_equiv_input_grid = ESMF_FieldCreate(input_grid, &
                                     typekind=ESMF_TYPEKIND_R8, &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID T2M."
 t2m_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID Q2M."
 q2m_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID TPRCP."
 tprcp_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID F10M."
 f10m_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID USTAR."
 ustar_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID FFMM."
 ffmm_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SRFLAG."
 srflag_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT SKIN TEMPERATURE."
 skin_temp_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT SOIL TYPE."
 soil_type_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT TERRAIN."
 terrain_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT SOIL TEMPERATURE."
 soil_temp_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT TOTAL SOIL MOISTURE."
 soilm_tot_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT LIQUID SOIL MOISTURE."
 soilm_liq_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lsoil_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

!-------------------------------------------------------------------------------
! Read the tiled 'warm' restart files.
!-------------------------------------------------------------------------------

 if (trim(input_type) == "restart") then

   call read_input_sfc_restart_file(localpet)

!-------------------------------------------------------------------------------
! Read the tiled or gaussian history files in netcdf format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "history" .or. trim(input_type) ==  &
         "gaussian_netcdf") then

   call read_input_sfc_netcdf_file(localpet)

!-------------------------------------------------------------------------------
! Read the gaussian history files in nemsio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gaussian_nemsio") then

   call read_input_sfc_gaussian_nemsio_file(localpet)

!-------------------------------------------------------------------------------
! Read the spectral gfs gaussian history files in nemsio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gfs_gaussian_nemsio") then

   call read_input_sfc_gfs_gaussian_nemsio_file(localpet)

!-------------------------------------------------------------------------------
! Read the spectral gfs gaussian history files in sfcio format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "gfs_sigio") then

   call read_input_sfc_gfs_sfcio_file(localpet)

!-------------------------------------------------------------------------------
! Read fv3gfs surface data in grib2 format.
!-------------------------------------------------------------------------------

 elseif (trim(input_type) == "grib2") then

   call read_input_sfc_grib2_file(localpet)

 endif

 end subroutine read_input_sfc_data

!---------------------------------------------------------------------------
! Create atmospheric esmf fields.
!---------------------------------------------------------------------------

 subroutine init_atm_esmf_fields
 
 implicit none

 integer                                  :: i, rc

 print*,"- INITIALIZE ATMOSPHERIC ESMF FIELDS."

 print*,"- CALL FieldCreate FOR INPUT GRID 3-D WIND."
 wind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1,1/), &
                                   ungriddedUBound=(/lev_input,3/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID SURFACE PRESSURE."
 ps_input_grid = ESMF_FieldCreate(input_grid, &
                                  typekind=ESMF_TYPEKIND_R8, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID TERRAIN."
 terrain_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID TEMPERATURE."
 temp_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 allocate(tracers_input_grid(num_tracers))

 do i = 1, num_tracers
   print*,"- CALL FieldCreate FOR INPUT GRID TRACER ", trim(tracers_input(i))
   tracers_input_grid(i) = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldCreate", rc)
 enddo

 print*,"- CALL FieldCreate FOR INPUT GRID DZDT."
 dzdt_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID U."
 u_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID V."
 v_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID PRESSURE."
 pres_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)
 
 end subroutine init_atm_esmf_fields

!---------------------------------------------------------------------------
! Read input atmospheric data from spectral gfs (old sigio format).
! Used prior to July 19, 2017.
!---------------------------------------------------------------------------

 subroutine read_input_atm_gfs_sigio_file(localpet)

 use sigio_module

 implicit none

 integer, intent(in)                   :: localpet

 character(len=300)                    :: the_file

 integer(sigio_intkind)                :: iret
 integer                               :: rc, i, j, k
 integer                               :: clb(3), cub(3)

 real(esmf_kind_r8)                    :: ak, bk
 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d2(:,:,:)
 real(esmf_kind_r8), pointer           :: pptr(:,:,:), psptr(:,:)
 real(esmf_kind_r8), allocatable       :: pi(:,:,:)

 type(sigio_head)                      :: sighead
 type(sigio_dbta)                      :: sigdata

 the_file = trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))

 print*,"- ATMOSPHERIC DATA IN SIGIO FORMAT."
 print*,"- OPEN AND READ: ", trim(the_file)

 call sigio_sropen(21, trim(the_file), iret)
 if (iret /= 0) then
   rc = iret
   call error_handler("OPENING SPECTRAL GFS SIGIO FILE.", rc)
 endif
 call sigio_srhead(21, sighead, iret)
 if (iret /= 0) then
   rc = iret
   call error_handler("READING SPECTRAL GFS SIGIO FILE.", rc)
 endif

 lev_input = sighead%levs
 levp1_input = lev_input + 1

 if (num_tracers /= sighead%ntrac) then
   call error_handler("WRONG NUMBER OF TRACERS EXPECTED.", 99)
 endif

 if (sighead%idvt == 0 .or. sighead%idvt == 21) then
   if (trim(tracers_input(1)) /= 'spfh'  .or.  &
       trim(tracers_input(2)) /= 'o3mr'   .or.  &
       trim(tracers_input(3)) /= 'clwmr') then 
     call error_handler("TRACERS SELECTED DO NOT MATCH FILE CONTENTS.", 99)
   endif
 else
   print*,'- UNRECOGNIZED IDVT: ', sighead%idvt
   call error_handler("UNRECOGNIZED IDVT", 99)
 endif

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 if (localpet == 0) then
   allocate(dummy2d(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lev_input))
   allocate(dummy3d2(i_input,j_input,lev_input))
 else
   allocate(dummy2d(0,0))
   allocate(dummy3d(0,0,0))
   allocate(dummy3d2(0,0,0))
 endif

 if (localpet == 0) then
   call sigio_aldbta(sighead, sigdata, iret)
   if (iret /= 0) then
     rc = iret
     call error_handler("ALLOCATING SIGDATA.", rc)
   endif
   call sigio_srdbta(21, sighead, sigdata, iret)
   if (iret /= 0) then
     rc = iret
     call error_handler("READING SIGDATA.", rc)
   endif
   call sptez(0,sighead%jcap,4,i_input, j_input, sigdata%ps, dummy2d, 1)
   dummy2d = exp(dummy2d) * 1000.0
   print*,'surface pres ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR SURFACE PRESSURE."
 call ESMF_FieldScatter(ps_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   call sptez(0,sighead%jcap,4,i_input, j_input, sigdata%hs, dummy2d, 1)
   print*,'terrain ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 do k = 1, num_tracers

   if (localpet == 0) then
     call sptezm(0,sighead%jcap,4,i_input, j_input, lev_input, sigdata%q(:,:,k), dummy3d, 1)
     print*,trim(tracers_input(k)),maxval(dummy3d),minval(dummy3d)
   endif

   print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input(k))
   call ESMF_FieldScatter(tracers_input_grid(k), dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo

 if (localpet == 0) then
   call sptezm(0,sighead%jcap,4,i_input, j_input, lev_input, sigdata%t, dummy3d, 1)
   print*,'temp ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
 call ESMF_FieldScatter(temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

!---------------------------------------------------------------------------
! The spectral gfs files have omega, not vertical velocity.  Set to
! zero for now.  Convert from omega to vv in the future?
!---------------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- NO VERTICAL VELOCITY RECORD.  SET TO ZERO."
   dummy3d = 0.0
 endif

 print*,"- CALL FieldScatter FOR INPUT DZDT."
 call ESMF_FieldScatter(dzdt_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   call sptezmv(0, sighead%jcap, 4, i_input, j_input, lev_input, sigdata%d, sigdata%z, dummy3d, dummy3d2, 1)
   print*,'u ',maxval(dummy3d),minval(dummy3d)
   print*,'v ',maxval(dummy3d2),minval(dummy3d2)
 endif

 print*,"- CALL FieldScatter FOR INPUT U-WIND."
 call ESMF_FieldScatter(u_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL FieldScatter FOR INPUT V-WIND."
 call ESMF_FieldScatter(v_input_grid, dummy3d2, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d, dummy3d, dummy3d2)

 if (localpet == 0) call sigio_axdbta(sigdata, iret)

 call sigio_sclose(21, iret)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d component winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute 3-d pressure from 'ak' and 'bk'.
!---------------------------------------------------------------------------

 print*,"- COMPUTE 3-D PRESSURE."

 print*,"- CALL FieldGet FOR 3-D PRES."
 nullify(pptr)
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 nullify(psptr)
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!---------------------------------------------------------------------------
! First, compute interface pressure.
!---------------------------------------------------------------------------

 allocate(pi(clb(1):cub(1),clb(2):cub(2),1:levp1_input),stat=rc)

 do k=1,levp1_input
   ak = sighead%vcoord(k,1)
   bk = sighead%vcoord(k,2)
   do i= clb(1), cub(1)
     do j= clb(2), cub(2)
       pi(i,j,k) = ak + bk*psptr(i,j)
     enddo
   enddo
 enddo

 if (localpet == 0) then
   print*,'pres int ',psptr(clb(1),clb(2)),pi(clb(1),clb(2),:)
 endif

!---------------------------------------------------------------------------
! Now comput mid-layer pressure from interface pressure.
!---------------------------------------------------------------------------

 do k=1,lev_input
   do i= clb(1), cub(1)
     do j= clb(2), cub(2)
       pptr(i,j,k) = (pi(i,j,k)+pi(i,j,k+1))/2.0_esmf_kind_r8
     enddo
   enddo
 enddo

 deallocate(pi)

 if (localpet == 0) then
   print*,'pres ',psptr(clb(1),clb(2)),pptr(clb(1),clb(2),:)
 endif

 end subroutine read_input_atm_gfs_sigio_file

!---------------------------------------------------------------------------
! Read input atmospheric data from spectral gfs (global gaussian in
! nemsio format. Starting July 19, 2017).  
!---------------------------------------------------------------------------

 subroutine read_input_atm_gfs_gaussian_nemsio_file(localpet)

 implicit none

 integer, intent(in)                   :: localpet

 character(len=300)                    :: the_file
 character(len=20)                     :: vlevtyp, vname

 integer(nemsio_intkind)               :: vlev, iret
 integer                               :: i, j, k, n, rc
 integer                               :: clb(3), cub(3)

 real(nemsio_realkind), allocatable    :: vcoord(:,:,:)
 real(nemsio_realkind), allocatable    :: dummy(:)
 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)
 real(esmf_kind_r8)                    :: ak, bk
 real(esmf_kind_r8), allocatable       :: pi(:,:,:)
 real(esmf_kind_r8), pointer           :: pptr(:,:,:), psptr(:,:)

 type(nemsio_gfile)                    :: gfile

 the_file = trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))

 print*,"- READ ATMOS DATA FROM SPECTRAL GFS NEMSIO FILE: ", trim(the_file)

 print*,"- OPEN FILE."
 call nemsio_open(gfile, the_file, "read", iret=iret)
 if (iret /= 0) call error_handler("OPENING SPECTRAL GFS NEMSIO ATM FILE.", iret)

 print*,"- READ NUMBER OF VERTICAL LEVELS."
 call nemsio_getfilehead(gfile, iret=iret, dimz=lev_input)
 if (iret /= 0) call error_handler("READING NUMBER OF VERTICAL LEVLES.", iret)

 levp1_input = lev_input + 1

 allocate(vcoord(levp1_input,3,2))

 print*,"- READ VERTICAL COORDINATE INFO."
 call nemsio_getfilehead(gfile, iret=iret, vcoord=vcoord)
 if (iret /= 0) call error_handler("READING VERTICAL COORDINATE INFO.", iret)

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 if (localpet == 0) then
   allocate(dummy(i_input*j_input))
   allocate(dummy2d(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lev_input))
 else
   allocate(dummy(0))
   allocate(dummy2d(0,0))
   allocate(dummy3d(0,0,0))
 endif

!-----------------------------------------------------------------------
! 3-d fields in gaussian files increment from bottom to model top.
! That is what is expected by this program, so no need to flip indices.
!-----------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- READ TEMPERATURE."
   vname = "tmp"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING TEMPERATURE RECORD.", iret)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
!    print*,'temp check after read ',vlev, dummy3d(1,1,vlev)
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
 call ESMF_FieldScatter(temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 do n = 1, num_tracers

   if (localpet == 0) then
     print*,"- READ ", trim(tracers_input(n))
     vname = trim(tracers_input(n))
     vlevtyp = "mid layer"
     do vlev = 1, lev_input
       call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
       if (iret /= 0) call error_handler("READING TRACER RECORD.", iret)
!      print*,'tracer ',vlev, maxval(dummy),minval(dummy)
       dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
     enddo
   endif

   print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input(n))
   call ESMF_FieldScatter(tracers_input_grid(n), dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo

 if (localpet == 0) then
   print*,"- READ U-WINDS."
   vname = "ugrd"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING U-WIND RECORD.", iret)
!    print*,'ugrd ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT U-WIND."
 call ESMF_FieldScatter(u_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ V-WINDS."
   vname = "vgrd"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING V-WIND RECORD.", iret)
!    print*,'vgrd ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT V-WIND."
 call ESMF_FieldScatter(v_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

!---------------------------------------------------------------------------
! The spectral gfs nemsio files do not have a vertical velocity or
! omega record.  So set to zero for now.
!---------------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- NO VERTICAL VELOCITY RECORD.  SET TO ZERO."
   dummy3d = 0.0
 endif

 print*,"- CALL FieldScatter FOR INPUT DZDT."
 call ESMF_FieldScatter(dzdt_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ HGT."
   vname = "hgt"
   vlevtyp = "sfc"
   vlev = 1
   call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
   if (iret /= 0) call error_handler("READING HGT RECORD.", iret)
!  print*,'hgt ',vlev, maxval(dummy),minval(dummy)
   dummy2d = reshape(dummy, (/i_input,j_input/))
 endif

 print*,"- CALL FieldScatter FOR TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ PRES."
   vname = "pres"
   vlevtyp = "sfc"
   vlev = 1
   call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
   if (iret /= 0) call error_handler("READING PRES RECORD.", iret)
!  print*,'pres ',vlev, maxval(dummy),minval(dummy)
   dummy2d = reshape(dummy, (/i_input,j_input/))
 endif

 print*,"- CALL FieldScatter FOR SURFACE PRESSURE."
 call ESMF_FieldScatter(ps_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 call nemsio_close(gfile)

 deallocate(dummy, dummy2d, dummy3d)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d component winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute 3-d pressure from 'ak' and 'bk'.
!---------------------------------------------------------------------------

 print*,"- COMPUTE 3-D PRESSURE."

 print*,"- CALL FieldGet FOR 3-D PRES."
 nullify(pptr)
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=pptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 nullify(psptr)
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

!---------------------------------------------------------------------------
! First, compute interface pressure.
!---------------------------------------------------------------------------

 allocate(pi(clb(1):cub(1),clb(2):cub(2),1:levp1_input))

 do k=1,levp1_input
   ak = vcoord(k,1,1)
   bk = vcoord(k,2,1)
   do i= clb(1), cub(1)
     do j= clb(2), cub(2)
       pi(i,j,k) = ak + bk*psptr(i,j)
     enddo
   enddo
 enddo

 deallocate(vcoord)

!---------------------------------------------------------------------------
! Now comput mid-layer pressure from interface pressure.
!---------------------------------------------------------------------------

 do k=1,lev_input
   do i= clb(1), cub(1)
     do j= clb(2), cub(2)
       pptr(i,j,k) = (pi(i,j,k)+pi(i,j,k+1))/2.0
     enddo
   enddo
 enddo

 deallocate(pi)

 end subroutine read_input_atm_gfs_gaussian_nemsio_file

!---------------------------------------------------------------------------
! Read input grid atmospheric fv3 gaussian nemsio files.
!---------------------------------------------------------------------------

 subroutine read_input_atm_gaussian_nemsio_file(localpet)

 implicit none

 integer, intent(in)                   :: localpet

 character(len=300)                    :: the_file
 character(len=20)                     :: vlevtyp, vname

 integer                               :: i, j, k, n
 integer                               :: rc, clb(3), cub(3)
 integer(nemsio_intkind)               :: vlev, iret

 real(nemsio_realkind), allocatable    :: vcoord(:,:,:)
 real(nemsio_realkind), allocatable    :: dummy(:)
 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)
 real(esmf_kind_r8), pointer           :: presptr(:,:,:), psptr(:,:)
 real(esmf_kind_r8), pointer           :: dpresptr(:,:,:)
 real(esmf_kind_r8), allocatable       :: pres_interface(:)

 type(nemsio_gfile)                    :: gfile

 the_file = trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))

 print*,"- READ ATMOS DATA FROM GAUSSIAN NEMSIO FILE: ", trim(the_file)

 print*,"- OPEN FILE."
 call nemsio_open(gfile, the_file, "read", iret=iret)
 if (iret /= 0) call error_handler("OPENING GAUSSIAN NEMSIO ATM FILE.", iret)

 print*,"- READ NUMBER OF VERTICAL LEVELS."
 call nemsio_getfilehead(gfile, iret=iret, dimz=lev_input)
 if (iret /= 0) call error_handler("READING NUMBER OF VERTICAL LEVLES.", iret)

 levp1_input = lev_input + 1

 allocate(vcoord(levp1_input,3,2))

 print*,"- READ VERTICAL COORDINATE INFO."
 call nemsio_getfilehead(gfile, iret=iret, vcoord=vcoord)
 if (iret /= 0) call error_handler("READING VERTICAL COORDINATE INFO.", iret)

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 print*,"- CALL FieldCreate FOR INPUT DPRES."
 dpres_input_grid = ESMF_FieldCreate(input_grid, &
                                 typekind=ESMF_TYPEKIND_R8, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 ungriddedLBound=(/1/), &
                                 ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (localpet == 0) then
   allocate(dummy(i_input*j_input))
   allocate(dummy2d(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lev_input))
 else
   allocate(dummy(0))
   allocate(dummy2d(0,0))
   allocate(dummy3d(0,0,0))
 endif

!-----------------------------------------------------------------------
! 3-d fields in gaussian files increment from bottom to model top.
! That is what is expected by this program, so no need to flip indices.
!-----------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- READ TEMPERATURE."
   vname = "tmp"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING TEMPERATURE RECORD.", iret)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
     print*,'temp check after read ',vlev, dummy3d(1,1,vlev)
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
 call ESMF_FieldScatter(temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 do n = 1, num_tracers

   if (localpet == 0) then
     print*,"- READ ", trim(tracers_input(n))
     vname = trim(tracers_input(n))
     vlevtyp = "mid layer"
     do vlev = 1, lev_input
       call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
       if (iret /= 0) call error_handler("READING TRACER RECORD.", iret)
       print*,'tracer ',vlev, maxval(dummy),minval(dummy)
       dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
     enddo
   endif

   print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input(n))
   call ESMF_FieldScatter(tracers_input_grid(n), dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo

 if (localpet == 0) then
   print*,"- READ U-WINDS."
   vname = "ugrd"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING U-WIND RECORD.", iret)
     print*,'ugrd ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT U-WIND."
 call ESMF_FieldScatter(u_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ V-WINDS."
   vname = "vgrd"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING V-WIND RECORD.", iret)
     print*,'vgrd ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT V-WIND."
 call ESMF_FieldScatter(v_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ DPRES."
   vname = "dpres"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING DPRES RECORD.", iret)
     print*,'dpres ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT DPRES."
 call ESMF_FieldScatter(dpres_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ DZDT."
   vname = "dzdt"
   vlevtyp = "mid layer"
   do vlev = 1, lev_input
     call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
     if (iret /= 0) call error_handler("READING DZDT RECORD.", iret)
     print*,'dzdt ',vlev, maxval(dummy),minval(dummy)
     dummy3d(:,:,vlev) = reshape(dummy, (/i_input,j_input/))
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT DZDT."
 call ESMF_FieldScatter(dzdt_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ HGT."
   vname = "hgt"
   vlevtyp = "sfc"
   vlev = 1
   call nemsio_readrecv(gfile, vname, vlevtyp, vlev, dummy, 0, iret)
   if (iret /= 0) call error_handler("READING HGT RECORD.", iret)
   print*,'hgt ',vlev, maxval(dummy),minval(dummy)
   dummy2d = reshape(dummy, (/i_input,j_input/))
 endif

 print*,"- CALL FieldScatter FOR TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 call nemsio_close(gfile)

 deallocate(dummy, dummy2d, dummy3d)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d component winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute 3-d pressure.  Mid-layer and surface pressure are computed
! from delta p.  The surface pressure in the file is not used.  After
! the model's write component interpolates from the cubed-sphere grid
! to the gaussian grid, the surface pressure is no longer consistent
! with the delta p (per Jun Wang).
!---------------------------------------------------------------------------

 print*,"- COMPUTE 3-D PRESSURE."

 print*,"- CALL FieldGet FOR DELTA PRESSURE."
 nullify(dpresptr)
 call ESMF_FieldGet(dpres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=dpresptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR 3-D PRESSURE."
 nullify(presptr)
 call ESMF_FieldGet(pres_input_grid, &
                    farrayPtr=presptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 nullify(psptr)
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(pres_interface(levp1_input))

 if (localpet == 0) then
   do k = clb(3), cub(3)
   print*,'dpres is ',cub(1),cub(2),k, dpresptr(cub(1),cub(2),k)
   enddo
 endif

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     pres_interface(levp1_input) = vcoord(levp1_input,1,1)
     do k = lev_input, 1, -1
       pres_interface(k) = pres_interface(k+1) + dpresptr(i,j,k)
     enddo
     psptr(i,j) = pres_interface(1)
     do k = 1, lev_input
       presptr(i,j,k) = (pres_interface(k) + pres_interface(k+1)) / 2.0_8
     enddo
   enddo
 enddo

 deallocate(vcoord)

 if (localpet == 0) then
   print*,'psfc is ',clb(1),clb(2),psptr(clb(1),clb(2))
   print*,'pres is ',clb(1),clb(2),presptr(clb(1),clb(2),:)
 endif

 print*,'pres check 1',localpet,maxval(presptr(:,:,1)),minval(presptr(:,:,1))
 print*,'pres check lev',localpet,maxval(presptr(:,:,lev_input)),minval(presptr(:,:,lev_input))

 deallocate(pres_interface)

 call ESMF_FieldDestroy(dpres_input_grid, rc=rc)

 end subroutine read_input_atm_gaussian_nemsio_file

!---------------------------------------------------------------------------
! Read input grid fv3 atmospheric data 'warm' restart files.
!
! Routine reads tiled files in parallel.  Tile 1 is read by 
! localpet 0; tile 2 by localpet 1, etc.  The number of pets
! must be equal to or greater than the number of tiled files.  
! Logic only tested with global input data of six tiles.
!---------------------------------------------------------------------------

 subroutine read_input_atm_restart_file(localpet)

 implicit none

 integer, intent(in)             :: localpet

 character(len=500)              :: tilefile

 integer                         :: i, j, k
 integer                         :: clb(3), cub(3)
 integer                         :: rc, tile, ncid, id_var
 integer                         :: error, id_dim

 real(esmf_kind_r8), allocatable :: ak(:)
 real(esmf_kind_r8), pointer     :: presptr(:,:,:), psptr(:,:)
 real(esmf_kind_r8), pointer     :: dpresptr(:,:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile_3d(:,:,:)
 real(esmf_kind_r8), allocatable :: pres_interface(:)

!---------------------------------------------------------------------------
! Get number of vertical levels and model top pressure.
!---------------------------------------------------------------------------

 tilefile = trim(data_dir_input_grid) // "/" // trim(atm_core_files_input_grid(7))
 print*,"- READ ATM VERTICAL LEVELS FROM: ", trim(tilefile)
 error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(tilefile) )

 error=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
 call netcdf_err(error, 'reading xaxis_1 id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=levp1_input)
 call netcdf_err(error, 'reading xaxis_1 value' )

 lev_input = levp1_input - 1

 allocate(ak(levp1_input))

 error=nf90_inq_varid(ncid, 'ak', id_var)
 call netcdf_err(error, 'reading field id' )
 error=nf90_get_var(ncid, id_var, ak)
 call netcdf_err(error, 'reading ak' )

 error = nf90_close(ncid)

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 print*,"- CALL FieldCreate FOR INPUT GRID DELTA PRESSURE."
 dpres_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (localpet < num_tiles_input_grid) then
   allocate(data_one_tile_3d(i_input,j_input,lev_input))
   allocate(data_one_tile(i_input,j_input))
 else
   allocate(data_one_tile_3d(0,0,0))
   allocate(data_one_tile(0,0))
 endif

 if (localpet < num_tiles_input_grid) then
   tile = localpet+1
   tilefile= trim(data_dir_input_grid) // "/" // trim(atm_core_files_input_grid(tile))
   print*,"- READ ATMOSPHERIC CORE FILE: ", trim(tilefile)
   error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
   call netcdf_err(error, 'opening: '//trim(tilefile) )
 endif

 if (localpet < num_tiles_input_grid) then
   error=nf90_inq_varid(ncid, 'phis', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile)
   call netcdf_err(error, 'reading field' )
   data_one_tile = data_one_tile / 9.806_8  ! geopotential height
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN for tile ",tile
   call ESMF_FieldScatter(terrain_input_grid, data_one_tile, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
!  error=nf90_inq_varid(ncid, 'W', id_var)
!  call netcdf_err(error, 'reading field id' )
!  error=nf90_get_var(ncid, id_var, data_one_tile_3d)
!  call netcdf_err(error, 'reading field' )
!  data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)

! Using 'w' from restart files has caused problems.  Set to zero.
   data_one_tile_3d = 0.0_8
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID VERTICAL VELOCITY for tile ",tile
   call ESMF_FieldScatter(dzdt_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   error=nf90_inq_varid(ncid, 'T', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
   call ESMF_FieldScatter(temp_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   error=nf90_inq_varid(ncid, 'delp', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT DELTA PRESSURE."
   call ESMF_FieldScatter(dpres_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   error=nf90_inq_varid(ncid, 'ua', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID U."
   call ESMF_FieldScatter(u_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   error=nf90_inq_varid(ncid, 'va', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID V."
   call ESMF_FieldScatter(v_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid)  error = nf90_close(ncid)

 if (localpet < num_tiles_input_grid) then
   tile = localpet+1
   tilefile= trim(data_dir_input_grid) // "/" // trim(atm_tracer_files_input_grid(tile))
   print*,"- READ ATMOSPHERIC TRACER FILE: ", trim(tilefile)
   error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
   call netcdf_err(error, 'opening: '//trim(tilefile) )
 endif

 do i = 1, num_tracers

   if (localpet < num_tiles_input_grid) then
     error=nf90_inq_varid(ncid, tracers_input(i), id_var)
     call netcdf_err(error, 'reading field id' )
     error=nf90_get_var(ncid, id_var, data_one_tile_3d)
     call netcdf_err(error, 'reading field' )
     data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
   endif

   do tile = 1, num_tiles_input_grid
     print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input(i))
     call ESMF_FieldScatter(tracers_input_grid(i), data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)
   enddo

 enddo

 if (localpet < num_tiles_input_grid) error=nf90_close(ncid)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d cartesian winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute pressures
!---------------------------------------------------------------------------

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR PRESSURE."
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=presptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR DELTA PRESSURE."
 call ESMF_FieldGet(dpres_input_grid, &
                    farrayPtr=dpresptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(pres_interface(levp1_input))

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     pres_interface(levp1_input) = ak(1)  ! model top in Pa
     do k = (levp1_input-1), 1, -1
       pres_interface(k) = pres_interface(k+1) + dpresptr(i,j,k)
     enddo
     do k = 1, lev_input
       presptr(i,j,k) = (pres_interface(k) + pres_interface(k+1)) / 2.0_8
     enddo
     psptr(i,j) = pres_interface(1)
   enddo
 enddo

 deallocate(ak)
 deallocate(pres_interface)

 call ESMF_FieldDestroy(dpres_input_grid, rc=rc)

 deallocate(data_one_tile_3d, data_one_tile)

 end subroutine read_input_atm_restart_file

!---------------------------------------------------------------------------
! Read fv3 netcdf gaussian history file.  Each task reads a horizontal
! slice.
!---------------------------------------------------------------------------

 subroutine read_input_atm_gaussian_netcdf_file(localpet)

 implicit none

 include 'mpif.h'

 integer, intent(in)               :: localpet

 character(len=500)                :: tilefile

 integer                           :: start(3), count(3), iscnt
 integer                           :: error, ncid, num_tracers_file
 integer                           :: id_dim, idim_input, jdim_input
 integer                           :: id_var, rc, nprocs, max_procs
 integer                           :: kdim, remainder, myrank, i, j, k, n
 integer                           :: clb(3), cub(3)
 integer, allocatable              :: kcount(:), startk(:), displ(:)
 integer, allocatable              :: ircnt(:)

 real(esmf_kind_r8), allocatable   :: phalf(:)
 real(esmf_kind_r8), allocatable   :: pres_interface(:)
 real(kind=4), allocatable         :: dummy3d(:,:,:)
 real(kind=4), allocatable         :: dummy3dall(:,:,:)
 real(esmf_kind_r8), allocatable   :: dummy3dflip(:,:,:)
 real(esmf_kind_r8), allocatable   :: dummy(:,:)
 real(esmf_kind_r8), pointer       :: presptr(:,:,:), dpresptr(:,:,:)
 real(esmf_kind_r8), pointer       :: psptr(:,:)

 print*,"- READ INPUT ATMOS DATA FROM GAUSSIAN NETCDF FILE."

 tilefile = trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))
 error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(tilefile) )

 error=nf90_inq_dimid(ncid, 'grid_xt', id_dim)
 call netcdf_err(error, 'reading grid_xt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=idim_input)
 call netcdf_err(error, 'reading grid_xt value' )

 error=nf90_inq_dimid(ncid, 'grid_yt', id_dim)
 call netcdf_err(error, 'reading grid_yt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=jdim_input)
 call netcdf_err(error, 'reading grid_yt value' )

 if (idim_input /= i_input .or. jdim_input /= j_input) then
   call error_handler("DIMENSION MISMATCH BETWEEN SFC AND OROG FILES.", 2)
 endif

 error=nf90_inq_dimid(ncid, 'pfull', id_dim)
 call netcdf_err(error, 'reading pfull id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=lev_input)
 call netcdf_err(error, 'reading pfull value' )

 error=nf90_inq_dimid(ncid, 'phalf', id_dim)
 call netcdf_err(error, 'reading phalf id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=levp1_input)
 call netcdf_err(error, 'reading phalf value' )
 allocate(phalf(levp1_input))
 error=nf90_inq_varid(ncid, 'phalf', id_var)
 call netcdf_err(error, 'getting phalf varid' )
 error=nf90_get_var(ncid, id_var, phalf)
 call netcdf_err(error, 'reading phalf varid' )

 error=nf90_get_att(ncid, nf90_global, 'ncnsto', num_tracers_file)
 call netcdf_err(error, 'reading ntracer value' )

 call mpi_comm_size(mpi_comm_world, nprocs, error)
 print*,'- Running with ', nprocs, ' processors'

 call mpi_comm_rank(mpi_comm_world, myrank, error)
 print*,'- myrank/localpet is ',myrank,localpet

 max_procs = nprocs
 if (nprocs > lev_input) then
   max_procs = lev_input
 endif

 kdim = lev_input / max_procs
 remainder = lev_input - (max_procs*kdim)

 allocate(kcount(0:nprocs-1))
 kcount=0
 allocate(startk(0:nprocs-1))
 startk=0
 allocate(displ(0:nprocs-1))
 displ=0
 allocate(ircnt(0:nprocs-1))
 ircnt=0

 do k = 0, max_procs-2
   kcount(k) = kdim
 enddo
 kcount(max_procs-1) = kdim + remainder

 startk(0) = 1
 do k = 1, max_procs-1
   startk(k) = startk(k-1) + kcount(k-1)
 enddo

 ircnt(:) = idim_input * jdim_input * kcount(:)

 displ(0) = 0
 do k = 1, max_procs-1
   displ(k) = displ(k-1) + ircnt(k-1)
 enddo

 iscnt=idim_input*jdim_input*kcount(myrank)

! Account for case if number of tasks exceeds the number of vert levels.

 if (myrank <= max_procs-1) then
   allocate(dummy3d(idim_input,jdim_input,kcount(myrank)))
 else
   allocate(dummy3d(0,0,0))
 endif

 if (myrank == 0) then
  allocate(dummy3dall(idim_input,jdim_input,lev_input))
  dummy3dall = 0.0
  allocate(dummy3dflip(idim_input,jdim_input,lev_input))
  dummy3dflip = 0.0
  allocate(dummy(idim_input,jdim_input))
  dummy = 0.0
 else
  allocate(dummy3dall(0,0,0))
  allocate(dummy3dflip(0,0,0))
  allocate(dummy(0,0))
 endif

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 print*,"- CALL FieldCreate FOR INPUT GRID DELTA PRESSURE."
 dpres_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

! Temperature 

 if (myrank <= max_procs-1) then
   start = (/1,1,startk(myrank)/)
   count = (/idim_input,jdim_input,kcount(myrank)/)
   error=nf90_inq_varid(ncid, 'tmp', id_var)
   call netcdf_err(error, 'reading tmp field id' )
   error=nf90_get_var(ncid, id_var, dummy3d, start=start, count=count)
   call netcdf_err(error, 'reading tmp field' )
 endif

 call mpi_gatherv(dummy3d, iscnt, mpi_real, &
                  dummy3dall, ircnt, displ, mpi_real, &
                  0, mpi_comm_world, error)
 if (error /= 0) call error_handler("IN mpi_gatherv of temperature", error)

 if (myrank == 0) then
   dummy3dflip(:,:,1:lev_input) = dummy3dall(:,:,lev_input:1:-1)
 endif
   
 print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE "
 call ESMF_FieldScatter(temp_input_grid, dummy3dflip, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

! dpres

 if (myrank <= max_procs-1) then
   error=nf90_inq_varid(ncid, 'dpres', id_var)
   call netcdf_err(error, 'reading dpres field id' )
   error=nf90_get_var(ncid, id_var, dummy3d, start=start, count=count)
   call netcdf_err(error, 'reading dpres field' )
 endif

 call mpi_gatherv(dummy3d, iscnt, mpi_real, &
                  dummy3dall, ircnt, displ, mpi_real, &
                  0, mpi_comm_world, error)
 if (error /= 0) call error_handler("IN mpi_gatherv of dpres", error)

 if (myrank == 0) then
   dummy3dflip(:,:,1:lev_input) = dummy3dall(:,:,lev_input:1:-1)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID DPRES "
 call ESMF_FieldScatter(dpres_input_grid, dummy3dflip, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

! ugrd

 if (myrank <= max_procs-1) then
   error=nf90_inq_varid(ncid, 'ugrd', id_var)
   call netcdf_err(error, 'reading ugrd field id' )
   error=nf90_get_var(ncid, id_var, dummy3d, start=start, count=count)
   call netcdf_err(error, 'reading ugrd field' )
 endif

 call mpi_gatherv(dummy3d, iscnt, mpi_real, &
                  dummy3dall, ircnt, displ, mpi_real, &
                  0, mpi_comm_world, error)
 if (error /= 0) call error_handler("IN mpi_gatherv of ugrd", error)

 if (myrank == 0) then
   dummy3dflip(:,:,1:lev_input) = dummy3dall(:,:,lev_input:1:-1)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID UGRD "
 call ESMF_FieldScatter(u_input_grid, dummy3dflip, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)

! vgrd

 if (myrank <= max_procs-1) then
   error=nf90_inq_varid(ncid, 'vgrd', id_var)
   call netcdf_err(error, 'reading vgrd field id' )
   error=nf90_get_var(ncid, id_var, dummy3d, start=start, count=count)
   call netcdf_err(error, 'reading vgrd field' )
 endif

 call mpi_gatherv(dummy3d, iscnt, mpi_real, &
                  dummy3dall, ircnt, displ, mpi_real, &
                  0, mpi_comm_world, error)
 if (error /= 0) call error_handler("IN mpi_gatherv of vgrd", error)

 if (myrank == 0) then
   dummy3dflip(:,:,1:lev_input) = dummy3dall(:,:,lev_input:1:-1)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID VGRD "
 call ESMF_FieldScatter(v_input_grid, dummy3dflip, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)

! tracers

 do n = 1, num_tracers

   if (myrank <= max_procs-1) then
     error=nf90_inq_varid(ncid, tracers_input(n), id_var)
     call netcdf_err(error, 'reading tracer field id' )
     error=nf90_get_var(ncid, id_var, dummy3d, start=start, count=count)
     call netcdf_err(error, 'reading tracer field' )
   endif

   call mpi_gatherv(dummy3d, iscnt, mpi_real, &
                    dummy3dall, ircnt, displ, mpi_real, &
                    0, mpi_comm_world, error)
   if (error /= 0) call error_handler("IN mpi_gatherv of tracer", error)

   if (myrank == 0) then
     dummy3dflip(:,:,1:lev_input) = dummy3dall(:,:,lev_input:1:-1)
     where(dummy3dflip < 0.0) dummy3dflip = 0.0
   endif

   print*,"- CALL FieldScatter FOR INPUT GRID ", tracers_input(n)
   call ESMF_FieldScatter(tracers_input_grid(n), dummy3dflip, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)

 enddo

! dzdt   set to zero for now.

 if (myrank == 0) then
   dummy3dflip = 0.0
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID DZDT"
 call ESMF_FieldScatter(dzdt_input_grid, dummy3dflip, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 
 deallocate(dummy3dflip, dummy3dall, dummy3d)

! terrain 

 if (myrank==0) then
   print*,"- READ TERRAIN."
   error=nf90_inq_varid(ncid, 'hgtsfc', id_var)
   call netcdf_err(error, 'reading hgtsfc field id' )
   error=nf90_get_var(ncid, id_var, dummy)
   call netcdf_err(error, 'reading hgtsfc field' )
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

! surface pressure

 if (myrank==0) then
   print*,"- READ SURFACE P."
   error=nf90_inq_varid(ncid, 'pressfc', id_var)
   call netcdf_err(error, 'reading pressfc field id' )
   error=nf90_get_var(ncid, id_var, dummy)
   call netcdf_err(error, 'reading pressfc field' )
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SURFACE P."
 call ESMF_FieldScatter(ps_input_grid, dummy, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 deallocate(kcount, startk, displ, ircnt, dummy)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d cartesian winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute pressure.
!---------------------------------------------------------------------------

 print*,"- CALL FieldGet FOR PRESSURE."
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=presptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR DELTA PRESSURE."
 call ESMF_FieldGet(dpres_input_grid, &
                    farrayPtr=dpresptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(pres_interface(levp1_input))

!---------------------------------------------------------------------------
! Compute 3-d pressure.
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!  When ingesting gaussian netcdf files, the mid-layer
!  surface pressure are computed top down from delta-p
!  The surface pressure in the file is not used.  According
!  to Jun Wang, after the model's write component interpolates from the
!  cubed-sphere grid to the gaussian grid, the surface pressure is
!  no longer consistent with the delta p.
!---------------------------------------------------------------------------

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     pres_interface(levp1_input) = phalf(1) * 100.0_8
     do k = lev_input, 1, -1
       pres_interface(k) = pres_interface(k+1) + dpresptr(i,j,k)
     enddo
     psptr(i,j) = pres_interface(1)
     do k = 1, lev_input
       presptr(i,j,k) = (pres_interface(k) + pres_interface(k+1)) / 2.0_8
     enddo
   enddo
 enddo

 deallocate(pres_interface, phalf)

 call ESMF_FieldDestroy(dpres_input_grid, rc=rc)

 end subroutine read_input_atm_gaussian_netcdf_file

!---------------------------------------------------------------------------
! Read input grid fv3 atmospheric tiled history files in netcdf format.
!
! Routine reads tiled files in parallel.  Tile 1 is read by 
! localpet 0; tile 2 by localpet 1, etc.  The number of pets
! must be equal to or greater than the number of tiled files.  
!---------------------------------------------------------------------------

 subroutine read_input_atm_tiled_history_file(localpet)

 implicit none

 include 'mpif.h'

 integer, intent(in)             :: localpet

 character(len=500)              :: tilefile

 integer                         :: error, ncid, rc, tile
 integer                         :: id_dim, idim_input, jdim_input
 integer                         :: id_var, i, j, k, n
 integer                         :: clb(3), cub(3), num_tracers_file

 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile_3d(:,:,:)
 real(esmf_kind_r8), pointer     :: presptr(:,:,:), dpresptr(:,:,:)
 real(esmf_kind_r8), pointer     :: psptr(:,:)
 real(esmf_kind_r8), allocatable :: pres_interface(:), phalf(:)

 print*,"- READ INPUT ATMOS DATA FROM TILED HISTORY FILES."

 tilefile = trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(1))
 error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(tilefile) )

 error=nf90_inq_dimid(ncid, 'grid_xt', id_dim)
 call netcdf_err(error, 'reading grid_xt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=idim_input)
 call netcdf_err(error, 'reading grid_xt value' )

 error=nf90_inq_dimid(ncid, 'grid_yt', id_dim)
 call netcdf_err(error, 'reading grid_yt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=jdim_input)
 call netcdf_err(error, 'reading grid_yt value' )

 if (idim_input /= i_input .or. jdim_input /= j_input) then
   call error_handler("DIMENSION MISMATCH BETWEEN SFC AND OROG FILES.", 2)
 endif

 error=nf90_inq_dimid(ncid, 'pfull', id_dim)
 call netcdf_err(error, 'reading pfull id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=lev_input)
 call netcdf_err(error, 'reading pfull value' )

 error=nf90_inq_dimid(ncid, 'phalf', id_dim)
 call netcdf_err(error, 'reading phalf id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=levp1_input)
 call netcdf_err(error, 'reading phalf value' )
 allocate(phalf(levp1_input))
 error=nf90_inq_varid(ncid, 'phalf', id_var)
 call netcdf_err(error, 'getting phalf varid' )
 error=nf90_get_var(ncid, id_var, phalf)
 call netcdf_err(error, 'reading phalf varid' )

 error=nf90_get_att(ncid, nf90_global, 'ncnsto', num_tracers_file)
 call netcdf_err(error, 'reading ntracer value' )

 error = nf90_close(ncid)

 print*,'- FILE HAS ', num_tracers_file, ' TRACERS.'
 print*,'- WILL PROCESS ', num_tracers, ' TRACERS.'

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 print*,"- CALL FieldCreate FOR INPUT GRID DELTA PRESSURE."
 dpres_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 if (localpet < num_tiles_input_grid) then
   allocate(data_one_tile(i_input,j_input))
   allocate(data_one_tile_3d(i_input,j_input,lev_input))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_3d(0,0,0))
 endif

 if (localpet < num_tiles_input_grid) then
   tile = localpet+1
   tilefile= trim(data_dir_input_grid) // "/" // trim(atm_files_input_grid(tile))
   print*,"- READ ATMOSPHERIC DATA FROM: ", trim(tilefile)
   error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
   call netcdf_err(error, 'opening: '//trim(tilefile) )
 endif

 if (localpet < num_tiles_input_grid) then
!  print*,"- READ VERTICAL VELOCITY."
!  error=nf90_inq_varid(ncid, 'dzdt', id_var)
!  call netcdf_err(error, 'reading field id' )
!  error=nf90_get_var(ncid, id_var, data_one_tile_3d)
!  call netcdf_err(error, 'reading field' )
!  data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)

! Using w from the tiled history files has caused problems.  
! Set to zero.
   data_one_tile_3d = 0.0_8
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID VERTICAL VELOCITY."
   call ESMF_FieldScatter(dzdt_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 do n = 1, num_tracers

   if (localpet < num_tiles_input_grid) then
     print*,"- READ ", trim(tracers_input(n))
     error=nf90_inq_varid(ncid, tracers_input(n), id_var)
     call netcdf_err(error, 'reading field id' )
     error=nf90_get_var(ncid, id_var, data_one_tile_3d)
     call netcdf_err(error, 'reading field' )
     data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
   endif

   do tile = 1, num_tiles_input_grid
     print*,"- CALL FieldScatter FOR INPUT GRID TRACER ", trim(tracers_input(n))
     call ESMF_FieldScatter(tracers_input_grid(n), data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)
   enddo

 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ TEMPERATURE."
   error=nf90_inq_varid(ncid, 'tmp', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
   call ESMF_FieldScatter(temp_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ U-WIND."
   error=nf90_inq_varid(ncid, 'ugrd', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID U."
   call ESMF_FieldScatter(u_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ V-WIND."
   error=nf90_inq_varid(ncid, 'vgrd', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID V."
   call ESMF_FieldScatter(v_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ SURFACE PRESSURE."
   error=nf90_inq_varid(ncid, 'pressfc', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile)
   call netcdf_err(error, 'reading field' )
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID SURFACE PRESSURE."
   call ESMF_FieldScatter(ps_input_grid, data_one_tile, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ TERRAIN."
   error=nf90_inq_varid(ncid, 'hgtsfc', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile)
   call netcdf_err(error, 'reading field' )
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN."
   call ESMF_FieldScatter(terrain_input_grid, data_one_tile, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) then
   print*,"- READ DELTA PRESSURE."
   error=nf90_inq_varid(ncid, 'dpres', id_var)
   call netcdf_err(error, 'reading field id' )
   error=nf90_get_var(ncid, id_var, data_one_tile_3d)
   call netcdf_err(error, 'reading field' )
   data_one_tile_3d(:,:,1:lev_input) = data_one_tile_3d(:,:,lev_input:1:-1)
 endif

 do tile = 1, num_tiles_input_grid
   print*,"- CALL FieldScatter FOR INPUT DELTA PRESSURE."
   call ESMF_FieldScatter(dpres_input_grid, data_one_tile_3d, rootpet=tile-1, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)
 enddo

 if (localpet < num_tiles_input_grid) error = nf90_close(ncid)

 deallocate(data_one_tile_3d, data_one_tile)

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d cartesian winds.
!---------------------------------------------------------------------------

 call convert_winds

!---------------------------------------------------------------------------
! Compute pressure.
!---------------------------------------------------------------------------

 print*,"- CALL FieldGet FOR PRESSURE."
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=presptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR DELTA PRESSURE."
 call ESMF_FieldGet(dpres_input_grid, &
                    farrayPtr=dpresptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 call ESMF_FieldGet(ps_input_grid, &
                    farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 allocate(pres_interface(levp1_input))

!---------------------------------------------------------------------------
! Compute 3-d pressure.
!---------------------------------------------------------------------------

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     pres_interface(1) = psptr(i,j)
     do k = 2, levp1_input
       pres_interface(k) = pres_interface(k-1) - dpresptr(i,j,k-1)
     enddo
     do k = 1, lev_input
       presptr(i,j,k) = (pres_interface(k) + pres_interface(k+1)) / 2.0_8
     enddo
   enddo
 enddo

 deallocate(pres_interface, phalf)

 call ESMF_FieldDestroy(dpres_input_grid, rc=rc)

 end subroutine read_input_atm_tiled_history_file
 
!---------------------------------------------------------------------------
! Read input grid atmospheric fv3gfs grib2 files.
!---------------------------------------------------------------------------

 subroutine read_input_atm_grib2_file(localpet)

 use wgrib2api
 
 use grib2_util, only                   : rh2spfh, convert_omega

 implicit none

 integer, intent(in)                   :: localpet
 
 integer, parameter                    :: ntrac_max=14

 character(len=300)                    :: the_file
 character(len=20)                     :: vlevtyp, vname, lvl_str,lvl_str_space, &
                                          trac_names_grib_1(ntrac_max), &
                                          trac_names_grib_2(ntrac_max), &
                                          trac_names_vmap(ntrac_max), &
                                          tracers_input_grib_1(num_tracers), &
                                          tracers_input_grib_2(num_tracers), &
                                          tmpstr, & 
                                          method, tracers_input_vmap(num_tracers), &
                                          tracers_default(ntrac_max), vname2
 character (len=500)                   :: metadata

 integer                               :: i, j, k, n, lvl_str_space_len
 integer                               :: rc, clb(3), cub(3)
 integer                               :: vlev, iret,varnum

 integer                               :: len_str
 logical                               :: lret

 logical                               :: conv_omega=.false., &
                                          hasspfh=.true.

 real(esmf_kind_r8), allocatable       :: rlevs(:)
 real(esmf_kind_r4), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:), dummy2d_8(:,:)
 real(esmf_kind_r8), pointer           :: presptr(:,:,:), psptr(:,:),tptr(:,:,:), &
                                          qptr(:,:,:), wptr(:,:,:),  &
                                          uptr(:,:,:), vptr(:,:,:)
 real(esmf_kind_r4)                    :: value
 real(esmf_kind_r8), parameter         :: p0 = 100000.0
 
 
 tracers(:) = "NULL"
 !trac_names_grib = (/":SPFH:",":CLWR:", "O3MR",":CICE:", ":RWMR:",":SNMR:",":GRLE:", &
 !              ":TCDC:", ":NCCICE:",":SPNCR:", ":NCONCD:",":PMTF:",":PMTC:",":TKE:"/)
 trac_names_grib_1 = (/":var0_2", ":var0_2",  ":var0_2",  ":var0_2",  ":var0_2",":var0_2", \
                       ":var0_2", ":var0_2", ":var0_2", ":var0_2", ":var0_2",":var0_2", \
                       ":var0_2", ":var0_2"/)
 trac_names_grib_2 = (/"_1_0:   ", "_1_22:  ",  "_14_192:", "_1_23:  ", "_1_24:  ","_1_25:  ", \
                       "_1_32:  ", "_6_1:   ",  "_6_29:  ", "_1_100: ", "_6_28:  ","_13_193:", \
                       "_13_192:", "_2_2:   "/)
 trac_names_vmap = (/"sphum   ", "liq_wat ", "o3mr    ", "ice_wat ", &
                     "rainwat ", "snowwat ", "graupel ", "cld_amt ", "ice_nc  ", &
                     "rain_nc ", "water_nc", "liq_aero", "ice_aero", &
                     "sgs_tke "/)
 tracers_default = (/"sphum   ", "liq_wat ", "o3mr    ", "ice_wat ", &
                     "rainwat ", "snowwat ", "graupel ", "cld_amt ", "ice_nc  ", &
                     "rain_nc ", "water_nc", "liq_aero", "ice_aero", &
                     "sgs_tke "/)

 the_file = trim(data_dir_input_grid) // "/" // trim(grib2_file_input_grid)

 print*,"- READ ATMOS DATA FROM GRIB2 FILE: ", trim(the_file)
 print*,"- USE INVENTORY FILE ", inv_file

 print*,"- OPEN FILE."
 inquire(file=the_file,exist=lret)
 if (.not.lret) call error_handler("OPENING GRIB2 ATM FILE.", iret)

 print*,"- READ VERTICAL COORDINATE."
 iret = grb2_inq(the_file,inv_file,":var_0_2","_0_0:"," hybrid level:")
  
 if (iret <= 0) then
   lvl_str = "mb:" 
   lvl_str_space = " mb:"
   lvl_str_space_len = 4
   iret = grb2_inq(the_file,inv_file,":UGRD:",lvl_str_space)
   lev_input=iret
   if (localpet == 0) print*,"- DATA IS ON ", lev_input, " ISOBARIC LEVELS."
 else
   call error_handler("HYBRID VERTICAL COORD DATA NOT SUPPORTED", -1)
 endif

 allocate(slevs(lev_input))
 allocate(rlevs(lev_input))
 levp1_input = lev_input + 1
    
! Get the vertical levels, and search string by sequential reads

 do i = 1,lev_input
   iret=grb2_inq(the_file,inv_file,':UGRD:',trim(lvl_str),sequential=i-1,desc=metadata)
   if (iret.ne.1) call error_handler(" IN SEQUENTIAL FILE READ.", iret)
    
   j = index(metadata,':UGRD:') + len(':UGRD:')
   k = index(metadata,trim(lvl_str_space)) + len(trim(lvl_str_space))-1

   read(metadata(j:k),*) rlevs(i)

   slevs(i) = metadata(j-1:k)	
   rlevs(i) = rlevs(i) * 100.0
   if (localpet==0) print*, "- LEVEL = ", slevs(i)
 enddo

! Jili Dong add sort to re-order isobaric levels.

 call quicksort(rlevs,1,lev_input)

 do i = 1,lev_input
   write(slevs(i),"(F20.10)") rlevs(i)/100.0
   len_str = len_trim(slevs(i))

   do while (slevs(i)(len_str:len_str) .eq. '0')
        slevs(i) = slevs(i)(:len_str-1)
        len_str = len_str - 1
   end do

   if (slevs(i)(len_str:len_str) .eq. '.') then
     slevs(i) = slevs(i)(:len_str-1)
     len_str = len_str - 1
   end if

   slevs(i) = trim(slevs(i))

   slevs(i) = ":"//trim(adjustl(slevs(i)))//" mb:"
   if (localpet==0) print*, "- LEVEL AFTER SORT = ",slevs(i)
 enddo

 if (localpet == 0) print*,"- FIND SPFH OR RH IN FILE"
 iret = grb2_inq(the_file,inv_file,trac_names_grib_1(1),trac_names_grib_2(1),lvl_str_space)

 if (iret <= 0) then
   iret = grb2_inq(the_file,inv_file, ':var0_2','_1_1:',lvl_str_space)
   if (iret <= 0) call error_handler("READING ATMOSPHERIC WATER VAPOR VARIABLE.", iret)
   hasspfh = .false.
   trac_names_grib_2(1)='_1_1:'
   if (localpet == 0) print*,"- FILE CONTAINS RH."
 else
   if (localpet == 0) print*,"- FILE CONTAINS SPFH."
 endif
   
 print*,"- COUNT NUMBER OF TRACERS TO BE READ IN BASED ON PHYSICS SUITE TABLE"
 do n = 1, num_tracers

   vname = tracers_input(n)

   i = maxloc(merge(1.,0.,trac_names_vmap == vname),dim=1)

   tracers_input_grib_1(n) = trac_names_grib_1(i)
   tracers_input_grib_2(n) = trac_names_grib_2(i)
   tracers_input_vmap(n)=trac_names_vmap(i)
   tracers(n)=tracers_default(i)

 enddo

 if (localpet==0) print*, "- NUMBER OF TRACERS IN FILE = ", num_tracers

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 if (localpet == 0) then
   allocate(dummy2d(i_input,j_input))
   allocate(dummy2d_8(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lev_input))
 else
   allocate(dummy2d(0,0))
   allocate(dummy2d_8(0,0))
   allocate(dummy3d(0,0,0))
 endif

!-----------------------------------------------------------------------
! Fields in non-native files read in from top to bottom. We will
! flip indices later.  This program expects bottom to top.
!-----------------------------------------------------------------------
 
 if (localpet == 0) then
   print*,"- READ TEMPERATURE."
   vname = ":TMP:"   
    do vlev = 1, lev_input
      iret = grb2_inq(the_file,inv_file,vname,slevs(vlev),data2=dummy2d)
      if (iret<=0) then 
        call error_handler("READING IN TEMPERATURE AT LEVEL "//trim(slevs(vlev)),iret)
      endif
      dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)
      print*,'temp check after read ',vlev, dummy3d(1,1,vlev)
    enddo
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
 call ESMF_FieldScatter(temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 do n = 1, num_tracers

   if (localpet == 0) print*,"- READ ", trim(tracers_input_vmap(n))
   vname = tracers_input_vmap(n)
   call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
                       this_field_var_name=tmpstr,loc=varnum)
   if (n==1 .and. .not. hasspfh) then 
        print*,"- CALL FieldGather TEMPERATURE." 
        call ESMF_FieldGather(temp_input_grid,dummy3d,rootPet=0, tile=1, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGet", rc) 
   endif
   if (localpet == 0) then
     vname = trim(tracers_input_grib_1(n))
     vname2 = trim(tracers_input_grib_2(n))
     
     do vlev = 1, lev_input
      iret = grb2_inq(the_file,inv_file,vname,slevs(vlev),vname2,data2=dummy2d)
     
      if (iret <= 0) then
        call handle_grib_error(vname, slevs(vlev),method,value,varnum,iret,var=dummy2d)
        if (iret==1) then ! missing_var_method == skip or no entry
          if (trim(vname2)=="_1_0:" .or. trim(vname2) == "_1_1:" .or.  &
              trim(vname2) == ":14:192:") then
            call error_handler("READING IN "//trim(vname)//" AT LEVEL "//trim(slevs(vlev))&
                      //". SET A FILL VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
          else
            exit
          endif
        endif
      endif
      
      if (n==1 .and. .not. hasspfh) then 
        call rh2spfh(dummy2d,rlevs(vlev),dummy3d(:,:,vlev))
      endif

       print*,'tracer ',vlev, maxval(dummy2d),minval(dummy2d)
       dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)
     enddo
   endif

   if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input_vmap(n))
   call ESMF_FieldScatter(tracers_input_grid(n), dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo
 
 if (localpet==0) then
   do vlev = 1, lev_input
 
     vname = ":var0_2"
     vname2 = "_2_2:"
     iret = grb2_inq(the_file,inv_file,vname,vname2,slevs(vlev),data2=dummy2d)
     if (iret<=0) then 
       call error_handler("READING UWIND AT LEVEL "//trim(slevs(vlev)),iret)
     endif

     print*, 'max, min U ', minval(dummy2d), maxval(dummy2d)
     dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)

   enddo
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT U-WIND."
 call ESMF_FieldScatter(u_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet==0) then
   do vlev = 1, lev_input
 
     vname = ":var0_2"
     vname2 = "_2_3:"
     iret = grb2_inq(the_file,inv_file,vname,vname2,slevs(vlev),data2=dummy2d)
     if (iret<=0) then 
       call error_handler("READING VWIND AT LEVEL "//trim(slevs(vlev)),iret)
     endif

     print*, 'max, min V ', minval(dummy2d), maxval(dummy2d)
     dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)

   enddo
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT V-WIND."
 call ESMF_FieldScatter(v_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SURFACE PRESSURE."
   !vname = ":PRES:"
   vname = ":var0_2"
   vname2 = "_3_0:"
   vlevtyp = ":surface:"
   iret = grb2_inq(the_file,inv_file,vname,vname2,vlevtyp,data2=dummy2d)
   if (iret <= 0) call error_handler("READING SURFACE PRESSURE RECORD.", iret)
   dummy2d_8 = real(dummy2d,esmf_kind_r8)
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID SURFACE PRESSURE."
 call ESMF_FieldScatter(ps_input_grid, dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ DZDT."
   vname = "dzdt"
   call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
                         loc=varnum)
   !vname = ":DZDT:"
   vname = ":var0_2"
   vname2 = "_2_9:"
   do vlev = 1, lev_input
     iret = grb2_inq(the_file,inv_file,vname,vname2,slevs(vlev),data2=dummy2d)
     if (iret <= 0 ) then
       print*,"DZDT not available at level ", trim(slevs(vlev)), " so checking for VVEL"
       !vname = ":VVEL:"
       vname2 = "_2_8:"
       iret = grb2_inq(the_file,inv_file,vname,vname2,slevs(vlev),data2=dummy2d)
       
       
       if (iret <= 0) then
        call handle_grib_error(vname, slevs(vlev),method,value,varnum,iret,var=dummy2d)
        if (iret==1) then ! missing_var_method == skip 
          cycle
        endif
       else
        conv_omega = .true.
       endif
       
     endif
     print*,'dzdt ',vlev, maxval(dummy2d),minval(dummy2d)
     dummy3d(:,:,vlev) = dummy2d
   enddo
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT DZDT."
 call ESMF_FieldScatter(dzdt_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ TERRAIN."
   !vname = ":HGT:"
   vname = ":var0_2"
    vname2 = "_3_5:"
   vlevtyp = ":surface:"
   iret = grb2_inq(the_file,inv_file,vname,vname2,vlevtyp,data2=dummy2d)
   if (iret <= 0) call error_handler("READING TERRAIN HEIGHT RECORD.", iret)
   dummy2d_8 = real(dummy2d,esmf_kind_r8)
 endif

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d, dummy3d, dummy2d_8)
 
!---------------------------------------------------------------------------
! Flip 'z' indices to all 3-d variables.  Data is read in from model
! top to surface.  This program expects surface to model top.
!---------------------------------------------------------------------------
    
 if (localpet == 0) print*,"- CALL FieldGet FOR SURFACE PRESSURE."
 nullify(psptr)
 call ESMF_FieldGet(ps_input_grid, &
                      farrayPtr=psptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)
      
 nullify(presptr)
 if (localpet == 0) print*,"- CALL FieldGet FOR 3-D PRESSURE."
 call ESMF_FieldGet(pres_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=presptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 nullify(tptr)
 if (localpet == 0) print*,"- CALL FieldGet TEMPERATURE."  
 call ESMF_FieldGet(temp_input_grid, &
                    farrayPtr=tptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc) 
  
 nullify(uptr)
 if (localpet == 0) print*,"- CALL FieldGet FOR U"
 call ESMF_FieldGet(u_input_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
 nullify(vptr)
 if (localpet == 0) print*,"- CALL FieldGet FOR V"
 call ESMF_FieldGet(v_input_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
  
 nullify(wptr)
 if (localpet == 0) print*,"- CALL FieldGet FOR W"
 call ESMF_FieldGet(dzdt_input_grid, &
                    farrayPtr=wptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
 
  if (localpet == 0) print*,"- CALL FieldGet FOR TRACERS."
  do n=1,num_tracers
    nullify(qptr)
    call ESMF_FieldGet(tracers_input_grid(n), &
                    farrayPtr=qptr, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)
    do i = clb(1),cub(1)
      do j = clb(2),cub(2)
        qptr(i,j,:) = qptr(i,j,lev_input:1:-1)
      end do
    end do
  end do
  
  do i = clb(1),cub(1)
    do j = clb(2),cub(2)
      presptr(i,j,:) = rlevs(lev_input:1:-1)
      tptr(i,j,:) = tptr(i,j,lev_input:1:-1)
      uptr(i,j,:) = uptr(i,j,lev_input:1:-1)
      vptr(i,j,:) = vptr(i,j,lev_input:1:-1)
      wptr(i,j,:) = wptr(i,j,lev_input:1:-1)
    end do
  end do

 if (localpet == 0) then
   print*,'psfc is ',clb(1),clb(2),psptr(clb(1),clb(2))
   print*,'pres is ',cub(1),cub(2),presptr(cub(1),cub(2),:) 
     
   print*,'pres check 1',localpet,maxval(presptr(clb(1):cub(1),clb(2):cub(2),1)), &
              minval(presptr(clb(1):cub(1),clb(2):cub(2),1))
   print*,'pres check lev',localpet,maxval(presptr(clb(1):cub(1),clb(2):cub(2), &
            lev_input)),minval(presptr(clb(1):cub(1),clb(2):cub(2),lev_input))
 endif
 
!---------------------------------------------------------------------------
! Convert from 2-d to 3-d component winds.
!---------------------------------------------------------------------------

 call convert_winds
 
!---------------------------------------------------------------------------
! Convert dpdt to dzdt if needed
!---------------------------------------------------------------------------

 if (conv_omega) then

  if (localpet == 0)  print*,"- CONVERT FROM OMEGA TO DZDT."

  nullify(tptr)
  if (localpet == 0) print*,"- CALL FieldGet TEMPERATURE."  
  call ESMF_FieldGet(temp_input_grid, &
                    farrayPtr=tptr, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc) 
    
  nullify(qptr)
  if (localpet == 0) print*,"- CALL FieldGet SPECIFIC HUMIDITY."  
  call ESMF_FieldGet(tracers_input_grid(1), &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=qptr, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
  nullify(wptr)
  if (localpet == 0) print*,"- CALL FieldGet DZDT." 
  call ESMF_FieldGet(dzdt_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=wptr, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
  
  nullify(presptr)
  call ESMF_FieldGet(pres_input_grid, &
                    farrayPtr=presptr, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)
    
  call convert_omega(wptr,presptr,tptr,qptr,clb,cub)
  
 endif
 
 end subroutine read_input_atm_grib2_file

!---------------------------------------------------------------------------
! Read input grid surface data from a spectral gfs gaussian sfcio file.
! Prior to July 19, 2017.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_gfs_sfcio_file(localpet)
 
 use sfcio_module

 implicit none

 integer, intent(in)                   :: localpet

 character(len=300)                    :: the_file

 integer(sfcio_intkind)                :: iret
 integer                               :: rc

 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)

 type(sfcio_head)                      :: sfchead
 type(sfcio_dbta)                      :: sfcdata

 the_file = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))

 print*,"- READ SURFACE DATA IN SFCIO FORMAT."
 print*,"- OPEN AND READ: ",trim(the_file)
 call sfcio_sropen(23, trim(the_file), iret)
 if (iret /= 0) then
   rc=iret
   call error_handler("OPENING FILE", rc)
 endif

 call sfcio_srhead(23, sfchead, iret)
 if (iret /= 0) then
   rc=iret
   call error_handler("READING HEADER", rc)
 endif

 if (localpet == 0) then
   call sfcio_aldbta(sfchead, sfcdata, iret)
   if (iret /= 0) then
     rc=iret
     call error_handler("ALLOCATING DATA.", rc)
   endif
   call sfcio_srdbta(23, sfchead, sfcdata, iret)
   if (iret /= 0) then
     rc=iret
     call error_handler("READING DATA.", rc)
   endif
   allocate(dummy2d(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lsoil_input))
 else
   allocate(dummy2d(0,0))
   allocate(dummy3d(0,0,0))
 endif

 if (localpet == 0) dummy2d = sfcdata%slmsk

 print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
 call ESMF_FieldScatter(landsea_mask_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%zorl

 print*,"- CALL FieldScatter FOR INPUT Z0."
 call ESMF_FieldScatter(z0_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = nint(sfcdata%vtype)

 print*,"- CALL FieldScatter FOR INPUT VEG TYPE."
 call ESMF_FieldScatter(veg_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

! Prior to July, 2017, gfs used zobler soil types.  '13' indicates permanent land ice.
 veg_type_landice_input = 13

 if (localpet == 0) dummy2d = sfcdata%canopy

 print*,"- CALL FieldScatter FOR INPUT CANOPY MC."
 call ESMF_FieldScatter(canopy_mc_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%fice

 print*,"- CALL FieldScatter FOR INPUT ICE FRACTION."
 call ESMF_FieldScatter(seaice_fract_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%hice

 print*,"- CALL FieldScatter FOR INPUT ICE DEPTH."
 call ESMF_FieldScatter(seaice_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%tisfc

 print*,"- CALL FieldScatter FOR INPUT ICE SKIN TEMP."
 call ESMF_FieldScatter(seaice_skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%snwdph ! mm (expected by program)

 print*,"- CALL FieldScatter FOR INPUT SNOW DEPTH."
 call ESMF_FieldScatter(snow_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%sheleg

 print*,"- CALL FieldScatter FOR INPUT SNOW LIQUID EQUIV."
 call ESMF_FieldScatter(snow_liq_equiv_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%t2m

 print*,"- CALL FieldScatter FOR INPUT T2M."
 call ESMF_FieldScatter(t2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%q2m

 print*,"- CALL FieldScatter FOR INPUT Q2M."
 call ESMF_FieldScatter(q2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%tprcp

 print*,"- CALL FieldScatter FOR INPUT TPRCP."
 call ESMF_FieldScatter(tprcp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%f10m

 print*,"- CALL FieldScatter FOR INPUT F10M."
 call ESMF_FieldScatter(f10m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%uustar

 print*,"- CALL FieldScatter FOR INPUT USTAR."
 call ESMF_FieldScatter(ustar_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%ffmm

 print*,"- CALL FieldScatter FOR INPUT FFMM."
 call ESMF_FieldScatter(ffmm_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%srflag

 print*,"- CALL FieldScatter FOR INPUT SRFLAG."
 call ESMF_FieldScatter(srflag_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%tsea

 print*,"- CALL FieldScatter FOR INPUT SKIN TEMP."
 call ESMF_FieldScatter(skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = nint(sfcdata%stype)

 print*,"- CALL FieldScatter FOR INPUT SOIL TYPE."
 call ESMF_FieldScatter(soil_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = sfcdata%orog

 print*,"- CALL FieldScatter FOR INPUT TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy3d = sfcdata%slc

 print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_liq_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy3d = sfcdata%smc

 print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_tot_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy3d = sfcdata%stc

 print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
 call ESMF_FieldScatter(soil_temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d, dummy3d)
 call sfcio_axdbta(sfcdata, iret)

 call sfcio_sclose(23, iret)

 end subroutine read_input_sfc_gfs_sfcio_file

!---------------------------------------------------------------------------
! Read input grid surface data from a spectral gfs gaussian nemsio file.
! Format used by gfs starting July 19, 2017.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_gfs_gaussian_nemsio_file(localpet)
 
 implicit none

 integer, intent(in)                   :: localpet

 character(len=300)                    :: the_file

 integer                               :: rc

 real(nemsio_realkind), allocatable    :: dummy(:)
 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)

 type(nemsio_gfile)                    :: gfile

 the_file = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))

 if (localpet == 0) then
   allocate(dummy3d(i_input,j_input,lsoil_input))
   allocate(dummy2d(i_input,j_input))
   allocate(dummy(i_input*j_input))
   print*,"- OPEN FILE ", trim(the_file)
   call nemsio_open(gfile, the_file, "read", iret=rc)
   if (rc /= 0) call error_handler("OPENING FILE.", rc)
 else
   allocate(dummy3d(0,0,0))
   allocate(dummy2d(0,0))
   allocate(dummy(0))
 endif

 if (localpet == 0) then
   print*,"- READ TERRAIN."
   call nemsio_readrecv(gfile, "orog", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING TERRAIN.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'orog ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ LANDSEA MASK."
   call nemsio_readrecv(gfile, "land", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LANDSEA MASK.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'landmask ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
 call ESMF_FieldScatter(landsea_mask_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ SEAICE FRACTION."
   call nemsio_readrecv(gfile, "icec", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE FRACTION.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'icec ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
 call ESMF_FieldScatter(seaice_fract_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SEAICE DEPTH."
   call nemsio_readrecv(gfile, "icetk", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE DEPTH.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'icetk ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE DEPTH."
 call ESMF_FieldScatter(seaice_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SEAICE SKIN TEMPERATURE."
   call nemsio_readrecv(gfile, "tisfc", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE SKIN TEMP.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'ti ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
 call ESMF_FieldScatter(seaice_skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SNOW LIQUID EQUIVALENT."
   call nemsio_readrecv(gfile, "weasd", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SNOW LIQUID EQUIVALENT.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'weasd ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
 call ESMF_FieldScatter(snow_liq_equiv_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SNOW DEPTH."
   call nemsio_readrecv(gfile, "snod", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SNOW DEPTH.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'snod ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
 call ESMF_FieldScatter(snow_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ VEG TYPE."
   call nemsio_readrecv(gfile, "vtype", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING VEG TYPE", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'vtype ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID VEG TYPE."
 call ESMF_FieldScatter(veg_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SOIL TYPE."
   call nemsio_readrecv(gfile, "sotyp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SOIL TYPE.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'sotype ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
 call ESMF_FieldScatter(soil_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ T2M."
   call nemsio_readrecv(gfile, "tmp", "2 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING T2M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'t2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID T2M."
 call ESMF_FieldScatter(t2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Q2M."
   call nemsio_readrecv(gfile, "spfh", "2 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING Q2M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'q2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
 call ESMF_FieldScatter(q2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ TPRCP."
   call nemsio_readrecv(gfile, "tprcp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING TPRCP.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'tprcp ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TPRCP."
 call ESMF_FieldScatter(tprcp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ FFMM."
   call nemsio_readrecv(gfile, "ffmm", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING FFMM.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'ffmm ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID FFMM"
 call ESMF_FieldScatter(ffmm_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ USTAR."
   call nemsio_readrecv(gfile, "fricv", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING USTAR.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'fricv ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID USTAR"
 call ESMF_FieldScatter(ustar_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = 0.0
 print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
 call ESMF_FieldScatter(srflag_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SKIN TEMPERATURE."
   call nemsio_readrecv(gfile, "tmp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SKIN TEMPERATURE.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'tmp ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
 call ESMF_FieldScatter(skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ F10M."
   call nemsio_readrecv(gfile, "f10m", "10 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING F10M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'f10m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID F10M."
 call ESMF_FieldScatter(f10m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ CANOPY MOISTURE CONTENT."
   call nemsio_readrecv(gfile, "cnwat", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING CANOPY MOISTURE CONTENT.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'cnwat ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID CANOPY MOISTURE CONTENT."
 call ESMF_FieldScatter(canopy_mc_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Z0."
   call nemsio_readrecv(gfile, "sfcr", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING Z0.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'sfcr ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Z0."
 call ESMF_FieldScatter(z0_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d)

 if (localpet == 0) then
   print*,"- READ LIQUID SOIL MOISTURE."
   call nemsio_readrecv(gfile, "slc", "soil layer", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "slc", "soil layer", 2, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "slc", "soil layer", 3, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "slc", "soil layer", 4, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'slc ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_liq_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ TOTAL SOIL MOISTURE."
   call nemsio_readrecv(gfile, "smc", "soil layer", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "smc", "soil layer", 2, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "smc", "soil layer", 3, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "smc", "soil layer", 4, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'smc ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_tot_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SOIL TEMPERATURE."
   call nemsio_readrecv(gfile, "stc", "soil layer", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 SOIL TEMP.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "stc", "soil layer", 2, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 SOIL TEMP.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "stc", "soil layer", 3, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 SOIL TEMP.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "stc", "soil layer", 4, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 SOIL TEMP.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'stc ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
 call ESMF_FieldScatter(soil_temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy3d, dummy)

 if (localpet == 0) call nemsio_close(gfile)

 end subroutine read_input_sfc_gfs_gaussian_nemsio_file

!---------------------------------------------------------------------------
! Read input grid surface data from an fv3 gaussian nemsio file.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_gaussian_nemsio_file(localpet)

 implicit none

 integer, intent(in)                   :: localpet

 character(len=250)                    :: the_file

 integer                               :: rc

 real(nemsio_realkind), allocatable    :: dummy(:)
 real(esmf_kind_r8), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)

 type(nemsio_gfile)                    :: gfile

 the_file = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))

 if (localpet == 0) then
   allocate(dummy3d(i_input,j_input,lsoil_input))
   allocate(dummy2d(i_input,j_input))
   allocate(dummy(i_input*j_input))
   print*,"- OPEN FILE ", trim(the_file)
   call nemsio_open(gfile, the_file, "read", iret=rc)
   if (rc /= 0) call error_handler("OPENING FILE.", rc)
 else
   allocate(dummy3d(0,0,0))
   allocate(dummy2d(0,0))
   allocate(dummy(0))
 endif

 if (localpet == 0) then
   print*,"- READ TERRAIN."
   call nemsio_readrecv(gfile, "orog", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING TERRAIN.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'orog ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ LANDSEA MASK."
   call nemsio_readrecv(gfile, "land", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LANDSEA MASK.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'landmask ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
 call ESMF_FieldScatter(landsea_mask_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ SEAICE FRACTION."
   call nemsio_readrecv(gfile, "icec", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE FRACTION.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'icec ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
 call ESMF_FieldScatter(seaice_fract_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SEAICE DEPTH."
   call nemsio_readrecv(gfile, "icetk", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE DEPTH.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'icetk ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE DEPTH."
 call ESMF_FieldScatter(seaice_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SEAICE SKIN TEMPERATURE."
   call nemsio_readrecv(gfile, "ti", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SEAICE SKIN TEMP.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'ti ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
 call ESMF_FieldScatter(seaice_skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SNOW LIQUID EQUIVALENT."
   call nemsio_readrecv(gfile, "weasd", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SNOW LIQUID EQUIVALENT.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'weasd ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
 call ESMF_FieldScatter(snow_liq_equiv_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SNOW DEPTH."
   call nemsio_readrecv(gfile, "snod", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SNOW DEPTH.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/)) * 1000.0_8
   print*,'snod ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
 call ESMF_FieldScatter(snow_depth_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ VEG TYPE."
   call nemsio_readrecv(gfile, "vtype", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING VEG TYPE", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'vtype ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID VEG TYPE."
 call ESMF_FieldScatter(veg_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SOIL TYPE."
   call nemsio_readrecv(gfile, "sotyp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SOIL TYPE.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'sotype ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
 call ESMF_FieldScatter(soil_type_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ T2M."
   call nemsio_readrecv(gfile, "tmp", "2 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING T2M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'t2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID T2M."
 call ESMF_FieldScatter(t2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Q2M."
   call nemsio_readrecv(gfile, "spfh", "2 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING Q2M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'q2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
 call ESMF_FieldScatter(q2m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ TPRCP."
   call nemsio_readrecv(gfile, "tprcp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING TPRCP.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'tprcp ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TPRCP."
 call ESMF_FieldScatter(tprcp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ FFMM."
   call nemsio_readrecv(gfile, "ffmm", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING FFMM.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'ffmm ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID FFMM"
 call ESMF_FieldScatter(ffmm_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ USTAR."
   call nemsio_readrecv(gfile, "fricv", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING USTAR.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'fricv ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID USTAR"
 call ESMF_FieldScatter(ustar_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) dummy2d = 0.0
 print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
 call ESMF_FieldScatter(srflag_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SKIN TEMPERATURE."
   call nemsio_readrecv(gfile, "tmp", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING SKIN TEMPERATURE.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'tmp ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
 call ESMF_FieldScatter(skin_temp_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ F10M."
   call nemsio_readrecv(gfile, "f10m", "10 m above gnd", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING F10M.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'f10m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID F10M."
 call ESMF_FieldScatter(f10m_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ CANOPY MOISTURE CONTENT."
   call nemsio_readrecv(gfile, "cnwat", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING CANOPY MOISTURE CONTENT.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'cnwat ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID CANOPY MOISTURE CONTENT."
 call ESMF_FieldScatter(canopy_mc_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Z0."
   call nemsio_readrecv(gfile, "sfcr", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING Z0.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/)) * 100.0_8 ! convert to cm
   print*,'sfcr ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Z0."
 call ESMF_FieldScatter(z0_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d)

 if (localpet == 0) then
   print*,"- READ LIQUID SOIL MOISTURE."
   call nemsio_readrecv(gfile, "soill", "0-10 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soill", "10-40 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soill", "40-100 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soill", "100-200 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 LIQUID SOIL MOIST.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'soill ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_liq_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ TOTAL SOIL MOISTURE."
   call nemsio_readrecv(gfile, "soilw", "0-10 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soilw", "10-40 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soilw", "40-100 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "soilw", "100-200 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 TOTAL SOIL MOIST.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'soilm ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_tot_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SOIL TEMPERATURE."
   call nemsio_readrecv(gfile, "tmp", "0-10 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 1 SOIL TEMP.", rc)
   dummy3d(:,:,1) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "tmp", "10-40 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 2 SOIL TEMP.", rc)
   dummy3d(:,:,2) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "tmp", "40-100 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 3 SOIL TEMP.", rc)
   dummy3d(:,:,3) = reshape(dummy, (/i_input,j_input/))
   call nemsio_readrecv(gfile, "tmp", "100-200 cm down", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING LAYER 4 SOIL TEMP.", rc)
   dummy3d(:,:,4) = reshape(dummy, (/i_input,j_input/))
   print*,'soilt ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
 call ESMF_FieldScatter(soil_temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy3d, dummy)

 if (localpet == 0) call nemsio_close(gfile)

 end subroutine read_input_sfc_gaussian_nemsio_file

!---------------------------------------------------------------------------
! Read input grid surface data tiled warm 'restart' files.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_restart_file(localpet)

 implicit none

 integer, intent(in)             :: localpet

 character(len=500)              :: tilefile

 integer                         :: error, rc
 integer                         :: id_dim, idim_input, jdim_input
 integer                         :: ncid, tile, id_var

 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile_3d(:,:,:)

!---------------------------------------------------------------------------
! Get i/j dimensions and number of soil layers from first surface file.
! Do dimensions match those from the orography file?
!---------------------------------------------------------------------------

 tilefile = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))
 print*,"- READ GRID DIMENSIONS FROM: ", trim(tilefile)
 error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(tilefile) )

 error=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
 call netcdf_err(error, 'reading xaxis_1 id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=idim_input)
 call netcdf_err(error, 'reading xaxis_1 value' )

 error=nf90_inq_dimid(ncid, 'yaxis_1', id_dim)
 call netcdf_err(error, 'reading yaxis_1 id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=jdim_input)
 call netcdf_err(error, 'reading yaxis_1 value' )

 if (idim_input /= i_input .or. jdim_input /= j_input) then
   call error_handler("DIMENSION MISMATCH BETWEEN SFC AND OROG FILES.", 1)
 endif

 error = nf90_close(ncid)

 if (localpet == 0) then
   allocate(data_one_tile(idim_input,jdim_input))
   allocate(data_one_tile_3d(idim_input,jdim_input,lsoil_input))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_3d(0,0,0))
 endif

 TERRAIN_LOOP: do tile = 1, num_tiles_input_grid

   if (localpet == 0) then
     tilefile = trim(orog_dir_input_grid) // trim(orog_files_input_grid(tile))
     print*,'- OPEN OROGRAPHY FILE: ', trim(tilefile)
     error=nf90_open(tilefile,nf90_nowrite,ncid)
     call netcdf_err(error, 'OPENING OROGRAPHY FILE' )
     error=nf90_inq_varid(ncid, 'orog_raw', id_var)
     call netcdf_err(error, 'READING OROG RECORD ID' )
     error=nf90_get_var(ncid, id_var, data_one_tile)
     call netcdf_err(error, 'READING OROG RECORD' )
     print*,'terrain check ',tile, maxval(data_one_tile)
     error=nf90_close(ncid)
   endif

   print*,"- CALL FieldScatter FOR INPUT TERRAIN."
   call ESMF_FieldScatter(terrain_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 enddo TERRAIN_LOOP

 TILE_LOOP : do tile = 1, num_tiles_input_grid

! liquid soil moisture

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('slc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata_3d=data_one_tile_3d)
  endif

  print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
  call ESMF_FieldScatter(soilm_liq_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('smc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata_3d=data_one_tile_3d)
  endif

  print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
  call ESMF_FieldScatter(soilm_tot_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('stc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata_3d=data_one_tile_3d)
  endif

  print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
  call ESMF_FieldScatter(soil_temp_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! land mask

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('slmsk', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
  call ESMF_FieldScatter(landsea_mask_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice fraction

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('fice', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
  call ESMF_FieldScatter(seaice_fract_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('hice', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE DEPTH."
  call ESMF_FieldScatter(seaice_depth_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice skin temperature

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tisfc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
  call ESMF_FieldScatter(seaice_skin_temp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! liquid equivalent snow depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('sheleg', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
  call ESMF_FieldScatter(snow_liq_equiv_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! physical snow depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('snwdph', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile = data_one_tile
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
  call ESMF_FieldScatter(snow_depth_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Vegetation type

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('vtype', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID VEG TYPE."
  call ESMF_FieldScatter(veg_type_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Soil type

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('stype', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
  call ESMF_FieldScatter(soil_type_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Two-meter temperature

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('t2m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID T2M."
  call ESMF_FieldScatter(t2m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Two-meter q

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('q2m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
  call ESMF_FieldScatter(q2m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tprcp', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID TPRCP."
  call ESMF_FieldScatter(tprcp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('f10m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID F10M"
  call ESMF_FieldScatter(f10m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('ffmm', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID FFMM"
  call ESMF_FieldScatter(ffmm_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('uustar', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID USTAR"
  call ESMF_FieldScatter(ustar_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('srflag', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
  call ESMF_FieldScatter(srflag_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tsea', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
  call ESMF_FieldScatter(skin_temp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('canopy', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID CANOPY MOISTURE CONTENT."
  call ESMF_FieldScatter(canopy_mc_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('zorl', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID Z0."
  call ESMF_FieldScatter(z0_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 enddo TILE_LOOP

 deallocate(data_one_tile, data_one_tile_3d)

 end subroutine read_input_sfc_restart_file

!---------------------------------------------------------------------------
! Read input grid surface data from tiled 'history' files (netcdf) or 
! gaussian netcdf files.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_netcdf_file(localpet)

 implicit none

 integer, intent(in)             :: localpet

 character(len=500)              :: tilefile

 integer                         :: error, id_var
 integer                         :: id_dim, idim_input, jdim_input
 integer                         :: ncid, rc, tile

 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable :: data_one_tile_3d(:,:,:)

!---------------------------------------------------------------------------
! Get i/j dimensions and number of soil layers from first surface file.
! Do dimensions match those from the orography file?
!---------------------------------------------------------------------------

 tilefile = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))
 print*,"- READ GRID DIMENSIONS FROM: ", trim(tilefile)
 error=nf90_open(trim(tilefile),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(tilefile) )

 error=nf90_inq_dimid(ncid, 'grid_xt', id_dim)
 call netcdf_err(error, 'reading grid_xt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=idim_input)
 call netcdf_err(error, 'reading grid_xt value' )

 error=nf90_inq_dimid(ncid, 'grid_yt', id_dim)
 call netcdf_err(error, 'reading grid_yt id' )
 error=nf90_inquire_dimension(ncid,id_dim,len=jdim_input)
 call netcdf_err(error, 'reading grid_yt value' )

 if (idim_input /= i_input .or. jdim_input /= j_input) then
   call error_handler("DIMENSION MISMATCH BETWEEN SFC AND OROG FILES.", 3)
 endif

 error = nf90_close(ncid)

 if (localpet == 0) then
   allocate(data_one_tile(idim_input,jdim_input))
   allocate(data_one_tile_3d(idim_input,jdim_input,lsoil_input))
 else
   allocate(data_one_tile(0,0))
   allocate(data_one_tile_3d(0,0,0))
 endif

 TERRAIN_LOOP: do tile = 1, num_tiles_input_grid

   if (trim(input_type) == "gaussian_netcdf") then
    if (localpet == 0) then
      call read_fv3_grid_data_netcdf('orog', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    endif

  else
   
   if (localpet == 0) then
     tilefile = trim(orog_dir_input_grid) // trim(orog_files_input_grid(tile))
     print*,'- OPEN OROGRAPHY FILE: ', trim(tilefile)
     error=nf90_open(tilefile,nf90_nowrite,ncid)
     call netcdf_err(error, 'OPENING OROGRAPHY FILE.' )
     error=nf90_inq_varid(ncid, 'orog_raw', id_var)
     call netcdf_err(error, 'READING OROGRAPHY RECORD ID.' )
     error=nf90_get_var(ncid, id_var, data_one_tile)
     call netcdf_err(error, 'READING OROGRAPHY RECORD.' )
     print*,'terrain check history ',tile, maxval(data_one_tile)
     error=nf90_close(ncid)
   endif

   endif

   print*,"- CALL FieldScatter FOR INPUT TERRAIN."
   call ESMF_FieldScatter(terrain_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 enddo TERRAIN_LOOP

 TILE_LOOP : do tile = 1, num_tiles_input_grid

! liquid soil moisture

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('soill1', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,1) = data_one_tile
    call read_fv3_grid_data_netcdf('soill2', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,2) = data_one_tile
    call read_fv3_grid_data_netcdf('soill3', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,3) = data_one_tile
    call read_fv3_grid_data_netcdf('soill4', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,4) = data_one_tile
  endif

  print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
  call ESMF_FieldScatter(soilm_liq_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! total soil moisture

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('soilw1', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,1) = data_one_tile
    call read_fv3_grid_data_netcdf('soilw2', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,2) = data_one_tile
    call read_fv3_grid_data_netcdf('soilw3', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,3) = data_one_tile
    call read_fv3_grid_data_netcdf('soilw4', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,4) = data_one_tile
  endif

  print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
  call ESMF_FieldScatter(soilm_tot_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! soil tempeature (ice temp at land ice points)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('soilt1', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,1) = data_one_tile
    call read_fv3_grid_data_netcdf('soilt2', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,2) = data_one_tile
    call read_fv3_grid_data_netcdf('soilt3', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,3) = data_one_tile
    call read_fv3_grid_data_netcdf('soilt4', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile_3d(:,:,4) = data_one_tile
  endif

  print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
  call ESMF_FieldScatter(soil_temp_input_grid, data_one_tile_3d, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! land mask

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('land', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
  call ESMF_FieldScatter(landsea_mask_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice fraction

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('icec', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
  call ESMF_FieldScatter(seaice_fract_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('icetk', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE DEPTH."
  call ESMF_FieldScatter(seaice_depth_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! sea ice skin temperature

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tisfc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
  call ESMF_FieldScatter(seaice_skin_temp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! liquid equivalent snow depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('weasd', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
  call ESMF_FieldScatter(snow_liq_equiv_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! physical snow depth

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('snod', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
    data_one_tile = data_one_tile * 1000.0  ! convert from meters to mm.
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
  call ESMF_FieldScatter(snow_depth_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Vegetation type

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('vtype', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID VEG TYPE."
  call ESMF_FieldScatter(veg_type_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Soil type

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('sotyp', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
  call ESMF_FieldScatter(soil_type_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Two-meter temperature

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tmp2m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID T2M."
  call ESMF_FieldScatter(t2m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! Two-meter q

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('spfh2m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
  call ESMF_FieldScatter(q2m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tprcp', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID TPRCP."
  call ESMF_FieldScatter(tprcp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('f10m', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID F10M"
  call ESMF_FieldScatter(f10m_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('ffmm', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID FFMM"
  call ESMF_FieldScatter(ffmm_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('fricv', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID USTAR"
  call ESMF_FieldScatter(ustar_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
!   call read_fv3_grid_data_netcdf('srflag', tile, idim_input, jdim_input, &
!                                  lsoil_input, sfcdata=data_one_tile)
    data_one_tile = 0.0
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
  call ESMF_FieldScatter(srflag_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tmpsfc', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
  call ESMF_FieldScatter(skin_temp_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('cnwat', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID CANOPY MOISTURE CONTENT."
  call ESMF_FieldScatter(canopy_mc_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('sfcr', tile, idim_input, jdim_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT GRID Z0."
  call ESMF_FieldScatter(z0_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 enddo TILE_LOOP

 deallocate(data_one_tile, data_one_tile_3d)

 end subroutine read_input_sfc_netcdf_file

!---------------------------------------------------------------------------
! Read surface data from an fv3gfs grib2 file.
!---------------------------------------------------------------------------

 subroutine read_input_sfc_grib2_file(localpet)

   use wgrib2api

   implicit none

   integer, intent(in)                   :: localpet

   character(len=250)                    :: the_file
   character(len=20)                     :: vname, vname_file,slev

   character(len=50)                      :: method

   integer                               :: rc, varnum, iret, i, j,k
   integer, parameter                    :: icet_default = 265.0

   logical                               :: exist

   real(esmf_kind_r4)                    :: value

   real(esmf_kind_r4), allocatable       :: dummy2d(:,:),tsk_save(:,:),icec_save(:,:)
   real(esmf_kind_r8), allocatable       :: dummy2d_8(:,:)
   real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:)
   integer(esmf_kind_i4), allocatable    :: slmsk_save(:,:)
   

   the_file = trim(data_dir_input_grid) // "/" // trim(grib2_file_input_grid)

   print*,"- READ SFC DATA FROM GRIB2 FILE: ", trim(the_file)
   inquire(file=the_file,exist=exist)
   if (.not.exist) then
     iret = 1
     call error_handler("OPENING GRIB2 FILE.", iret)
   end if

   lsoil_input = grb2_inq(the_file, inv_file, ':TSOIL:',' below ground:')
   print*, "- FILE HAS ", lsoil_input, " SOIL LEVELS"
   if (lsoil_input <= 0) call error_handler("COUNTING SOIL LEVELS.", rc)

 if (localpet == 0) then
   allocate(dummy2d(i_input,j_input))
   allocate(slmsk_save(i_input,j_input))
   allocate(tsk_save(i_input,j_input))
   allocate(icec_save(i_input,j_input))
   allocate(dummy2d_8(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lsoil_input))
 else
   allocate(dummy3d(0,0,0))
   allocate(dummy2d_8(0,0))
   allocate(dummy2d(0,0))

 endif
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! These variables are always in grib files, or are required, so no need to check for them 
 ! in the varmap table. If they can't be found in the input file, then stop the program.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (localpet == 0) then
   print*,"- READ TERRAIN."
   rc = grb2_inq(the_file, inv_file, ':HGT:',':surface:', data2=dummy2d)
   if (rc /= 1) call error_handler("READING TERRAIN.", rc)
   print*,'orog ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
if (localpet == 0) then
   print*,"- READ SEAICE FRACTION."
   rc = grb2_inq(the_file, inv_file, ':ICEC:',':surface:', data2=dummy2d)
   if (rc /= 1) call error_handler("READING SEAICE FRACTION.", rc)
   !dummy2d = dummy2d(i_input:1:-1,j_input:1:-1)
   print*,'icec ',maxval(dummy2d),minval(dummy2d)
   icec_save = dummy2d
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
 call ESMF_FieldScatter(seaice_fract_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

!----------------------------------------------------------------------------------
! GFS v14 and v15.2 grib data has two land masks.  LANDN is created by
! nearest neighbor interpolation.  LAND is created by bilinear interpolation.
! LANDN matches the bitmap.  So use it first.  For other GFS versions, use LAND.
! Mask in grib file is '1' (land), '0' (not land).  Add sea/lake ice category
! '2' based on ice concentration.
!----------------------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- READ LANDSEA MASK."
   rc = grb2_inq(the_file, inv_file, ':LANDN:',':surface:', data2=dummy2d)

   if (rc /= 1) then 
     rc = grb2_inq(the_file, inv_file, ':LAND:',':surface:', data2=dummy2d)
     if (rc /= 1) call error_handler("READING LANDSEA MASK.", rc)
   endif

   do j = 1, j_input
     do i = 1, i_input
       if(dummy2d(i,j) < 0.5_esmf_kind_r4) dummy2d(i,j)=0.0_esmf_kind_r4
       if(icec_save(i,j) > 0.15_esmf_kind_r4) then 
         !if (dummy2d(i,j) == 0.0_esmf_kind_r4) print*, "CONVERTING WATER TO SEA/LAKE ICE AT ", i, j
         dummy2d(i,j) = 2.0_esmf_kind_r4
       endif
     enddo
   enddo

   slmsk_save = nint(dummy2d)
  
   deallocate(icec_save)
 endif

 print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
 call ESMF_FieldScatter(landsea_mask_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SEAICE SKIN TEMPERATURE."
   rc = grb2_inq(the_file, inv_file, ':TMP:',':surface:', data2=dummy2d)
   if (rc /= 1) call error_handler("READING SEAICE SKIN TEMP.", rc)
   print*,'ti ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
 call ESMF_FieldScatter(seaice_skin_temp_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

!----------------------------------------------------------------------------------
! Read snow fields.  Zero out at non-land points and undefined points (points
! removed using the bitmap).  Program expects depth and liquid equivalent
! in mm.
!----------------------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- READ SNOW LIQUID EQUIVALENT."
   rc = grb2_inq(the_file, inv_file, ':WEASD:',':surface:',':anl:',data2=dummy2d)
   if (rc /= 1) then 
     rc = grb2_inq(the_file, inv_file, ':WEASD:',':surface:','hour fcst:',data2=dummy2d)
     if (rc /= 1) call error_handler("READING SNOW LIQUID EQUIVALENT.", rc)
   endif
   do j = 1, j_input
     do i = 1, i_input
       if(slmsk_save(i,j) == 0) dummy2d(i,j) = 0.0_esmf_kind_r4
       if(dummy2d(i,j) == grb2_UNDEFINED) dummy2d(i,j) = 0.0_esmf_kind_r4
     enddo
   enddo
!  print*,'weasd ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
 call ESMF_FieldScatter(snow_liq_equiv_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ SNOW DEPTH."
   rc = grb2_inq(the_file, inv_file, ':SNOD:',':surface:', data2=dummy2d)
   if (rc /= 1) call error_handler("READING SNOW DEPTH.", rc)
   where(dummy2d == grb2_UNDEFINED) dummy2d = 0.0_esmf_kind_r4
   dummy2d = dummy2d*1000.0 ! Grib2 files have snow depth in (m), fv3 expects it in mm
   where(slmsk_save == 0) dummy2d = 0.0_esmf_kind_r4
!  print*,'snod ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
 call ESMF_FieldScatter(snow_depth_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 if (localpet == 0) then
   print*,"- READ T2M."
   rc = grb2_inq(the_file, inv_file, ':TMP:',':2 m above ground:',data2=dummy2d)
   if (rc <= 0) call error_handler("READING T2M.", rc)

   print*,'t2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID T2M."
 call ESMF_FieldScatter(t2m_input_grid,real(dummy2d,esmf_kind_r8), rootpet=0,rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Q2M."
   rc = grb2_inq(the_file, inv_file, ':SPFH:',':2 m above ground:',data2=dummy2d)
   if (rc <=0) call error_handler("READING Q2M.", rc)
   print*,'q2m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
 call ESMF_FieldScatter(q2m_input_grid,real(dummy2d,esmf_kind_r8), rootpet=0,rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 if (localpet == 0) then
   print*,"- READ SKIN TEMPERATURE."
   rc = grb2_inq(the_file, inv_file, ':TMP:',':surface:', data2=dummy2d)
   if (rc <= 0 ) call error_handler("READING SKIN TEMPERATURE.", rc)
   tsk_save(:,:) = dummy2d
   dummy2d_8 = real(dummy2d,esmf_kind_r8)
   do j = 1, j_input
     do i = 1, i_input
       if(slmsk_save(i,j) == 0 .and. dummy2d(i,j) < 271.2) then
!        print*,'too cool SST ',i,j,dummy2d(i,j)
         dummy2d(i,j) = 271.2
       endif
       if(slmsk_save(i,j) == 0 .and. dummy2d(i,j) > 310.) then
!        print*,'too hot SST ',i,j,dummy2d(i,j)
         dummy2d(i,j) = 310.0
       endif
     enddo
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
 call ESMF_FieldScatter(skin_temp_input_grid,real(dummy2d,esmf_kind_r8),rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 if (localpet == 0) dummy2d = 0.0
 
 print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
 call ESMF_FieldScatter(srflag_input_grid,real(dummy2d,esmf_kind_r8), rootpet=0,rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

! Soil type is not available.  Set to a large negative fill value.

 dummy2d_8 = -99999.0_esmf_kind_r8

 print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
 call ESMF_FieldScatter(soil_type_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
 ! Begin variables whose presence in grib2 files varies, but no climatological
 ! data is 
 ! available, so we have to account for values in the varmap table
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (localpet == 0) then
   print*,"- READ SEAICE DEPTH."
   vname="hice"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)                 
   vname=":ICETK:"
   rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
   if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL BE"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'hice ',maxval(dummy2d),minval(dummy2d)

 endif

 print*,"- CALL FieldScatter FOR INPUT GRID SEAICE DEPTH."
 call ESMF_FieldScatter(seaice_depth_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 if (localpet == 0) then
   print*,"- READ TPRCP."
   vname="tprcp"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":TPRCP:"              
   rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
   if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                   " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'tprcp ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID TPRCP."
 call ESMF_FieldScatter(tprcp_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ FFMM."
   vname="ffmm"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":FFMM:"               
    rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
    if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                   " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'ffmm ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID FFMM"
 call ESMF_FieldScatter(ffmm_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 if (localpet == 0) then
   print*,"- READ USTAR."
   vname="fricv"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":FRICV:"              
    rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
    if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL "//&
                   "REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'fricv ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID USTAR"
 call ESMF_FieldScatter(ustar_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ F10M."
   vname="f10m"
   slev=":10 m above ground:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":F10M:"               
    rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
    if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                   " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'f10m ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID F10M."
 call ESMF_FieldScatter(f10m_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ CANOPY MOISTURE CONTENT."
   vname="cnwat"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":CNWAT:"              
    rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
    if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'cnwat ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID CANOPY MOISTURE CONTENT."
 call ESMF_FieldScatter(canopy_mc_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ Z0."
   vname="sfcr"
   slev=":surface:" 
   call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)  
    vname=":SFCR:"               
    rc= grb2_inq(the_file, inv_file, vname,slev, data2=dummy2d)
    if (rc <= 0) then
      call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
      if (rc==1) then ! missing_var_method == skip or no entry in varmap table
        print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL BE"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
        dummy2d(:,:) = 0.0_esmf_kind_r4
      endif
    else
      ! Grib files have z0 (m), but fv3 expects z0(cm)
      dummy2d(:,:) = dummy2d(:,:)*10.0
    endif
   dummy2d_8= real(dummy2d,esmf_kind_r8)
   print*,'sfcr ',maxval(dummy2d),minval(dummy2d)
   
 endif

 print*,"- CALL FieldScatter FOR INPUT GRID Z0."
 call ESMF_FieldScatter(z0_input_grid,dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
    
 deallocate(dummy2d)
 
 if (localpet == 0) then
   print*,"- READ LIQUID SOIL MOISTURE."
   vname = "soill"
   vname_file = ":SOILL:"
   call read_grib_soil(the_file,inv_file,vname,vname_file,dummy3d,rc) !!! NEEDTO HANDLE 
                                                                      !!! SOIL LEVELS
   print*,'soill ',maxval(dummy3d),minval(dummy3d)
 endif

 print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_liq_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)
 
 if (localpet == 0) then
   print*,"- READ TOTAL SOIL MOISTURE."
   vname = "soilw"
   !vname_file = "var2_2_1_7_0_192"  !Some files don't recognize this as soilw,so use
   vname_file = "var2_2_1_"         ! the var number instead
   call read_grib_soil(the_file,inv_file,vname,vname_file,dummy3d,rc)
   print*,'soilm ',maxval(dummy3d),minval(dummy3d)
 endif

!-----------------------------------------------------------------------
! Vegetation type is not available.  However, it is needed to identify
! permanent land ice points.  At land ice, the total soil moisture
! is a flag value of '1'.  Use this flag as a temporary solution.
!-----------------------------------------------------------------------

 if (localpet == 0) then
   dummy2d_8(:,:) = 0.0_esmf_kind_r8
   do j = 1, j_input
       do i = 1, i_input
         if(slmsk_save(i,j) == 1_esmf_kind_i4 .and. dummy3d(i,j,1) > 0.99) &
            dummy2d_8(i,j) = real(veg_type_landice_input,esmf_kind_r8)
       enddo
   enddo
 endif

 print*,"- CALL FieldScatter FOR INPUT VEG TYPE."
 call ESMF_FieldScatter(veg_type_input_grid, dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
 call ESMF_FieldScatter(soilm_tot_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

!---------------------------------------------------------------------------------
! At open water (slmsk==0), the soil temperature array is not used and set
! to the filler value of SST.  At lake/sea ice points (slmsk=2), the soil 
! temperature array holds ice column temperature.  That field is not available
! in GFS grib data, so set to a default value.
!---------------------------------------------------------------------------------

 if (localpet == 0) then
   print*,"- READ SOIL TEMPERATURE."
   vname = "soilt"
   vname_file = ":TSOIL:"
   call read_grib_soil(the_file,inv_file,vname,vname_file,dummy3d,rc)
   do k=1,lsoil_input
     do j = 1, j_input
       do i = 1, i_input
         if (slmsk_save(i,j) == 0_esmf_kind_i4 ) dummy3d(i,j,k) = tsk_save(i,j)
         if (slmsk_save(i,j) == 2_esmf_kind_i4 ) dummy3d(i,j,k) = icet_default
       enddo
     enddo
   enddo
   print*,'soilt ',maxval(dummy3d),minval(dummy3d)

   deallocate(tsk_save, slmsk_save)
 endif

 print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
 call ESMF_FieldScatter(soil_temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy3d)
 deallocate(dummy2d_8)
 
 end subroutine read_input_sfc_grib2_file
   
!---------------------------------------------------------------------------
! Read nst data from these netcdf formatted fv3 files: tiled history,
! tiled warm restart, and gaussian history.
!---------------------------------------------------------------------------

 subroutine read_input_nst_netcdf_file(localpet)

 implicit none

 integer, intent(in)             :: localpet

 character(len=10)               :: field

 integer                         :: rc, tile

 real(esmf_kind_r8), allocatable :: data_one_tile(:,:)

 if (localpet == 0) then
   allocate(data_one_tile(i_input,j_input))
 else
   allocate(data_one_tile(0,0))
 endif

 TILE_LOOP : do tile = 1, num_tiles_input_grid

! c_d

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='c_d'
    else
      field='cd'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT C_D"
  call ESMF_FieldScatter(c_d_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! c_0

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='c_0'
    else
      field='c0'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT C_0"
  call ESMF_FieldScatter(c_0_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! d_conv

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='d_conv'
    else
      field='dconv'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT D_CONV."
  call ESMF_FieldScatter(d_conv_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! dt_cool

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='dt_cool'
    else
      field='dtcool'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT DT_COOL."
  call ESMF_FieldScatter(dt_cool_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! ifd - xu li said initialize to '1'.

  if (localpet == 0) then
    data_one_tile = 1.0
  endif

  print*,"- CALL FieldScatter FOR INPUT IFD."
  call ESMF_FieldScatter(ifd_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! qrain

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('qrain', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT QRAIN."
  call ESMF_FieldScatter(qrain_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! tref

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('tref', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT TREF"
  call ESMF_FieldScatter(tref_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! w_d

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='w_d'
    else
      field='wd'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT W_D"
  call ESMF_FieldScatter(w_d_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! w_0

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='w_0'
    else
      field='w0'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT W_0"
  call ESMF_FieldScatter(w_0_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xs

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xs', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XS"
  call ESMF_FieldScatter(xs_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xt

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xt', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XT"
  call ESMF_FieldScatter(xt_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xu

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xu', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XU"
  call ESMF_FieldScatter(xu_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xv

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xv', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XV"
  call ESMF_FieldScatter(xv_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xz

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xz', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XZ"
  call ESMF_FieldScatter(xz_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xtts

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xtts', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XTTS"
  call ESMF_FieldScatter(xtts_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! xzts

  if (localpet == 0) then
    call read_fv3_grid_data_netcdf('xzts', tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT XZTS"
  call ESMF_FieldScatter(xzts_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! z_c

  if (localpet == 0) then
    if (trim(input_type) == "restart") then
      field='z_c'
    else
      field='zc'
    endif
    call read_fv3_grid_data_netcdf(trim(field), tile, i_input, j_input, &
                                   lsoil_input, sfcdata=data_one_tile)
  endif

  print*,"- CALL FieldScatter FOR INPUT Z_C"
  call ESMF_FieldScatter(z_c_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

! zm - Not used yet. Xu li said set to '0'.

  if (localpet == 0) then
    data_one_tile = 0.0
  endif

  print*,"- CALL FieldScatter FOR INPUT ZM"
  call ESMF_FieldScatter(zm_input_grid, data_one_tile, rootpet=0, tile=tile, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldScatter", rc)

 enddo TILE_LOOP

 deallocate(data_one_tile)

 end subroutine read_input_nst_netcdf_file

!--------------------------------------------------------------------------
! Read input grid nst data from fv3 gaussian nemsio history file or
! spectral GFS nemsio file.  The spectral GFS nst data is in a separate
! file from the surface data.  The fv3 surface and nst data are in a
! single file.
!--------------------------------------------------------------------------

 subroutine read_input_nst_nemsio_file(localpet)

 implicit none

 integer, intent(in)                    :: localpet

 character(len=300)                     :: the_file

 integer                                :: rc

 real(nemsio_realkind), allocatable     :: dummy(:)
 real(esmf_kind_r8), allocatable        :: dummy2d(:,:)

 type(nemsio_gfile)                     :: gfile

 if (trim(input_type) == "gfs_gaussian_nemsio") then ! spectral gfs nemsio in
                                                     ! separate file.
   the_file = trim(data_dir_input_grid) // "/" // trim(nst_files_input_grid)
 else
   the_file = trim(data_dir_input_grid) // "/" // trim(sfc_files_input_grid(1))
 endif

 print*,"- READ NST DATA FROM: ", trim(the_file)

 if (localpet == 0) then
   allocate(dummy(i_input*j_input))
   allocate(dummy2d(i_input,j_input))
   call nemsio_open(gfile, the_file, "read", iret=rc)
 else
   allocate(dummy(0))
   allocate(dummy2d(0,0))
 endif

 if (localpet == 0) then
   print*,"- READ TREF"
   call nemsio_readrecv(gfile, "tref", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING TREF.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'tref ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT TREF."
 call ESMF_FieldScatter(tref_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ CD"
   call nemsio_readrecv(gfile, "cd", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING CD.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'cd ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT C_D."
 call ESMF_FieldScatter(c_d_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ C0"
   call nemsio_readrecv(gfile, "c0", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING C0.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'c0 ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT C_0."
 call ESMF_FieldScatter(c_0_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ DCONV"
   call nemsio_readrecv(gfile, "dconv", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING DCONV.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'dconv ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT D_CONV."
 call ESMF_FieldScatter(d_conv_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ DTCOOL"
   call nemsio_readrecv(gfile, "dtcool", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING DTCOOL.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'dtcool ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT DT_COOL."
 call ESMF_FieldScatter(dt_cool_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   dummy2d = 1.0  ! IFD not in file.  Set to '1' per Xu Li.
 endif

 print*,"- CALL FieldScatter FOR INPUT IFD."
 call ESMF_FieldScatter(ifd_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ QRAIN"
   call nemsio_readrecv(gfile, "qrain", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING QRAIN.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'qrain ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT QRAIN."
 call ESMF_FieldScatter(qrain_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ WD"
   call nemsio_readrecv(gfile, "wd", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING WD.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'wd ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT WD."
 call ESMF_FieldScatter(w_d_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ W0"
   call nemsio_readrecv(gfile, "w0", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING W0.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'w0 ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT W0."
 call ESMF_FieldScatter(w_0_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XS"
   call nemsio_readrecv(gfile, "xs", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XS.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xs ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XS."
 call ESMF_FieldScatter(xs_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XT"
   call nemsio_readrecv(gfile, "xt", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XT.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xt ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XT."
 call ESMF_FieldScatter(xt_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XU"
   call nemsio_readrecv(gfile, "xu", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XU.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xu ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XU."
 call ESMF_FieldScatter(xu_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XV"
   call nemsio_readrecv(gfile, "xv", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XV.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xv ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XV."
 call ESMF_FieldScatter(xv_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XZ"
   call nemsio_readrecv(gfile, "xz", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XZ.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xz ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XZ."
 call ESMF_FieldScatter(xz_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XTTS"
   call nemsio_readrecv(gfile, "xtts", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XTTS.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xtts ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XTTS."
 call ESMF_FieldScatter(xtts_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ XZTS"
   call nemsio_readrecv(gfile, "xzts", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING XZTS.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'xzts ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT XZTS."
 call ESMF_FieldScatter(xzts_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ ZC"
   call nemsio_readrecv(gfile, "zc", "sfc", 1, dummy, 0, iret=rc)
   if (rc /= 0) call error_handler("READING ZC.", rc)
   dummy2d = reshape(dummy, (/i_input,j_input/))
   print*,'zc ',maxval(dummy2d),minval(dummy2d)
 endif

 print*,"- CALL FieldScatter FOR INPUT Z_C."
 call ESMF_FieldScatter(z_c_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   dummy2d = 0.0 ! zm not used yet. Set to zero per Xu Li.
 endif

 print*,"- CALL FieldScatter FOR INPUT ZM."
 call ESMF_FieldScatter(zm_input_grid, dummy2d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy, dummy2d)

 if (localpet == 0) call nemsio_close(gfile)

 end subroutine read_input_nst_nemsio_file

 SUBROUTINE READ_FV3_GRID_DATA_NETCDF(FIELD,TILE_NUM,IMO,JMO,LMO, &
                                      SFCDATA, SFCDATA_3D)

 IMPLICIT NONE

 CHARACTER(LEN=*),INTENT(IN)      :: FIELD

 INTEGER, INTENT(IN)   :: IMO, JMO, LMO, TILE_NUM

 REAL(ESMF_KIND_R8), INTENT(OUT), OPTIONAL     :: SFCDATA(IMO,JMO)
 REAL(ESMF_KIND_R8), INTENT(OUT), OPTIONAL     :: SFCDATA_3D(IMO,JMO,LMO)

 CHARACTER(LEN=256)    :: TILEFILE

 INTEGER               :: ERROR, NCID, ID_VAR

 TILEFILE = TRIM(DATA_DIR_INPUT_GRID) // "/" // TRIM(SFC_FILES_INPUT_GRID(TILE_NUM))

 PRINT*,'WILL READ ',TRIM(FIELD), ' FROM: ', TRIM(TILEFILE)

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERR(ERROR, 'OPENING: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_VARID(NCID, FIELD, ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING FIELD ID' )

 IF (PRESENT(SFCDATA_3D)) THEN
   ERROR=NF90_GET_VAR(NCID, ID_VAR, SFCDATA_3D)
   CALL NETCDF_ERR(ERROR, 'READING FIELD' )
 ELSE
   ERROR=NF90_GET_VAR(NCID, ID_VAR, SFCDATA)
   CALL NETCDF_ERR(ERROR, 'READING FIELD' )
 ENDIF

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_FV3_GRID_DATA_NETCDF

!---------------------------------------------------------------------------
! Convert from 2-d to 3-d winds.
!---------------------------------------------------------------------------

 subroutine convert_winds

 implicit none

 integer                         :: clb(4), cub(4)
 integer                         :: i, j, k, rc

 real(esmf_kind_r8)              :: latrad, lonrad
 real(esmf_kind_r8), pointer     :: windptr(:,:,:,:)
 real(esmf_kind_r8), pointer     :: uptr(:,:,:)
 real(esmf_kind_r8), pointer     :: vptr(:,:,:)
 real(esmf_kind_r8), pointer     :: latptr(:,:)
 real(esmf_kind_r8), pointer     :: lonptr(:,:)

 print*,"- CALL FieldGet FOR 3-D WIND."
 call ESMF_FieldGet(wind_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=windptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR U."
 call ESMF_FieldGet(u_input_grid, &
                    farrayPtr=uptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR V."
 call ESMF_FieldGet(v_input_grid, &
                    farrayPtr=vptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LATITUDE."
 call ESMF_FieldGet(latitude_input_grid, &
                    farrayPtr=latptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR LONGITUDE."
 call ESMF_FieldGet(longitude_input_grid, &
                    farrayPtr=lonptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 do i = clb(1), cub(1)
   do j = clb(2), cub(2)
     latrad = latptr(i,j) * acos(-1.) / 180.0
     lonrad = lonptr(i,j) * acos(-1.) / 180.0
     do k = clb(3), cub(3)
       windptr(i,j,k,1) = uptr(i,j,k) * cos(lonrad) - vptr(i,j,k) * sin(latrad) * sin(lonrad)
       windptr(i,j,k,2) = uptr(i,j,k) * sin(lonrad) + vptr(i,j,k) * sin(latrad) * cos(lonrad)
       windptr(i,j,k,3) = vptr(i,j,k) * cos(latrad)
     enddo
   enddo
 enddo

 call ESMF_FieldDestroy(u_input_grid, rc=rc)
 call ESMF_FieldDestroy(v_input_grid, rc=rc)

 end subroutine convert_winds
 
subroutine handle_grib_error(vname,lev,method,value,varnum, iret,var,var8,var3d)

  use, intrinsic :: ieee_arithmetic

  implicit none
  
  real(esmf_kind_r4), intent(in)    :: value
  real(esmf_kind_r4), intent(inout), optional :: var(:,:)
  real(esmf_kind_r8), intent(inout), optional :: var8(:,:)
  real(esmf_kind_r8), intent(inout), optional  :: var3d(:,:,:)
  
  character(len=20), intent(in)     :: vname, lev, method
  
  integer, intent(in)               :: varnum 
  integer, intent(inout)            :: iret
  
  iret = 0
  if (varnum == 9999) then
    print*, "WARNING: ", trim(vname), " NOT FOUND AT LEVEL ", lev, " IN EXTERNAL FILE ", &
            "AND NO ENTRY EXISTS IN VARMAP TABLE. VARIABLE WILL NOT BE USED."
    iret = 1

    return
  endif

  if (trim(method) == "skip" ) then
    print*, "WARNING: SKIPPING ", trim(vname), " IN FILE"
    read_from_input(varnum) = .false.
    iret = 1
  elseif (trim(method) == "set_to_fill") then
    print*, "WARNING: ,", trim(vname), " NOT AVILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO FILL VALUE OF ", value
    if(present(var)) var(:,:) = value
    if(present(var8)) var8(:,:) = value
    if(present(var3d)) var3d(:,:,:) = value
  elseif (trim(method) == "set_to_NaN") then
    print*, "WARNING: ,", trim(vname), " NOT AVILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO NaNs"
    if(present(var)) var(:,:) = ieee_value(var,IEEE_QUIET_NAN)
    if(present(var8)) var8(:,:) = ieee_value(var8,IEEE_QUIET_NAN)
    if(present(var3d)) var3d(:,:,:) = ieee_value(var3d,IEEE_QUIET_NAN)
  elseif (trim(method) == "stop") then
    call error_handler("READING "//trim(vname)// " at level "//lev//". TO MAKE THIS NON- &
                        FATAL, CHANGE STOP TO SKIP FOR THIS VARIABLE IN YOUR VARMAP &
                        FILE.", iret)
  else
    call error_handler("ERROR USING MISSING_VAR_METHOD. PLEASE SET VALUES IN" // &
                       " VARMAP TABLE TO ONE OF: set_to_fill, set_to_NaN,"// &
                       " , skip, or stop.", 1)
  endif

end subroutine handle_grib_error

subroutine read_grib_soil(the_file,inv_file,vname,vname_file,dummy3d,rc)
  
  use wgrib2api
  implicit none
  
  character(len=*), intent(in)            :: the_file, inv_file
  character(len=20), intent(in)           :: vname,vname_file
  
  integer, intent(out)                    :: rc
  
  real(esmf_kind_r8), intent(inout)       :: dummy3d(:,:,:)
  
  real(esmf_kind_r4), allocatable         :: dummy2d(:,:)
  real(esmf_kind_r4)                      :: value
  integer                                 :: varnum,i
  character(len=50)                       :: slevs(lsoil_input)
  character(len=50)                       :: method

  allocate(dummy2d(i_input,j_input))

  if(lsoil_input == 4) then
    slevs = (/character(24)::':0-0.1 m below ground:', ':0.1-0.4 m below ground:', &
                             ':0.4-1 m below ground:', ':1-2 m below ground:'/)
  else
    rc = -1
    call error_handler("reading soil levels. File must have 4 soil levels.", rc)
  endif
 
  call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)
  do i = 1,lsoil_input
    if (vname_file=="var2_2_1_") then
      rc = grb2_inq(the_file,inv_file,vname_file,"_0_192:",slevs(i),data2=dummy2d)
    else
      rc = grb2_inq(the_file,inv_file,vname_file,slevs(i),data2=dummy2d)
    endif
    if (rc <= 0) then
      call handle_grib_error(vname_file, slevs(i),method,value,varnum,rc,var=dummy2d)
      if (rc==1 .and. trim(vname) /= "soill") then 
        ! missing_var_method == skip or no entry in varmap table
        call error_handler("READING IN "//trim(vname)//". SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",rc)
      elseif (rc==1) then
        dummy3d(:,:,:) = 0.0_esmf_kind_r8
        exit
      endif
    endif

    dummy3d(:,:,i) = real(dummy2d,esmf_kind_r8)
  end do    

 deallocate(dummy2d)

 end subroutine read_grib_soil

 subroutine cleanup_input_atm_data

 implicit none

 integer                         :: rc, n

 print*,'- DESTROY ATMOSPHERIC INPUT DATA.'

 call ESMF_FieldDestroy(terrain_input_grid, rc=rc)
 call ESMF_FieldDestroy(pres_input_grid, rc=rc)
 call ESMF_FieldDestroy(dzdt_input_grid, rc=rc)
 call ESMF_FieldDestroy(temp_input_grid, rc=rc)
 call ESMF_FieldDestroy(wind_input_grid, rc=rc)
 call ESMF_FieldDestroy(ps_input_grid, rc=rc)

 do n = 1, num_tracers
   call ESMF_FieldDestroy(tracers_input_grid(n), rc=rc)
 enddo
 deallocate(tracers_input_grid)

 end subroutine cleanup_input_atm_data

 subroutine cleanup_input_nst_data

 implicit none

 integer                         :: rc

 print*,'- DESTROY NST INPUT DATA.'

 call ESMF_FieldDestroy(landsea_mask_input_grid, rc=rc)
 call ESMF_FieldDestroy(c_d_input_grid, rc=rc)
 call ESMF_FieldDestroy(c_0_input_grid, rc=rc)
 call ESMF_FieldDestroy(d_conv_input_grid, rc=rc)
 call ESMF_FieldDestroy(dt_cool_input_grid, rc=rc)
 call ESMF_FieldDestroy(ifd_input_grid, rc=rc)
 call ESMF_FieldDestroy(qrain_input_grid, rc=rc)
 call ESMF_FieldDestroy(tref_input_grid, rc=rc)
 call ESMF_FieldDestroy(w_d_input_grid, rc=rc)
 call ESMF_FieldDestroy(w_0_input_grid, rc=rc)
 call ESMF_FieldDestroy(xs_input_grid, rc=rc)
 call ESMF_FieldDestroy(xt_input_grid, rc=rc)
 call ESMF_FieldDestroy(xu_input_grid, rc=rc)
 call ESMF_FieldDestroy(xv_input_grid, rc=rc)
 call ESMF_FieldDestroy(xz_input_grid, rc=rc)
 call ESMF_FieldDestroy(xtts_input_grid, rc=rc)
 call ESMF_FieldDestroy(xzts_input_grid, rc=rc)
 call ESMF_FieldDestroy(z_c_input_grid, rc=rc)
 call ESMF_FieldDestroy(zm_input_grid, rc=rc)

 end subroutine cleanup_input_nst_data

 subroutine cleanup_input_sfc_data

 implicit none

 integer                         :: rc

 print*,"- CALL FieldDestroy FOR INPUT GRID FIELDS."

 call ESMF_FieldDestroy(canopy_mc_input_grid, rc=rc)
 call ESMF_FieldDestroy(f10m_input_grid, rc=rc)
 call ESMF_FieldDestroy(ffmm_input_grid, rc=rc)
 if (.not. convert_nst) then
   call ESMF_FieldDestroy(landsea_mask_input_grid, rc=rc)
 endif
 call ESMF_FieldDestroy(q2m_input_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_depth_input_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_fract_input_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_skin_temp_input_grid, rc=rc)
 call ESMF_FieldDestroy(skin_temp_input_grid, rc=rc)
 call ESMF_FieldDestroy(snow_depth_input_grid, rc=rc)
 call ESMF_FieldDestroy(snow_liq_equiv_input_grid, rc=rc)
 call ESMF_FieldDestroy(soil_temp_input_grid, rc=rc)
 call ESMF_FieldDestroy(soil_type_input_grid, rc=rc)
 call ESMF_FieldDestroy(soilm_liq_input_grid, rc=rc)
 call ESMF_FieldDestroy(soilm_tot_input_grid, rc=rc)
 call ESMF_FieldDestroy(srflag_input_grid, rc=rc)
 call ESMF_FieldDestroy(t2m_input_grid, rc=rc)
 call ESMF_FieldDestroy(tprcp_input_grid, rc=rc)
 call ESMF_FieldDestroy(ustar_input_grid, rc=rc)
 call ESMF_FieldDestroy(veg_type_input_grid, rc=rc)
 call ESMF_FieldDestroy(z0_input_grid, rc=rc)
 call ESMF_FieldDestroy(terrain_input_grid, rc=rc)

 end subroutine cleanup_input_sfc_data

! Jili Dong add sort subroutine 
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

 end module input_data
