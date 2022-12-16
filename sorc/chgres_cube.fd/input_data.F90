!> @file
!! @brief Read atmospheric and surface data from GRIB2, NEMSIO and NetCDF files.
!! @author George Gayno NCEP/EMC

!> Read atmospheric, surface and nst data on the input grid.
!! Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!! 'history' files, fv3 gaussian history files, spectral gfs
!! gaussian nemsio files, and spectral gfs sigio/sfcio files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author George Gayno NCEP/EMC
 module input_data

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
                                    tracers_input, num_tracers_input, &
                                    input_type, tracers, &
                                    get_var_cond, read_from_input, &
                                    geogrid_file_input_grid, &
                                    external_model, &
                                    vgfrc_from_climo, &
                                    minmax_vgfrc_from_climo, &
                                    lai_from_climo

 use model_grid, only             : input_grid,        &
                                    i_input, j_input,  &
                                    ip1_input, jp1_input,  &
                                    num_tiles_input_grid, &
                                    latitude_input_grid, &
                                    longitude_input_grid

 implicit none

 private

! Fields associated with the atmospheric model.

 type(esmf_field), public              :: dzdt_input_grid       !< vert velocity
 type(esmf_field)                      :: dpres_input_grid      !< pressure thickness
 type(esmf_field), public              :: pres_input_grid       !< 3-d pressure
 type(esmf_field), public              :: ps_input_grid         !< surface pressure
 type(esmf_field), public              :: terrain_input_grid    !< terrain height
 type(esmf_field), public              :: temp_input_grid       !< temperature

 type(esmf_field), public              :: u_input_grid          !< u/v wind at grid
 type(esmf_field), public              :: v_input_grid          !< box center
 type(esmf_field), public              :: wind_input_grid       !< 3-component wind
 type(esmf_field), allocatable, public :: tracers_input_grid(:) !< tracers

 integer, public                 :: lev_input      !< number of atmospheric layers
 integer, public                 :: levp1_input    !< number of atmos layer interfaces

! Fields associated with the land-surface model.

 integer, public                 :: veg_type_landice_input = 15 !< NOAH land ice option
                                                                !< defined at this veg type.
                                                                !< Default is igbp.
 integer, parameter              :: ICET_DEFAULT = 265.0    !< Default value of soil and skin
                                                            !< temperature (K) over ice.
 type(esmf_field), public :: canopy_mc_input_grid    !< canopy moist content
 type(esmf_field), public :: f10m_input_grid         !< log((z0+10)*1/z0)
 type(esmf_field), public :: ffmm_input_grid         !< log((z0+z1)*1/z0)
                                                            !! See sfc_diff.f for details.
 type(esmf_field), public :: landsea_mask_input_grid !< land sea mask;
                                                            !! 0-water, 1-land, 2-ice
 type(esmf_field), public :: q2m_input_grid          !< 2-m spec hum
 type(esmf_field), public :: seaice_depth_input_grid !< sea ice depth
 type(esmf_field), public :: seaice_fract_input_grid !< sea ice fraction
 type(esmf_field), public :: seaice_skin_temp_input_grid  !< sea ice skin temp
 type(esmf_field), public :: skin_temp_input_grid    !< skin temp/sst
 type(esmf_field), public :: snow_depth_input_grid   !< snow dpeth
 type(esmf_field), public :: snow_liq_equiv_input_grid !< snow liq equiv depth
 type(esmf_field), public :: soil_temp_input_grid    !< 3-d soil temp
 type(esmf_field), public :: soil_type_input_grid    !< soil type
 type(esmf_field), public :: soilm_liq_input_grid    !< 3-d liquid soil moisture
 type(esmf_field), public :: soilm_tot_input_grid    !< 3-d total soil moisture
 type(esmf_field), public :: srflag_input_grid       !< snow/rain flag
 type(esmf_field), public :: t2m_input_grid          !< 2-m temperature
 type(esmf_field), public :: tprcp_input_grid        !< precip
 type(esmf_field), public :: ustar_input_grid        !< fric velocity
 type(esmf_field), public :: veg_type_input_grid     !< vegetation type
 type(esmf_field), public :: z0_input_grid           !< roughness length
 type(esmf_field), public :: veg_greenness_input_grid !< vegetation fraction
 type(esmf_field), public :: lai_input_grid          !< leaf area index
 type(esmf_field), public :: max_veg_greenness_input_grid !< shdmax
 type(esmf_field), public :: min_veg_greenness_input_grid !< shdmin

 integer, public      :: lsoil_input=4  !< number of soil layers, no longer hardwired to allow
                                        !! for 7 layers of soil for the RUC LSM
 
 character(len=50), private, allocatable :: slevs(:) !< The atmospheric levels in the GRIB2 input file.

! Fields associated with the nst model.

 type(esmf_field), public :: c_d_input_grid   !< Coefficient 2 to calculate d(tz)/d(ts)
 type(esmf_field), public :: c_0_input_grid   !< Coefficient 1 to calculate d(tz)/d(ts)
 type(esmf_field), public :: d_conv_input_grid   !< Thickness of free convection layer
 type(esmf_field), public :: dt_cool_input_grid   !< Sub-layer cooling amount
 type(esmf_field), public :: ifd_input_grid   !< Model mode index. 0-diurnal model not
                                                     !< started; 1-diurnal model started.
 type(esmf_field), public :: qrain_input_grid   !< Sensible heat flux due to rainfall
 type(esmf_field), public :: tref_input_grid  !< Reference temperature
 type(esmf_field), public :: w_d_input_grid   !< Coefficient 4 to calculate d(tz)/d(ts)
 type(esmf_field), public :: w_0_input_grid   !< Coefficient 3 to calculate d(tz)/d(ts)
 type(esmf_field), public :: xs_input_grid   !< Salinity content in diurnal thermocline layer
 type(esmf_field), public :: xt_input_grid   !< Heat content in diurnal thermocline layer
 type(esmf_field), public :: xu_input_grid   !< u-current content in diurnal thermocline layer
 type(esmf_field), public :: xv_input_grid   !< v-current content in diurnal thermocline layer
 type(esmf_field), public :: xz_input_grid   !< Diurnal thermocline layer thickness
 type(esmf_field), public :: xtts_input_grid   !< d(xt)/d(ts)
 type(esmf_field), public :: xzts_input_grid   !< d(xz)/d(ts)
 type(esmf_field), public :: z_c_input_grid   !< Sub-layer cooling thickness
 type(esmf_field), public :: zm_input_grid   !< Oceanic mixed layer depth

 public :: read_input_atm_data
 public :: cleanup_input_atm_data
 public :: read_input_sfc_data
 public :: cleanup_input_sfc_data
 public :: read_input_nst_data
 public :: cleanup_input_nst_data
 public :: check_soilt
 public :: check_cnwat
 public :: quicksort
 public :: convert_winds
 public :: init_sfc_esmf_fields
 public :: dint2p
 
 contains

!> Read input grid atmospheric data driver.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Driver to read input grid nst data.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Driver to read input grid surface data.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
 subroutine read_input_sfc_data(localpet)

 implicit none

 integer, intent(in)             :: localpet

 call init_sfc_esmf_fields()

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

!> Create atmospheric esmf fields.
!!
!! @author George Gayno NCEP/EMC   
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

 allocate(tracers_input_grid(num_tracers_input))

 do i = 1, num_tracers_input
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

!> Create surface input grid esmf fields
!!
!! @author George Gayno NCEP/EMC 
 subroutine init_sfc_esmf_fields
 
 implicit none

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
    
 

 if (.not. vgfrc_from_climo) then
   print*,"- CALL FieldCreate FOR INPUT VEGETATION GREENNESS."
   veg_greenness_input_grid = ESMF_FieldCreate(input_grid, &
                     typekind=ESMF_TYPEKIND_R8, &
                     staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)
 endif
 
 if (.not. minmax_vgfrc_from_climo) then
   print*,"- CALL FieldCreate FOR INPUT MIN VEGETATION GREENNESS."
   min_veg_greenness_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)
    
    print*,"- CALL FieldCreate FOR INPUT MAX VEGETATION GREENNESS."
   max_veg_greenness_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)
 endif
 
 if (.not. lai_from_climo) then
    print*,"- CALL FieldCreate FOR INPUT LEAF AREA INDEX."
   lai_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", rc)
 endif
 end subroutine init_sfc_esmf_fields

!> Read input atmospheric data from spectral gfs (old sigio format).
!! 
!! @note Format used prior to July 19, 2017.
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

 if (num_tracers_input /= sighead%ntrac) then
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

 do k = 1, num_tracers_input

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

!> Read input atmospheric data from spectral gfs (global gaussian in
!! nemsio format. Starting July 19, 2017).
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

 do n = 1, num_tracers_input

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

!> Read input grid atmospheric fv3 gaussian nemsio files.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

 do n = 1, num_tracers_input

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

!> Read input grid fv3 atmospheric data 'warm' restart files.
!!
!! @note Routine reads tiled files in parallel.  Tile 1 is read by 
!! localpet 0; tile 2 by localpet 1, etc.  The number of pets
!! must be equal to or greater than the number of tiled files.  
!! Logic only tested with global input data of six tiles.
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

 do i = 1, num_tracers_input

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

!> Read fv3 netcdf gaussian history file.  Each task reads a horizontal
!! slice.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
 subroutine read_input_atm_gaussian_netcdf_file(localpet)

 use mpi

 implicit none

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

 do n = 1, num_tracers_input

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

!> Read input grid fv3 atmospheric tiled history files in netcdf
!! format.
!!
!! @note Routine reads tiled files in parallel.  Tile 1 is read by 
!! localpet 0; tile 2 by localpet 1, etc.  The number of pets
!! must be equal to or greater than the number of tiled files.  
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
 subroutine read_input_atm_tiled_history_file(localpet)

 use mpi

 implicit none

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
 print*,'- WILL PROCESS ', num_tracers_input, ' TRACERS.'

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

 do n = 1, num_tracers_input

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
 
!> Read input grid atmospheric fv3gfs grib2 files.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
 subroutine read_input_atm_grib2_file(localpet)

 use mpi
 use grib_mod
 
 use grib2_util, only                   : rh2spfh, rh2spfh_gfs, convert_omega

 implicit none

 integer, intent(in)                   :: localpet
 
 integer, parameter                    :: ntrac_max=14
 integer, parameter                    :: max_levs=1000

 character(len=300)                    :: the_file
 character(len=20)                     :: vname, &
                                          trac_names_vmap(ntrac_max), &
                                          tmpstr, & 
                                          method, tracers_input_vmap(num_tracers_input), &
                                          tracers_default(ntrac_max)

 integer                               :: i, j, k, n
 integer                               :: ii,jj
 integer                               :: rc, clb(3), cub(3)
 integer                               :: vlev, iret,varnum, o3n, pdt_num
 integer                               :: intrp_ier, done_print
 integer                               :: trac_names_oct10(ntrac_max)
 integer                               :: tracers_input_oct10(num_tracers_input)
 integer                               :: trac_names_oct11(ntrac_max)
 integer                               :: tracers_input_oct11(num_tracers_input)
 integer                               :: lugb, lugi, jdisc, jpdt(200), jgdt(200), iscale
 integer                               :: jids(200), jpdtn, jgdtn, octet_23, octet_29
 integer                               :: count_spfh, count_rh, count_icmr, count_scliwc
 integer                               :: count_cice, count_rwmr, count_scllwc, count

 logical                               :: conv_omega=.false., &
                                          hasspfh=.true., &
                                          isnative=.false., &
                                          use_rh=.false. , unpack, &
                                          all_empty, is_missing

 real(esmf_kind_r8), allocatable       :: dum2d_1(:,:)
                                          

 real(esmf_kind_r8)                    :: rlevs_hold(max_levs)
 real(esmf_kind_r8), allocatable       :: rlevs(:)
 real(esmf_kind_r4), allocatable       :: dummy2d(:,:)
 real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:), dummy2d_8(:,:),&
                                          u_tmp_3d(:,:,:), v_tmp_3d(:,:,:)
 real(esmf_kind_r8), pointer           :: presptr(:,:,:), psptr(:,:),tptr(:,:,:), &
                                          qptr(:,:,:), wptr(:,:,:),  &
                                          uptr(:,:,:), vptr(:,:,:)
 real(esmf_kind_r4)                    :: value
 real(esmf_kind_r8), parameter         :: p0 = 100000.0
 real(esmf_kind_r8), allocatable       :: dummy3d_col_in(:),dummy3d_col_out(:)
 real(esmf_kind_r8), parameter         :: intrp_missing = -999.0 
 real(esmf_kind_r4), parameter         :: lev_no_tr_fill = 20000.0
 real(esmf_kind_r4), parameter         :: lev_no_o3_fill = 40000.0

 type(gribfield)                       :: gfld
 
 tracers(:) = "NULL"
 
 trac_names_oct10 = (/1,  1,  14,  1,  1,  1,  1, 6,  6,   1,  6,  13,  13, 2 /)
 trac_names_oct11 = (/0, 22, 192, 23, 24, 25, 32, 1, 29, 100, 28, 193, 192, 2 /)

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

 if (localpet == 0) then

   lugb=14
   lugi=0
   call baopenr(lugb,the_file,iret)
   if (iret /= 0) call error_handler("ERROR OPENING GRIB2 FILE.", iret)

   jdisc   = 0     ! Search for discipline - meteorological products
   j = 0           ! Search at beginning of file.
   jpdt    = -9999  ! Array of values in product definition template, set to wildcard
   jids    = -9999  ! Array of values in identification section, set to wildcard
   jgdt    = -9999  ! Array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! Search for any grid definition number.
   jpdtn   = -1     ! Search for any product definition template number.
   unpack  =.false.

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

!----------------------------------------------------------------------
! Read first record and check if this is NCEP GEFS data. 
! This will determine what product definition template number to
! search for (Section 4/Octets 8-9).
!
! Section 1/Octets 6-7 is '7' (NCEP)
! Section 1/Octets 8-9 is '2' (NCEP Ensemble products).
!----------------------------------------------------------------------
  
   if (iret == 0) then
     if (gfld%idsect(1) == 7 .and. gfld%idsect(2) == 2) then
       print*,'- THIS IS NCEP GEFS DATA.'
       pdt_num = 1 ! Search for product definition template number 1.
                   ! Individual ensember forecast.
     else
       pdt_num = 0 ! Search for product definition template number 0.
                   ! Analysis or forecast.
     endif
   else
     call error_handler("READING GRIB2 FILE", iret)
   endif

!----------------------------------------------------------------------
! First, check for the vertical coordinate. If temperture at the 10 hybrid
! level is found, hybrid coordinates are assumed. Otherwise, data is on
! isobaric levels.
!----------------------------------------------------------------------

   j = 0
   jpdtn   = pdt_num  ! Search for the specific product definition template number.
   jpdt(1) = 0      ! Sect4/oct 10 - Parameter category - temperature field
   jpdt(2) = 0      ! Sect4/oct 11 - Parameter number - temperature
   jpdt(10) = 105   ! Sect4/oct 23 - Type of level - hybrid
   jpdt(12) = 10    ! Sect4/octs 25/28 - Value of hybrid level
   unpack=.false.

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
    
   if (iret == 0) then
     print*,'- DATA IS ON HYBRID LEVELS.'
     octet_23 = 105 ! Section 4/Oct 23 - type of first fixed surface.
     octet_29 = 255 ! Section 4/Oct 29 - type of second fixed surface (N/A).
     isnative=.true.
   else
     print*,'- DATA IS ON ISOBARIC LEVELS.'
     octet_23 = 100 ! Section 4/Oct 23 - type of first fixed surface.
     octet_29 = 255 ! Section 4/Oct 29 - type of second fixed surface (N/A).
     isnative=.false.
   endif

! Now count the number of vertical levels by searching for u-wind.
! Store the value of each level.

   rlevs_hold = -999.9
   lev_input = 0
   iret = 0
   j = 0
   jpdtn = -1
   jpdt = -9999

   do
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret /= 0) exit

     if (gfld%discipline == 0) then ! Discipline - meteorological products
       if (gfld%ipdtnum == pdt_num) then  ! Product definition template number.
         if (gfld%ipdtmpl(1) == 2 .and. gfld%ipdtmpl(2) == 2) then  ! u-wind
                                                                    ! Sect4/octs 10 and 11.
           if (gfld%ipdtmpl(10) == octet_23 .and. gfld%ipdtmpl(13) == octet_29) then  
                                                                    ! Sect4 octs 23 and 29.
                                                                    ! Hybrid or isobaric.
             lev_input = lev_input + 1
             iscale = 10 ** gfld%ipdtmpl(11)
             rlevs_hold(lev_input) = float(gfld%ipdtmpl(12))/float(iscale)
           endif
         endif
       endif
     endif
    
     j = k
   enddo

 endif ! read file on task 0.

 call mpi_barrier(MPI_COMM_WORLD, iret)
 call MPI_BCAST(isnative,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iret)
 call MPI_BCAST(lev_input,1,MPI_INTEGER,0,MPI_COMM_WORLD,iret)
 call MPI_BCAST(pdt_num,1,MPI_INTEGER,0,MPI_COMM_WORLD,iret)
 call MPI_BCAST(rlevs_hold, max_levs, MPI_INTEGER,0,MPI_COMM_WORLD,iret)

 allocate(slevs(lev_input))
 allocate(rlevs(lev_input))
 allocate(dummy3d_col_in(lev_input))
 allocate(dummy3d_col_out(lev_input))

 levp1_input = lev_input + 1

! Jili Dong add sort to re-order isobaric levels.

 do i = 1, lev_input
   rlevs(i) = rlevs_hold(i)
 enddo

 call quicksort(rlevs,1,lev_input)

 do i = 1, lev_input
   if (isnative) then
     write(slevs(i), '(i6)') nint(rlevs(i))
     slevs(i) = trim(slevs(i)) // " hybrid"
   else
     write(slevs(i), '(f11.2)') rlevs(i)
     slevs(i) = trim(slevs(i)) // " Pa"
   endif
 enddo

 if(localpet == 0) then
   do i = 1,lev_input
     print*, "- LEVEL AFTER SORT = ",trim(slevs(i))
   enddo
 endif

! Check to see if specfic humidity exists at all the same levels as ugrd.

 if (localpet == 0) then
   
   jpdtn = pdt_num ! Product definition template number.
   jpdt = -9999
   jpdt(1) = 1  ! Sect4/oct 10 - Parameter category - moisture
   jpdt(2) = 0  ! Sect4/oct 11 - Parameter number - specific humidity
   jpdt(10) =  octet_23 ! Sect4/oct 23 - type of level.
   unpack=.false.

   count_spfh=0

   do vlev = 1, lev_input
     j = 0
     jpdt(12) = nint(rlevs(vlev))

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_spfh = count_spfh + 1
     endif
   enddo

   jpdt(1) = 1  ! Sec4/oct 10 - Parameter category - moisture
   jpdt(2) = 1  ! Sec4/oct 11 - Parameter number - rel humidity
   count_rh=0

   do vlev = 1, lev_input
     j = 0
     jpdt(12) = nint(rlevs(vlev))

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_rh = count_rh + 1
     endif
   enddo

   if (count_spfh /= lev_input) then
     use_rh = .true.
   endif

   if (count_spfh == 0 .or. use_rh) then
     if (count_rh == 0) then
       call error_handler("READING ATMOSPHERIC WATER VAPOR VARIABLE.", 2)
     endif
     hasspfh = .false.  ! Will read rh and convert to specific humidity.
     trac_names_oct10(1) = 1
     trac_names_oct11(1) = 1
     print*,"- FILE CONTAINS RH."
   else
     print*,"- FILE CONTAINS SPFH."
   endif

 endif

 call MPI_BARRIER(MPI_COMM_WORLD, rc)
 call MPI_BCAST(hasspfh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,rc)
 
! Search for and count the number of tracers in the file.

 if (localpet == 0) then

   jpdtn = pdt_num ! Product definition template number.
   jpdt = -9999
   jpdt(10) =  octet_23 ! Sect4/oct 23 - type of level.
   unpack=.false.

   count_icmr=0
   count_scliwc=0
   count_cice=0
   count_rwmr=0
   count_scllwc=0

   do vlev = 1, lev_input

     j = 0
     jpdt(1) = 1  ! Sect4/oct 10 - Parameter category - moisture
     jpdt(2) = 23 ! Sect4/oct 11 - Parameter number - ice water mixing ratio
     jpdt(12) = nint(rlevs(vlev))

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_icmr = count_icmr + 1
     endif

     j = 0
     jpdt(1) = 1  ! Sect4/oct 10 - Parameter category - moisture
     jpdt(2) = 84 ! Sect4/oct 11 - Parameter number - cloud ice water content.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
           unpack, k, gfld, iret)

     if (iret == 0) then
       count_scliwc = count_scliwc + 1
     endif

     j = 0
     jpdt(1) = 6  ! Sect4/oct 10 - Parameter category - clouds
     jpdt(2) = 0  ! Sect4/oct 11 - Parameter number - cloud ice
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_cice = count_cice + 1
     endif

     j = 0
     jpdt(1) = 1   ! Sect4/oct 10 - Parameter category - moisture
     jpdt(2) = 24  ! Sect4/oct 11 - Parameter number - rain mixing ratio
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_rwmr = count_rwmr + 1
     endif

     j = 0
     jpdt(1) = 1   ! Sect4/oct 10 - Parameter category - moisture
     jpdt(2) = 83  ! Sect4/oct 11 - Parameter number - specific cloud liquid
                                  ! water content.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret == 0) then
       count_scllwc = count_scllwc + 1
     endif

   enddo

   if (count_icmr == 0) then
     if (count_scliwc == 0) then
       if (count_cice == 0) then
         print*,'- FILE DOES NOT CONTAIN CICE.'
       else
         trac_names_oct10(4) = 6 ! Sect4/oct 10 - Parameter category - clouds
         trac_names_oct11(4) = 0 ! Sect4/oct 11 - Parameter number - cloud ice
         print*,"- FILE CONTAINS CICE."
       endif
     else
       trac_names_oct10(4) = 1  ! Sect4/oct 10 - Parameter category - moisture
       trac_names_oct11(4) = 84 ! Sect4/oct 11 - Parameter number - cloud ice water content.
       print*,"- FILE CONTAINS SCLIWC."
     endif
   else
     print*,"- FILE CONTAINS ICMR."
   endif ! count of icmr

   if (count_rwmr == 0) then
     if (count_scllwc == 0) then
       print*,"- FILE DOES NOT CONTAIN SCLLWC."
     else
       trac_names_oct10(4) = 1  ! Sect4/oct 10 - Parameter category - moisture
       trac_names_oct11(4) = 83 ! Sect4/oct 11 - Parameter number - specific cloud liquid
                                ! water content.
       print*,"- FILE CONTAINS SCLLWC."
     endif
   else
     print*,"- FILE CONTAINS CLWMR."
   endif

 endif ! count of tracers/localpet = 0
   
 call MPI_BARRIER(MPI_COMM_WORLD, rc)
 call MPI_BCAST(trac_names_oct10,ntrac_max,MPI_INTEGER,0,MPI_COMM_WORLD,rc)
 call MPI_BCAST(trac_names_oct11,ntrac_max,MPI_INTEGER,0,MPI_COMM_WORLD,rc)
 
 print*,"- COUNT NUMBER OF TRACERS TO BE READ IN BASED ON PHYSICS SUITE TABLE"
 do n = 1, num_tracers_input

   vname = tracers_input(n)

   i = maxloc(merge(1.,0.,trac_names_vmap == vname),dim=1)

   tracers_input_vmap(n)=trac_names_vmap(i)
   tracers(n)=tracers_default(i)
   if(trim(tracers(n)) .eq. "o3mr") o3n = n

   tracers_input_oct10(n) = trac_names_oct10(i)
   tracers_input_oct11(n) = trac_names_oct11(i)

 enddo

!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------

 call init_atm_esmf_fields

 if (localpet == 0) then
   allocate(dummy2d(i_input,j_input))
   allocate(dummy2d_8(i_input,j_input))
   allocate(dummy3d(i_input,j_input,lev_input))
   allocate(dum2d_1(i_input,j_input))
 else
   allocate(dummy2d(0,0))
   allocate(dummy2d_8(0,0))
   allocate(dummy3d(0,0,0))
   allocate(dum2d_1(0,0))
 endif

!----------------------------------------------------------------------------------
! This program expects field levels from bottom to top. Fields in non-native 
! files read in from top to bottom. We will flip indices later. Fields on 
! native vertical coordinates read from bottom to top so those need no adjustments.
!----------------------------------------------------------------------------------
 
 if (localpet == 0) then

   print*,"- READ TEMPERATURE."

   jdisc   = 0     ! search for discipline - meteorological products
   j = 0           ! search at beginning of file.
   jpdt    = -9999  ! array of values in product definition template, set to wildcard
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! search for any grid definition number.
   jpdtn   =  pdt_num  ! Search for specific product definition template number.
   jpdt(1) = 0      ! Sect 4/oct 10 - parameter category - temperature
   jpdt(2) = 0      ! Sect 4/oct 11 - parameter number - temperature
   jpdt(10) = octet_23 ! Sect4/oct 23 - type of level.

   unpack=.true.

   do vlev = 1, lev_input

     jpdt(12) = nint(rlevs(vlev))

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
     if (iret /= 0) then 
       call error_handler("READING IN TEMPERATURE AT LEVEL "//trim(slevs(vlev)),iret)
     endif

     dum2d_1 = reshape(gfld%fld, (/i_input,j_input/) )

     dummy3d(:,:,vlev) = dum2d_1

   enddo

 endif ! Read of temperature

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID TEMPERATURE."
 call ESMF_FieldScatter(temp_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

! Read tracers

 do n = 1, num_tracers_input

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

     jdisc   = 0     ! search for discipline - meteorological products
     jpdt    = -9999  ! array of values in product definition template, set to wildcard
     jids    = -9999  ! array of values in identification section, set to wildcard
     jgdt    = -9999  ! array of values in grid definition template, set to wildcard
     jgdtn   = -1     ! search for any grid definition number.
     jpdtn   =  pdt_num  ! Search for the product definition template number.
     jpdt(10) = octet_23 ! Sect4/oct 23 - type of level.
     unpack = .false.

     count = 0

     do vlev = 1, lev_input

       j = 0
       jpdt(1) = tracers_input_oct10(n)
       jpdt(2) = tracers_input_oct11(n)
       jpdt(12) = nint(rlevs(vlev))

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

       if (iret == 0) then
         count = count + 1
       endif

     enddo
     iret=count

     ! Check to see if file has any data for this tracer
     if (iret == 0) then
       all_empty = .true.
     else
       all_empty = .false.
     endif
 
     is_missing = .false.

     do vlev = 1, lev_input

       unpack=.true.
       j = 0
       jpdt(1) = tracers_input_oct10(n)
       jpdt(2) = tracers_input_oct11(n)
       jpdt(12) = nint(rlevs(vlev) )

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

       if (iret == 0) then ! found data
         dummy2d = real((reshape(gfld%fld, (/i_input,j_input/) )), kind=esmf_kind_r4)
       else ! did not find data.
        if (trim(method) .eq. 'intrp' .and. .not.all_empty) then
          dummy2d = intrp_missing 
          is_missing = .true.
        else
          ! Abort if input data has some data for current tracer, but has
          ! missing data below 200 mb/ above 400mb
            if (.not.all_empty .and. n == o3n) then
              if (rlevs(vlev) .lt. lev_no_o3_fill) &
                call error_handler("TRACER "//trim(tracers(n))//" HAS MISSING DATA AT "//trim(slevs(vlev))//&
                  ". SET MISSING VARIABLE CONDITION TO 'INTRP' TO AVOID THIS ERROR", 1)
            elseif (.not.all_empty .and. n .ne. o3n) then 
              if (rlevs(vlev) .gt. lev_no_tr_fill) &
                call error_handler("TRACER "//trim(tracers(n))//" HAS MISSING DATA AT "//trim(slevs(vlev))//&
                  ". SET MISSING VARIABLE CONDITION TO 'INTRP' TO AVOID THIS ERROR.", 1)
            endif 
          ! If entire array is empty and method is set to intrp, switch method to fill
          if (trim(method) .eq. 'intrp' .and. all_empty) method='set_to_fill' 

          call handle_grib_error(vname, slevs(vlev),method,value,varnum,iret,var=dummy2d)
          if (iret==1) then ! missing_var_method == skip or no entry
            if ( (tracers_input_oct10(n) == 1 .and. tracers_input_oct11(n) == 0) .or. &  ! spec humidity
                 (tracers_input_oct10(n) == 1 .and. tracers_input_oct11(n) == 1) .or. &  ! rel humidity
                 (tracers_input_oct10(n) == 14 .and. tracers_input_oct11(n) == 192) ) then ! ozone
              call error_handler("READING IN "//trim(tracers(n))//" AT LEVEL "//trim(slevs(vlev))&
                        //". SET A FILL VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
            endif
          endif
        endif ! method intrp
      endif !iret<=0

      if (n==1 .and. .not. hasspfh) then 
        if (trim(external_model) .eq. 'GFS') then
          print *,'- CALL CALRH GFS'
          call rh2spfh_gfs(dummy2d,rlevs(vlev),dummy3d(:,:,vlev))
        else 
          print *,'- CALL CALRH non-GFS'
          call rh2spfh(dummy2d,rlevs(vlev),dummy3d(:,:,vlev))
        end if
      endif

       dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)

     enddo !vlev

! Jili Dong interpolation for missing levels 
     if (is_missing .and. trim(method) .eq. 'intrp') then
       print *,'- INTERPOLATE TRACER '//trim(tracers(n))
       done_print = 0
       do jj = 1, j_input
         do ii = 1, i_input
           dummy3d_col_in=dummy3d(ii,jj,:)
           call dint2p(rlevs,dummy3d_col_in,lev_input,rlevs,dummy3d_col_out,    &
                        lev_input, 2, intrp_missing, intrp_ier) 
           if (intrp_ier .gt. 0) call error_handler("Interpolation failed.",intrp_ier)
           dummy3d(ii,jj,:)=dummy3d_col_out
         enddo
       enddo
       do vlev=1,lev_input
         dummy2d = real(dummy3d(:,:,n) , kind=esmf_kind_r4)
         if (any(dummy2d .eq. intrp_missing)) then 
           ! If we're outside the appropriate region, don't fill but error instead
           if (n == o3n .and. rlevs(vlev) .lt. lev_no_o3_fill) then
             call error_handler("TRACER "//trim(tracers(n))//" HAS MISSING DATA AT "//trim(slevs(vlev)),1)
           elseif (n .ne. o3n .and. rlevs(vlev) .gt. lev_no_tr_fill) then
             call error_handler("TRACER "//trim(tracers(n))//" HAS MISSING DATA AT "//trim(slevs(vlev)),1)
           else ! we're okay to fill missing data with provided fill value
             if (done_print .eq. 0) then
               print*, "Pressure out of range of existing data. Defaulting to fill value."
               done_print = 1
             end if !done print
             where(dummy2d .eq. intrp_missing) dummy2d = value
             dummy3d(:,:,vlev) = dummy2d
           end if !n & lev
         endif ! intrp_missing
         ! zero out negative tracers from interpolation/extrapolation
         where(dummy3d(:,:,vlev) .lt. 0.0)  dummy3d(:,:,vlev) = 0.0
!        print*,'tracer af intrp',vlev, maxval(dummy3d(:,:,vlev)),minval(dummy3d(:,:,vlev))
       end do !nlevs do
     end if !if intrp
   endif !localpet == 0

   if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT ", trim(tracers_input_vmap(n))
   call ESMF_FieldScatter(tracers_input_grid(n), dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 enddo
 
 deallocate(dummy3d_col_in, dummy3d_col_out)
 
 call read_winds(u_tmp_3d,v_tmp_3d,localpet,octet_23,rlevs,lugb,pdt_num)

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT U-WIND."
 call ESMF_FieldScatter(u_input_grid, u_tmp_3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT V-WIND."
 call ESMF_FieldScatter(v_input_grid, v_tmp_3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then

   print*,"- READ SURFACE PRESSURE."
   jdisc   = 0     ! search for discipline - meteorological products
   j = 0           ! search at beginning of file.
   jpdt    = -9999  ! array of values in product definition template, set to wildcard
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! search for any grid definition number.
   jpdtn   =  pdt_num  ! Search for the product definition template number.
   jpdt(1) = 3      ! Sect4/oct 10 - param category - mass
   jpdt(2) = 0      ! Sect4/oct 11 - param number - pressure
   jpdt(10) = 1     ! Sect4/oct 23 - type of level - ground surface
   unpack=.true.

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
   if (iret /= 0) call error_handler("READING SURFACE PRESSURE RECORD.", iret)

   dummy2d_8 = reshape(gfld%fld, (/i_input,j_input/) )

 endif ! Read surface pressure

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID SURFACE PRESSURE."
 call ESMF_FieldScatter(ps_input_grid, dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

! Read dzdt.

 if (localpet == 0) then

   print*,"- READ DZDT."
   vname = "dzdt"
   call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
                         loc=varnum)

   jdisc   = 0     ! search for discipline - meteorological products
   j = 0           ! search at beginning of file.
   jpdt    = -9999  ! array of values in product definition template, set to wildcard
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! search for any grid definition number.
   jpdtn   =  pdt_num ! Search for the product definition template number.
   jpdt(1) = 2      ! Sect4/oct 10 - param category - momentum
   jpdt(2) = 9      ! Sect4/oct 11 - param number - dzdt
   jpdt(10) = octet_23 ! Sect4/oct 23 - type of level

   unpack=.true.

   do vlev = 1, lev_input

     jpdt(12) = nint(rlevs(vlev))

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret /= 0) then ! dzdt not found, look for omega.
       print*,"DZDT not available at level ", trim(slevs(vlev)), " so checking for VVEL"
       jpdt(2) = 8  ! Sect4/oct 11 - parameter number - omega
       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
       if (iret /= 0) then
         call handle_grib_error(vname, slevs(vlev),method,value,varnum,iret,var8=dum2d_1)
         if (iret==1) then ! missing_var_method == skip 
           cycle
         endif
       else
         conv_omega = .true.
         dum2d_1 = reshape(gfld%fld, (/i_input,j_input/) )
       endif
     else ! found dzdt
       dum2d_1 = reshape(gfld%fld, (/i_input,j_input/) )
     endif

     dummy3d(:,:,vlev) = dum2d_1

   enddo

 endif ! Read of dzdt

 call mpi_bcast(conv_omega,1,MPI_LOGICAL,0,MPI_COMM_WORLD,rc)

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT DZDT."
 call ESMF_FieldScatter(dzdt_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

! Read terrain

 if (localpet == 0) then

   print*,"- READ TERRAIN."
   jdisc   = 0     ! search for discipline - meteorological products
   j = 0           ! search at beginning of file.
   jpdt    = -9999  ! array of values in product definition template, set to wildcard
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! search for any grid definition number.
   jpdtn   =  pdt_num  ! Search for the product definition template number.
   jpdt(1) = 3      ! Sect4/oct 10 - param category - mass
   jpdt(2) = 5      ! Sect4/oct 11 - param number - geopotential height
   jpdt(10) = 1     ! Sect4/oct 23 - type of level - ground surface
   unpack=.true.

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
   if (iret /= 0) call error_handler("READING TERRAIN HEIGHT RECORD.", iret)

   dummy2d_8 = reshape(gfld%fld, (/i_input,j_input/) )

 endif ! read of terrain.

 if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID TERRAIN."
 call ESMF_FieldScatter(terrain_input_grid, dummy2d_8, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 deallocate(dummy2d, dummy2d_8)
 
if (.not. isnative) then

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
  do n=1,num_tracers_input
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
 
else ! is native coordinate (hybrid).

! For native files, read in pressure field directly from file but don't flip levels

   if (localpet == 0) then

    print*,"- READ PRESSURE."

    jdisc   = 0     ! search for discipline - meteorological products
    j = 0           ! search at beginning of file.
    jpdt    = -9999  ! array of values in product definition template, set to wildcard
    jids    = -9999  ! array of values in identification section, set to wildcard
    jgdt    = -9999  ! array of values in grid definition template, set to wildcard
    jgdtn   = -1     ! search for any grid definition number.
    jpdtn   =  pdt_num ! Search for the product definition template number.
    jpdt(1) = 3      ! Sect4/oct 10 - parameter category - mass
    jpdt(2) = 0      ! Sect4/oct 11 - parameter number - pressure
    jpdt(10) = octet_23 ! Sect4/oct 23 - type of level.
    unpack=.true.

    do vlev = 1, lev_input

      jpdt(12) = nint(rlevs(vlev))
      call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
      if (iret /= 0) then
        call error_handler("READING IN PRESSURE AT LEVEL "//trim(slevs(vlev)),iret)
      endif

      dum2d_1 = reshape(gfld%fld, (/i_input,j_input/) )

      dummy3d(:,:,vlev) = dum2d_1

    enddo

  endif  ! localpet == 0

  if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID PRESSURE."
  call ESMF_FieldScatter(pres_input_grid, dummy3d, rootpet=0, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldScatter", rc)

 endif

 deallocate(dummy3d, dum2d_1) 
 
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
 
 if (localpet == 0) call baclose(lugb, rc)

 end subroutine read_input_atm_grib2_file

!> Read input grid surface data from a spectral gfs gaussian sfcio
!! file.
!!
!! @note Prior to July 19, 2017.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid surface data from a spectral gfs gaussian nemsio
!! file.
!!
!! @note Format used by gfs starting July 19, 2017.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid surface data from an fv3 gaussian nemsio file.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid surface data from fv3 tiled warm 'restart' files.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid surface data from tiled 'history' files (netcdf) or 
!! gaussian netcdf files.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid surface data from a grib2 file.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames 
 subroutine read_input_sfc_grib2_file(localpet)

   use mpi
   use grib_mod
   use program_setup, only : vgtyp_from_climo, sotyp_from_climo
   use model_grid, only    : input_grid_type
   use search_util

   implicit none

   integer, intent(in)                   :: localpet

   character(len=250)                    :: the_file
   character(len=250)                    :: geo_file
   character(len=200)                    :: err_msg
   character(len=20)                     :: vname, vname_file, slev
   character(len=50)                     :: method
   character(len=20)                     :: to_upper
 
   integer                               :: rc, varnum, iret, i, j,k
   integer                               :: ncid2d, varid, varsize
   integer                               :: lugb, lugi
   integer                               :: jdisc, jgdtn, jpdtn, pdt_num
   integer                               :: jids(200), jgdt(200), jpdt(200)

   logical                               :: rap_latlon, unpack

   real(esmf_kind_r4)                    :: value
   real(esmf_kind_r4), allocatable       :: dummy2d(:,:)
   real(esmf_kind_r8), allocatable       :: icec_save(:,:)
   real(esmf_kind_r4), allocatable       :: dummy1d(:)
   real(esmf_kind_r8), allocatable       :: dummy2d_8(:,:),dummy2d_82(:,:),tsk_save(:,:)
   real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:), dummy3d_stype(:,:,:)
   integer(esmf_kind_i4), allocatable    :: slmsk_save(:,:)
   integer(esmf_kind_i8), allocatable    :: dummy2d_i(:,:)
   
   type(gribfield)                       :: gfld
    
   rap_latlon = trim(to_upper(external_model))=="RAP" .and. trim(input_grid_type) == "rotated_latlon"

   the_file = trim(data_dir_input_grid) // "/" // trim(grib2_file_input_grid)
   geo_file = trim(geogrid_file_input_grid)
   
   print*,"- READ SFC DATA FROM GRIB2 FILE: ", trim(the_file)

! Determine the number of soil layers in file.

   if (localpet == 0) then

     lugb=12
     call baopenr(lugb,the_file,rc)
     if (rc /= 0) call error_handler("ERROR OPENING GRIB2 FILE.", rc)

     j       = 0      ! search at beginning of file
     lugi    = 0      ! no grib index file
     jdisc   = -1     ! search for any discipline
     jpdtn   = -1     ! search for any product definition template number
     jgdtn   = -1     ! search for any grid definition template number
     jids    = -9999  ! array of values in identification section, set to wildcard
     jgdt    = -9999  ! array of values in grid definition template, set to wildcard
     jpdt    = -9999  ! array of values in product definition template, set to wildcard
     unpack  = .false. ! unpack data

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc == 0) then
       if (gfld%idsect(1) == 7 .and. gfld%idsect(2) == 2) then
         print*,'- THIS IS NCEP GEFS DATA.'
         pdt_num = 1
       else
         pdt_num = 0
       endif
     else
       if (rc /= 0) call error_handler("ERROR READING GRIB2 FILE.", rc)
     endif

     j = 0
     lsoil_input = 0

     do

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

       if (rc /= 0) exit

       if (gfld%discipline == 2) then ! discipline - land products
         if (gfld%ipdtnum == pdt_num) then  ! prod template number - analysis or forecast at single level.
           if (gfld%ipdtmpl(1) == 0 .and. gfld%ipdtmpl(2) == 2) then  ! soil temp
                                                                      ! Sect4/octs 10 and 11
             if (gfld%ipdtmpl(10) == 106 .and. gfld%ipdtmpl(13) == 106) then  ! Sect4/octs 23/29.
                                                                              ! Layer below ground.
               lsoil_input = lsoil_input + 1
             endif
           endif
         endif
       endif
    
       j = k

     enddo

     print*, "- FILE HAS ", lsoil_input, " SOIL LEVELS."
     if (lsoil_input == 0) call error_handler("COUNTING SOIL LEVELS.", rc)

   endif ! localpet == 0

   call MPI_BARRIER(MPI_COMM_WORLD, rc)
   call MPI_BCAST(lsoil_input,1,MPI_INTEGER,0,MPI_COMM_WORLD,rc)

 ! We need to recreate the soil fields if we have something other than 4 levels

   if (lsoil_input /= 4) then
   
     call ESMF_FieldDestroy(soil_temp_input_grid, rc=rc)
     call ESMF_FieldDestroy(soilm_tot_input_grid, rc=rc)
     call ESMF_FieldDestroy(soilm_liq_input_grid, rc=rc)
     
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
   
   endif
 
   if (localpet == 0) then
     allocate(dummy2d(i_input,j_input))
     allocate(slmsk_save(i_input,j_input))
     allocate(tsk_save(i_input,j_input))
     allocate(icec_save(i_input,j_input))
     allocate(dummy2d_8(i_input,j_input))
     allocate(dummy2d_82(i_input,j_input))
     allocate(dummy3d(i_input,j_input,lsoil_input))
   else
     allocate(dummy3d(0,0,0))
     allocate(dummy2d_8(0,0))
     allocate(dummy2d_82(0,0))
     allocate(dummy2d(0,0))
     allocate(slmsk_save(0,0))
   endif
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! These variables are always in grib files, or are required, so no need to check for them 
 ! in the varmap table. If they can't be found in the input file, then stop the program.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (localpet == 0) then

     print*,"- READ TERRAIN."
 
     j = 0
     jdisc   = 0  ! Search for discipline 0 - meteorological products
     jpdt    = -9999  ! array of values in product definition template, set to wildcard.
     jpdtn   = pdt_num       ! search for product definition template number 0 - anl or fcst.
     jpdt(1) = 3  ! Sec4/oct 10 - param cat - mass field
     jpdt(2) = 5  ! Sec4/oct 11 - param number - geopotential height
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, rc)
     if (rc /= 0) call error_handler("READING TERRAIN.", rc)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
!    print*,'orog ', maxval(dummy2d_8),minval(dummy2d_8)
   
   endif

   print*,"- CALL FieldScatter FOR INPUT TERRAIN."
   call ESMF_FieldScatter(terrain_input_grid, dummy2d_8, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
    
   if (localpet == 0) then

     print*,"- READ SEAICE FRACTION."
 
     jdisc   = 10 ! Search for discipline - ocean products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for product def template number 0 - anl or fcst.
     jpdt    = -9999  ! Array of values in Sec 4 product definition template;
                      ! Initialize to wildcard.
     jpdt(1) = 2  ! Sec4/oct 10 - parameter category - ice
     jpdt(2) = 0  ! Sec4/oct 11 - parameter number - ice cover
     unpack=.true.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, rc)
     if (rc /= 0) call error_handler("READING SEAICE FRACTION.", rc)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
!    print*,'icec ', maxval(dummy2d_8),minval(dummy2d_8)
 
     icec_save = dummy2d_8

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID SEAICE FRACTION."
   call ESMF_FieldScatter(seaice_fract_input_grid, dummy2d_8 ,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

!----------------------------------------------------------------------------------
! GFS v14 and v15.2 grib data has two land masks.  LANDN is created by
! nearest neighbor interpolation.  LAND is created by bilinear interpolation.
! LANDN matches the bitmap.  So use it first.  For other GFS versions or other models,
! use LAND. Mask in grib file is '1' (land), '0' (not land).  Add sea/lake ice category
! '2' based on ice concentration.
!----------------------------------------------------------------------------------

   if (localpet == 0) then

     print*,"- READ LANDSEA MASK."

     jdisc   = 2  ! Search for discipline - land products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for product definition template number 0 - anl or fcst.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec 4.
     jpdt(1) = 0      ! Sec4/oct 10 - parameter category - veg/biomass
     jpdt(2) = 218    ! Sec4/oct 11 - parameter number - land nearest neighbor
     unpack=.true.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc == 0) then

       print*,'landnn ', maxval(gfld%fld),minval(gfld%fld)

     else

       jdisc   = 2  ! Search for discipline - land products
       j = 0        ! Search at beginning of file.
       jpdtn   = pdt_num  ! Search for product def template number 0 - anl or fcst.
       jpdt    = -9999  ! Initialize array of values in product definition template - Sec 4.
       jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
       jpdt(2) = 0  ! Sec4/oct 11 - parameter number - land cover (fraction)
       unpack=.true.
       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
              unpack, k, gfld, rc)
       if (rc /= 0) call error_handler("READING LANDSEA MASK.", rc)
    
!      print*,'land ', maxval(gfld%fld),minval(gfld%fld)

     endif

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
 
     do j = 1, j_input
       do i = 1, i_input
         if(dummy2d_8(i,j) < 0.5_esmf_kind_r8) dummy2d_8(i,j)=0.0
         if(icec_save(i,j) > 0.15_esmf_kind_r8) then 
           dummy2d_8(i,j) = 2.0_esmf_kind_r8
         endif
       enddo
     enddo

     slmsk_save = nint(dummy2d_8)

     deallocate(icec_save)

   endif ! read land mask

   print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
   call ESMF_FieldScatter(landsea_mask_input_grid, dummy2d_8 ,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   if (localpet == 0) then

     print*,"- READ SEAICE SKIN TEMPERATURE."

     jdisc   = 0  ! Search for discipline - meteorological products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for product definition template number 0 - anl or fcst.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 0  ! Sec4/oct 10 - parameter category - temperature
     jpdt(2) = 0  ! Sec4/oct 11 - parameter number - temperature
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.
     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     if (rc /= 0) call error_handler("READING SEAICE SKIN TEMP.", rc)

!    print*,'ti ',maxval(gfld%fld),minval(gfld%fld)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID SEAICE SKIN TEMPERATURE."
   call ESMF_FieldScatter(seaice_skin_temp_input_grid, dummy2d_8 ,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
     call error_handler("IN FieldScatter", rc)

!----------------------------------------------------------------------------------
! Read snow fields.  Zero out at non-land points and undefined points (points
! removed using the bitmap).  Program expects depth and liquid equivalent
! in mm.
!----------------------------------------------------------------------------------

   if (localpet == 0) then

     print*,"- READ SNOW LIQUID EQUIVALENT."

     jdisc   = 0 ! Search for discipline - meteorological products
     j = 0       ! Search at beginning of file.
     jpdtn   = pdt_num ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 1  ! Sec4/oct 10 - parameter category - moisture
     jpdt(2) = 13 ! Sec4/oct 11 - parameter number - liquid equiv snow depth
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     if (rc /= 0) call error_handler("READING SNOW LIQUID EQUIVALENT.", rc)

!    print*,'weasd ', maxval(gfld%fld),minval(gfld%fld)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

     do j = 1, j_input
       do i = 1, i_input
         if(slmsk_save(i,j) == 0) dummy2d_8(i,j) = 0.0
       enddo
     enddo

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID SNOW LIQUID EQUIVALENT."
   call ESMF_FieldScatter(snow_liq_equiv_input_grid, dummy2d_8 ,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   if (localpet == 0) then

     print*,"- READ SNOW DEPTH."

     jdisc   = 0  ! Search for discipline - meteorological products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999 ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 1  ! Sec4/oct 10 - parameter category - moisture
     jpdt(2) = 11 ! Sec4/oct 11 - parameter number - snow depth
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0) then
       call error_handler("READING SNOW DEPTH.", rc)
     else
       gfld%fld = gfld%fld * 1000.0
!      print*,'snod ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
     endif

     do j = 1, j_input
     do i = 1, i_input
       if(slmsk_save(i,j) == 0) dummy2d_8(i,j) = 0.0
     enddo
     enddo

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID SNOW DEPTH."
   call ESMF_FieldScatter(snow_depth_input_grid,dummy2d_8,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
    
   if (localpet == 0) then

     print*,"- READ T2M."

     jdisc   = 0 ! Search for discipline - meteorological products
     j = 0       ! Search at beginning of file.
     jpdtn   = pdt_num ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 0    ! Sec4/oct 10 - parameter category - temperature
     jpdt(2) = 0    ! Sec4/oct 11 - parameter number - temperature
     jpdt(10) = 103 ! Sec4/oct 23 - type of level - height above ground surface
     jpdt(12) = 2   ! Sec4/octs 25-28 - 2 meters above ground.
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0) call error_handler("READING T2M.", rc)
!    print*,'t2m ', maxval(gfld%fld),minval(gfld%fld)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID T2M."
   call ESMF_FieldScatter(t2m_input_grid, dummy2d_8, rootpet=0,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   if (localpet == 0) then

     print*,"- READ Q2M."

     jdisc   = 0  ! Search for discipline - meteorological products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 1  ! Sec4/oct 10 - parameter category - moisture
     jpdt(2) = 0  ! Sec4/oct 11 - parameter number - specific humidity
     jpdt(10) = 103 ! Sec4/oct 23 - type of level - height above ground surface
     jpdt(12) = 2 ! Sec4/octs 25-28 - 2 meters above ground.
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     if (rc /=0) call error_handler("READING Q2M.", rc)

!    print*,'q2m ',maxval(gfld%fld),minval(gfld%fld)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID Q2M."
   call ESMF_FieldScatter(q2m_input_grid,dummy2d_8, rootpet=0,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
    
   if (localpet == 0) then

     print*,"- READ SKIN TEMPERATURE."

     jdisc   = 0  ! Search for discipline - meteorological products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 0  ! Sec4/oct 10 - parameter category - temperature
     jpdt(2) = 0  ! Sec4/oct 11 - parameter number - temperature
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0 ) call error_handler("READING SKIN TEMPERATURE.", rc)
!    print*,'skint ', maxval(gfld%fld),minval(gfld%fld)

     dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

     tsk_save(:,:) = dummy2d_8

     do j = 1, j_input
       do i = 1, i_input
        if(slmsk_save(i,j) == 0 .and. dummy2d_8(i,j) < 271.2) then
!         print*,'too cool SST ',i,j,dummy2d_8(i,j)
          dummy2d_8(i,j) = 271.2
        endif
        if(slmsk_save(i,j) == 0 .and. dummy2d_8(i,j) > 310.) then
!         print*,'too hot SST ',i,j,dummy2d_8(i,j)
          dummy2d_8(i,j) = 310.0
        endif
       enddo
     enddo

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID SKIN TEMPERATURE"
   call ESMF_FieldScatter(skin_temp_input_grid,dummy2d_8,rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
    
! srflag not in files. Set to zero.

   if (localpet == 0) dummy2d_8 = 0.0
 
   print*,"- CALL FieldScatter FOR INPUT GRID SRFLAG"
   call ESMF_FieldScatter(srflag_input_grid,dummy2d_8, rootpet=0,rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   if (localpet == 0) then

     print*,"- READ SOIL TYPE."

     jdisc   = 2  ! Search for discipline - land products
     j = 0        ! Search at beginning of file
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template - Sec4
     jpdt(1) = 3  ! Sec4/oct 10 - parameter category - soil products
     jpdt(2) = 0  ! Sec4/oct 11 - parameter number - soil type
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc == 0 ) then
!      print*,'soil type ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d = reshape(real(gfld%fld,kind=esmf_kind_r4) , (/i_input,j_input/))

     endif

     if (rc /= 0 .and. (trim(to_upper(external_model))=="HRRR" .or. rap_latlon) .and. geo_file .ne. "NULL")  then
     ! Some HRRR and RAP files don't have dominant soil type in the output, but the geogrid files
     ! do, so this gives users the option to provide the geogrid file and use input soil
     ! type 
       print*, "OPEN GEOGRID FILE ", trim(geo_file)
       rc = nf90_open(geo_file,NF90_NOWRITE,ncid2d)
       call netcdf_err(rc,"READING GEOGRID FILE")

       print*, "INQURE ABOUT DIM IDS"
       rc = nf90_inq_dimid(ncid2d,"west_east",varid)
       call netcdf_err(rc,"READING west_east DIMENSION FROM GEOGRID FILE")
     
       rc = nf90_inquire_dimension(ncid2d,varid,len=varsize)
       call netcdf_err(rc,"READING west_east DIMENSION SIZE")
       if (varsize .ne. i_input) call error_handler ("GEOGRID FILE GRID SIZE DIFFERS FROM INPUT DATA.", -1)
        
       print*, "INQUIRE ABOUT SOIL TYPE FROM GEOGRID FILE"
       rc = nf90_inq_varid(ncid2d,"SCT_DOM",varid)
       call netcdf_err(rc,"FINDING SCT_DOM IN GEOGRID FILE")
     
       print*, "READ SOIL TYPE FROM GEOGRID FILE "
       rc = nf90_get_var(ncid2d,varid,dummy2d)
       call netcdf_err(rc,"READING SCT_DOM FROM FILE")
       
       print*, "INQUIRE ABOUT SOIL TYPE FRACTIONS FROM GEOGRID FILE"
       rc = nf90_inq_varid(ncid2d,"SOILCTOP",varid)
       call netcdf_err(rc,"FINDING SOILCTOP IN GEOGRID FILE")
     
       allocate(dummy3d_stype(i_input,j_input,16))
       print*, "READ SOIL TYPE FRACTIONS FROM GEOGRID FILE "
       rc = nf90_get_var(ncid2d,varid,dummy3d_stype)
       call netcdf_err(rc,"READING SCT_DOM FROM FILE")

       print*, "CLOSE GEOGRID FILE "
       iret = nf90_close(ncid2d)
     
     ! There's an issue with the geogrid file containing soil type water at land points. 
     ! This correction replaces the soil type at these points with the soil type with
     ! the next highest fractional coverage.
       allocate(dummy1d(16))
       do j = 1, j_input
       do i = 1, i_input
         if(dummy2d(i,j) == 14.0_esmf_kind_r4 .and. slmsk_save(i,j) == 1) then
           dummy1d(:) = real(dummy3d_stype(i,j,:),kind=esmf_kind_r4)
           dummy1d(14) = 0.0_esmf_kind_r4
           dummy2d(i,j) = real(MAXLOC(dummy1d, 1),esmf_kind_r4)
         endif
       enddo
       enddo
       deallocate(dummy1d)
       deallocate(dummy3d_stype)
     endif ! failed
   
     if ((rc /= 0 .and. trim(to_upper(external_model)) /= "HRRR" .and. .not. rap_latlon) & 
       .or. (rc /= 0 .and. (trim(to_upper(external_model)) == "HRRR" .or. rap_latlon))) then
       if (.not. sotyp_from_climo) then
         call error_handler("COULD NOT FIND SOIL TYPE IN FILE. PLEASE SET SOTYP_FROM_CLIMO=.TRUE. . EXITING", rc)
       else
         vname = "sotyp"
         slev = "surface"
         call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
                             loc=varnum)  
         call handle_grib_error(vname, slev ,method,value,varnum,rc, var= dummy2d)
         if (rc == 1) then ! missing_var_method == skip or no entry in varmap table
            print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. WILL NOT "//&
                       "SCALE SOIL MOISTURE FOR DIFFERENCES IN SOIL TYPE. "
            dummy2d(:,:) = -99999.0_esmf_kind_r4
         endif
       endif
     endif
   
   ! In the event that the soil type on the input grid still contains mismatches between 
   ! soil type and landmask, this correction is a last-ditch effort to replace these points
   ! with soil type from a nearby land point.

     if (.not. sotyp_from_climo) then
       do j = 1, j_input
       do i = 1, i_input
         if(dummy2d(i,j) == 14.0_esmf_kind_r4 .and. slmsk_save(i,j) == 1) dummy2d(i,j) = -99999.9_esmf_kind_r4
       enddo
       enddo

       allocate(dummy2d_i(i_input,j_input))
       dummy2d_8 = real(dummy2d,esmf_kind_r8)
       dummy2d_i(:,:) = 0
       where(slmsk_save == 1) dummy2d_i = 1
   
       call search(dummy2d_8,dummy2d_i,i_input,j_input,1,230)
       deallocate(dummy2d_i)
     else
       dummy2d_8=real(dummy2d,esmf_kind_r8)
     endif
   
     print*,'sotype ',maxval(dummy2d_8),minval(dummy2d_8)

   endif ! read of soil type
  
   print*,"- CALL FieldScatter FOR INPUT GRID SOIL TYPE."
   call ESMF_FieldScatter(soil_type_input_grid,dummy2d_8, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   deallocate(dummy2d)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
 ! Begin variables whose presence in grib2 files varies, but no climatological
 ! data is available, so we have to account for values in the varmap table
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   if (.not. vgfrc_from_climo) then  

     if (localpet == 0) then

       print*,"- READ VEG FRACTION."

       jdisc   = 2  ! Search for discipline - land products
       j = 0        ! Search at beginning of file.
       jpdtn   = pdt_num  ! Search for the product definition template number.
       jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
       jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
       jpdt(2) = 4  ! Sec4/oct 11 - parameter number - vegetation
       jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
       unpack=.true.

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, rc)

       if (rc /= 0 )then
         err_msg="COULD NOT FIND VEGETATION FRACTION IN FILE. PLEASE SET VGFRC_FROM_CLIMO=.TRUE."
         call error_handler(err_msg, rc)
       else
         if (maxval(gfld%fld) > 2.0) gfld%fld = gfld%fld / 100.0
!        print*,'vfrac ', maxval(gfld%fld),minval(gfld%fld)
         dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

       endif

     endif ! localpet 0

     print*,"- CALL FieldScatter FOR INPUT GRID VEG GREENNESS."
     call ESMF_FieldScatter(veg_greenness_input_grid,dummy2d_8, rootpet=0, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter", rc)

   endif

   if (.not. minmax_vgfrc_from_climo) then

     if (localpet == 0) then

       print*,"- READ MIN VEG FRACTION."

       jdisc   = 2  ! Search for discipline - land products
       j = 1105 ! grib2 file does not distinguish between the various veg
                ! fractions. Need to search using record number.
       jpdtn   = pdt_num  ! Search for the product definition template number.
       jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
       jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
       jpdt(2) = 4  ! Sec4/oct 11 - parameter number - vegetation
       jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
       unpack=.true.

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, rc)

       if (rc /= 0) then
         j = 1101 ! Have to search by record number.
         call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
                unpack, k, gfld, rc)
         if (rc /= 0) then
           j = 1151 ! Have to search by record number.
           call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
                  unpack, k, gfld, rc)
           err_msg="COULD NOT FIND MIN VEGETATION FRACTION IN FILE. SET MINMAX_VGFRC_FROM_CLIMO=.TRUE."
           if (rc/=0) call error_handler(err_msg, rc)
         endif
       endif
    
       if (maxval(gfld%fld) > 2.0) gfld%fld = gfld%fld / 100.0
       print*,'vfrac min ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

     endif ! localpet == 0

     print*,"- CALL FieldScatter FOR INPUT GRID MIN VEG GREENNESS."
     call ESMF_FieldScatter(min_veg_greenness_input_grid,dummy2d_8, rootpet=0, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)
   
     if (localpet == 0) then

       print*,"- READ MAX VEG FRACTION."

       jdisc   = 2  ! Search for discipline - land products
       j = 1106 ! Have to search by record number.
       jpdtn   = pdt_num  ! Search for the product definition template number.
       jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
       jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
       jpdt(2) = 4  ! Sec4/oct 11 - parameter number - vegetation
       jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
       unpack=.true.

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, rc)
       if (rc /= 0) then
         j = 1102 ! Have to search by record number.
         call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
              unpack, k, gfld, rc)
         if (rc /= 0) then
           j = 1152 ! Have to search by record number.
           call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
                unpack, k, gfld, rc)
           err_msg="COULD NOT FIND MAX VEGETATION FRACTION IN FILE. SET MINMAX_VGFRC_FROM_CLIMO=.TRUE."
           if (rc <= 0) call error_handler(err_msg, rc)
         endif
       endif
    
       if (maxval(gfld%fld) > 2.0) gfld%fld = gfld%fld / 100.0
!      print*,'vfrac max ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

     endif !localpet==0

     print*,"- CALL FieldScatter FOR INPUT GRID MAX VEG GREENNESS."
     call ESMF_FieldScatter(max_veg_greenness_input_grid,dummy2d_8,rootpet=0, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)

   endif !minmax_vgfrc_from_climo
 
   if (.not. lai_from_climo) then

     if (localpet == 0) then

       print*,"- READ LAI."

       jdisc   = 0  ! Search for discipline - meteorological products
       j = 0        ! Search at beginning of file.
       jpdtn   = pdt_num  ! Search for the product definition template number.
       jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
       jpdt(1) = 7   ! Sec4/oct 10 - parameter category - thermo stability indices
       jpdt(2) = 198 ! Sec4/oct 11 - parameter number - leaf area index
       jpdt(10) = 1  ! Sec4/oct 23 - type of level - ground surface
       unpack=.true.

       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

       err_msg="COULD NOT FIND LAI IN FILE. SET LAI_FROM_CLIMO=.TRUE."
       if (rc /= 0) call error_handler(err_msg, rc)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))

     endif !localpet==0

     print*,"- CALL FieldScatter FOR INPUT GRID LAI."
     call ESMF_FieldScatter(lai_input_grid,dummy2d_8,rootpet=0, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldScatter", rc)

   endif ! lai

   if (localpet == 0) then

     print*,"- READ SEAICE DEPTH."
     vname="hice"
     slev=":surface:" 
     call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)                 

     jdisc   = 10  ! Search for discipline - ocean products
     j = 0         ! Search at beginning of file.
     jpdtn   = pdt_num ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
     jpdt(1) = 2  ! Sec4/oct 10 - parameter category - ice
     jpdt(2) = 1  ! Sec4/oct 11 - parameter number - thickness
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0 ) then
       call handle_grib_error(vname, slev ,method,value,varnum,rc,var8=dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL BE"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8(:,:) = 0.0
       endif
     else
!      print*,'hice ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
     endif

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

! No test data contained this field. So could not test with g2 library.
     rc = 1
     if (rc /= 0) then
        call handle_grib_error(vname, slev ,method,value,varnum,rc, var8=dummy2d_8)
        if (rc==1) then ! missing_var_method == skip or no entry in varmap table
          print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                     " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                        "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
          dummy2d_8 = 0.0
        endif
     endif
     print*,'tprcp ',maxval(dummy2d_8),minval(dummy2d_8)

   endif ! tprcp

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

! No sample data contained this field, so could not test g2lib.
     rc = 1
     if (rc /= 0) then
       call handle_grib_error(vname, slev ,method,value,varnum,rc, var8=dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                   " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8(:,:) = 0.0
       endif
     endif
     print*,'ffmm ',maxval(dummy2d_8),minval(dummy2d_8)

   endif ! ffmm

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

     jdisc   = 0  ! Search for discipline - meteorological products
     j = 0        ! Search at beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
     jpdt(1) = 2  ! Sec4/oct 10 - parameter category - momentum
     jpdt(2) = 30 ! Sec4/oct 11 - parameter number - friction velocity
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     if (rc /= 0) then
       jpdt(2) = 197  ! oct 11 - param number - friction vel.
       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     endif

     if (rc == 0) then
!      print*,'fricv ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
     else
       call handle_grib_error(vname, slev ,method,value,varnum,rc, var8=dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL "//&
                   "REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8(:,:) = 0.0
       endif
     endif

   endif ! ustar

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

     rc = -1 ! None of the test cases have this record. Can't test with g2lib.
     if (rc /= 0) then
       call handle_grib_error(vname, slev ,method,value,varnum,rc, var8=dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL NOT"//&
                   " BE WRITTEN TO THE INPUT FILE. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8(:,:) = 0.0
       endif
     endif
     print*,'f10m ',maxval(dummy2d_8),minval(dummy2d_8)

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

     jdisc   = 2  ! Search for discipline - land products
     j = 0        ! Search from beginning of file
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
     jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
     jpdt(2) = 13 ! Sec4/oct 11 - parameter number - canopy water
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0 ) then
       jpdt(2) = 196 ! Sec4/oct 11 - param number - canopy water
       call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)
     endif

     if (rc == 0 ) then
       print*,'cnwat ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
       call check_cnwat(dummy2d_8)
     else
       call handle_grib_error(vname, slev ,method,value,varnum,rc, var8=dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8 = 0.0
       endif
     endif

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

     jdisc   = 2  ! Search for discipline - land products
     j = 0        ! Search from beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
     jpdt(1) = 0  ! Sec4/oct 10 - parameter category - veg/biomass
     jpdt(2) = 1  ! Sec4/oct 11 - parameter number - surface roughness
     jpdt(10) = 1 ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0 ) then
       call handle_grib_error(vname, slev ,method,value,varnum,rc, var8= dummy2d_8)
       if (rc==1) then ! missing_var_method == skip or no entry in varmap table
         print*, "WARNING: "//trim(vname)//" NOT AVAILABLE IN FILE. THIS FIELD WILL BE"//&
                   " REPLACED WITH CLIMO. SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS IS NOT DESIRABLE."
         dummy2d_8(:,:) = 0.0
       endif
     else
       gfld%fld = gfld%fld * 10.0 ! Grib files have z0 (m), but fv3 expects z0(cm)
!      print*,'sfcr ', maxval(gfld%fld),minval(gfld%fld)
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
     endif

   endif

   print*,"- CALL FieldScatter FOR INPUT GRID Z0."
   call ESMF_FieldScatter(z0_input_grid,dummy2d_8, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
 
   if (localpet == 0) then
     print*,"- READ LIQUID SOIL MOISTURE."
     vname = "soill"
     vname_file = ":SOILL:"
     call read_grib_soil(vname,vname_file,lugb, pdt_num,dummy3d) !!! NEED TO HANDLE 
                                                         !!! SOIL LEVELS
   endif

   print*,"- CALL FieldScatter FOR INPUT LIQUID SOIL MOISTURE."
   call ESMF_FieldScatter(soilm_liq_input_grid, dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
 
   if (localpet == 0) then
     print*,"- READ TOTAL SOIL MOISTURE."
     vname = "soilw"
     vname_file = "var2_2_1_"         ! the var number instead
     call read_grib_soil(vname,vname_file,lugb, pdt_num,dummy3d)
   endif
 
   print*,"- CALL FieldScatter FOR INPUT TOTAL SOIL MOISTURE."
   call ESMF_FieldScatter(soilm_tot_input_grid, dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)
    
!----------------------------------------------------------------------------------------
! Vegetation type is not available in some files.  However, it is needed to identify
! permanent land ice points.  At land ice, the total soil moisture is a flag value of
! '1'. Use this flag as a temporary solution.
!----------------------------------------------------------------------------------------

   print*, "- CALL FieldGather for INPUT SOIL TYPE."
   call ESMF_FieldGather(soil_type_input_grid, dummy2d_82, rootPet=0, tile=1, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGather", rc)

   if (localpet == 0) then

     print*,"- READ VEG TYPE."
   
     jdisc   = 2  ! Search for discipline - land products
     j = 0        ! Search from beginning of file.
     jpdtn   = pdt_num  ! Search for the product definition template number.
     jpdt    = -9999  ! Initialize array of values in product definition template Sec4.
     jpdt(1) = 0   ! Sec4/oct 10 - parameter category - veg/biomass
     jpdt(2) = 198 ! Sec4/oct 11 - parameter number - vegetation type
     jpdt(10) = 1  ! Sec4/oct 23 - type of level - ground surface
     unpack=.true.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc)

     if (rc /= 0 ) then
       if (.not. vgtyp_from_climo) then
         call error_handler("COULD NOT FIND VEGETATION TYPE IN FILE. PLEASE SET VGTYP_FROM_CLIMO=.TRUE. . EXITING", rc)
       else ! Set input veg type at land ice from soil moisture flag (1.0).
         do j = 1, j_input
          do i = 1, i_input
            dummy2d_8(i,j) = 0.0
            if(slmsk_save(i,j) == 1 .and. dummy3d(i,j,1) > 0.99) &  ! land ice indicated by
                                                                    ! soil moisture flag of '1'.
            dummy2d_8(i,j) = real(veg_type_landice_input,esmf_kind_r8)
          enddo
         enddo    
       endif
     else  ! found vtype in file.
       dummy2d_8 = reshape(gfld%fld , (/i_input,j_input/))
     endif

     if (trim(external_model) .ne. "GFS") then
       do j = 1, j_input
       do i = 1,i_input
         if (dummy2d_8(i,j) == 15.0_esmf_kind_r8 .and. slmsk_save(i,j) == 1) then
           if (dummy3d(i,j,1) < 0.6) then 
             dummy2d_8(i,j) = real(veg_type_landice_input,esmf_kind_r8)
           elseif (dummy3d(i,j,1) > 0.99) then
             slmsk_save(i,j) = 0
             dummy2d_8(i,j) = 0.0_esmf_kind_r8
             dummy2d_82(i,j) = 0.0_esmf_kind_r8
           endif
         elseif (dummy2d_8(i,j) == 17.0_esmf_kind_r8 .and. slmsk_save(i,j)==0) then
           dummy2d_8(i,j) = 0.0_esmf_kind_r8
         endif
       enddo
       enddo
     endif     

!    print*,'vgtyp ',maxval(dummy2d_8),minval(dummy2d_8)

   endif ! read veg type

   print*,"- CALL FieldScatter FOR INPUT VEG TYPE."
   call ESMF_FieldScatter(veg_type_input_grid, dummy2d_8, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   print*,"- CALL FieldScatter FOR INPUT SOIL TYPE."
   call ESMF_FieldScatter(soil_type_input_grid, dummy2d_82, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   deallocate(dummy2d_82)

   print*,"- CALL FieldScatter FOR INPUT LANDSEA MASK."
   call ESMF_FieldScatter(landsea_mask_input_grid,real(slmsk_save,esmf_kind_r8),rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

!---------------------------------------------------------------------------------
! At open water (slmsk==0), the soil temperature array is not used and set
! to the filler value of SST.  At lake/sea ice points (slmsk=2), the soil 
! temperature array holds ice column temperature.  This field is not available
! in the grib data, so set to a default value.
!---------------------------------------------------------------------------------

   if (localpet == 0) then
     print*,"- READ SOIL TEMPERATURE."
     vname = "soilt"
     vname_file = ":TSOIL:"
     call read_grib_soil(vname,vname_file,lugb,pdt_num,dummy3d)
     call check_soilt(dummy3d,slmsk_save,tsk_save)
     deallocate(tsk_save)
   endif

   deallocate(slmsk_save)

   print*,"- CALL FieldScatter FOR INPUT SOIL TEMPERATURE."
   call ESMF_FieldScatter(soil_temp_input_grid, dummy3d, rootpet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldScatter", rc)

   deallocate(dummy3d)
   deallocate(dummy2d_8)
 
   if (localpet == 0) call baclose(lugb, rc)

 end subroutine read_input_sfc_grib2_file
   
!> Read nst data from these netcdf formatted fv3 files: tiled history,
!! tiled warm restart, and gaussian history.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read input grid nst data from fv3 gaussian nemsio history file or
!! spectral GFS nemsio file.
!!
!! @note The spectral GFS nst data is in a separate file from
!! the surface data.  The fv3 surface and nst data are in a
!! single file.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author George Gayno NCEP/EMC   
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

!> Read a record from a netcdf file
!!
!! @param [in] field  name of field to be read 
!! @param [in] tile_num  grid tile number
!! @param [in] imo i-dimension of field
!! @param [in] jmo j-dimension of field
!! @param [in] lmo number of vertical levels of field
!! @param [out] sfcdata 1-d array containing field data
!! @param [out] sfcdata_3d  3-d array containing field data
!! @author George Gayno NCEP/EMC   
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
 
!> Read winds from a grib2 file.  Rotate winds
!! to be earth relative if necessary.
!!
!! @param [inout] u  u-component wind
!! @param [inout] v  v-component wind
!! @param[in] localpet  ESMF local persistent execution thread
!! @param[in] octet_23 Section 4/Octet 23 - Type of first fixed surface.
!! @param[in] rlevs Array of atmospheric level values
!! @param[in] lugb Logical unit number of GRIB2 file.
!! @param[in] pdt_num Product definition template number.
!! @author Larissa Reames
 subroutine read_winds(u,v,localpet,octet_23,rlevs,lugb,pdt_num)

 use grib_mod
 use program_setup, only      : get_var_cond

 implicit none

 integer, intent(in)                                  :: localpet, lugb
 integer, intent(in)                                  :: pdt_num, octet_23

 real(esmf_kind_r8), intent(inout), allocatable       :: u(:,:,:),v(:,:,:)
 real(esmf_kind_r8), intent(in), dimension(lev_input) :: rlevs

 real(esmf_kind_r4), dimension(i_input,j_input)  :: alpha
 real(esmf_kind_r8), dimension(i_input,j_input)  :: lon, lat
 real(esmf_kind_r4), allocatable                 :: u_tmp(:,:),v_tmp(:,:)
 real(esmf_kind_r8), allocatable                 :: dum2d(:,:)
 real(esmf_kind_r4), dimension(i_input,j_input)  :: ws,wd
 real(esmf_kind_r4)                      :: value_u, value_v,lov,latin1,latin2
 real(esmf_kind_r8)                      :: d2r

 integer                                 :: varnum_u, varnum_v, vlev, &
                                            error, iret
 integer                                 :: j, k, lugi, jgdtn, jpdtn
 integer                                 :: jdisc, jids(200), jgdt(200), jpdt(200)

 character(len=20)                       :: vname
 character(len=50)                       :: method_u, method_v

 logical                                 :: unpack

 type(gribfield)                         :: gfld

 d2r=acos(-1.0_esmf_kind_r8) / 180.0_esmf_kind_r8
 if (localpet==0) then
   allocate(u(i_input,j_input,lev_input))
   allocate(v(i_input,j_input,lev_input))
 else
   allocate(u(0,0,0))
   allocate(v(0,0,0))
 endif

 vname = "u"
 call get_var_cond(vname,this_miss_var_method=method_u, this_miss_var_value=value_u, &
                       loc=varnum_u)
 vname = "v"
 call get_var_cond(vname,this_miss_var_method=method_v, this_miss_var_value=value_v, &
                       loc=varnum_v)

 print*,"- CALL FieldGather FOR INPUT GRID LONGITUDE"
 call ESMF_FieldGather(longitude_input_grid, lon, rootPet=0, tile=1, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)

 print*,"- CALL FieldGather FOR INPUT GRID LATITUDE"
 call ESMF_FieldGather(latitude_input_grid, lat, rootPet=0, tile=1, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)

 if (localpet==0) then

   lugi    = 0     ! index file unit number
   jdisc   = 0     ! search for discipline - meteorological products
   j = 0           ! search at beginning of file.
   jpdt    = -9999  ! array of values in product definition template, set to wildcard
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template, set to wildcard
   jgdtn   = -1     ! search for any grid definition number.
   jpdtn   =  pdt_num ! Search for the product definition template number.
   unpack=.false.

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)
    
   if (iret /= 0) call error_handler("ERROR READING GRIB2 FILE.", iret)

   if (gfld%igdtnum == 32769) then ! grid definition template number - rotated lat/lon grid

     latin1 = real(float(gfld%igdtmpl(15))/1.0E6, kind=esmf_kind_r4)
     lov = real(float(gfld%igdtmpl(16))/1.0E6, kind=esmf_kind_r4)

     print*, "- CALL CALCALPHA_ROTLATLON with center lat,lon = ",latin1,lov
     call calcalpha_rotlatlon(lat,lon,latin1,lov,alpha)

   elseif (gfld%igdtnum == 30) then ! grid definition template number - lambert conformal grid.

     lov = real(float(gfld%igdtmpl(14))/1.0E6, kind=esmf_kind_r4)
     latin1 = real(float(gfld%igdtmpl(19))/1.0E6, kind=esmf_kind_r4)
     latin2 = real(float(gfld%igdtmpl(20))/1.0E6, kind=esmf_kind_r4)

     print*, "- CALL GRIDROT for LC grid with lov,latin1/2 = ",lov,latin1,latin2
     call gridrot(lov,latin1,latin2,lon,alpha)

   endif

   jpdt(10) = octet_23 ! Sec4/oct 23 - type of level.

   unpack=.true.

   allocate(dum2d(i_input,j_input))
   allocate(u_tmp(i_input,j_input))
   allocate(v_tmp(i_input,j_input))

   do vlev = 1, lev_input

     vname = ":UGRD:"

     jpdt(1) = 2  ! Sec4/oct 10 - parameter category - momentum
     jpdt(2) = 2  ! Sec4/oct 11 - parameter number - u-wind
     jpdt(12) = nint(rlevs(vlev)) ! Sect4/octs 25-28 - scaled value of fixed surface.

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret /= 0) then
        call handle_grib_error(vname, slevs(vlev),method_u,value_u,varnum_u,iret,var=u_tmp)
        if (iret==1) then ! missing_var_method == skip
          call error_handler("READING IN U AT LEVEL "//trim(slevs(vlev))//". SET A FILL "// &
                        "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
        endif
     else
       dum2d = reshape(gfld%fld, (/i_input,j_input/) )
       u_tmp(:,:) = real(dum2d, kind=esmf_kind_r4)
     endif

     vname = ":VGRD:"

     jpdt(2) = 3  ! Sec4/oct 11 - parameter number - v-wind

     call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

     if (iret /= 0) then
        call handle_grib_error(vname, slevs(vlev),method_v,value_v,varnum_v,iret,var=v_tmp)
        if (iret==1) then ! missing_var_method == skip
          call error_handler("READING IN V AT LEVEL "//trim(slevs(vlev))//". SET A FILL "// &
                        "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
        endif
     else
       dum2d = reshape(gfld%fld, (/i_input,j_input/) )
       v_tmp(:,:) = real(dum2d, kind=esmf_kind_r4)
     endif

     deallocate(dum2d)

      if (gfld%igdtnum == 0) then ! grid definition template number - lat/lon grid
        if (external_model == 'UKMET') then
          u(:,:,vlev) = u_tmp
          v(:,:,vlev) = (v_tmp(:,2:jp1_input) + v_tmp(:,1:j_input))/2
        else
          u(:,:,vlev) = u_tmp
          v(:,:,vlev) = v_tmp
        endif
      else if (gfld%igdtnum == 32769) then ! grid definition template number - rotated lat/lon grid
        ws = sqrt(u_tmp**2 + v_tmp**2)
        wd = real((atan2(-u_tmp,-v_tmp) / d2r), kind=esmf_kind_r4) ! calculate grid-relative wind direction
        wd = real((wd + alpha + 180.0), kind=esmf_kind_r4) ! Rotate from grid- to earth-relative direction
        wd = real((270.0 - wd), kind=esmf_kind_r4) ! Convert from meteorological (true N) to mathematical direction
        u(:,:,vlev) = -ws*cos(wd*d2r)
        v(:,:,vlev) = -ws*sin(wd*d2r)
      else
        u(:,:,vlev) = real(u_tmp * cos(alpha) + v_tmp * sin(alpha),esmf_kind_r8)
        v(:,:,vlev) = real(v_tmp * cos(alpha) - u_tmp * sin(alpha),esmf_kind_r8)
      endif

      print*, 'max, min U ', minval(u(:,:,vlev)), maxval(u(:,:,vlev))
      print*, 'max, min V ', minval(v(:,:,vlev)), maxval(v(:,:,vlev))
    enddo
 endif

end subroutine read_winds

!> Convert winds from 2-d to 3-d components.
!!
!! @author George Gayno NCEP/EMC   
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
 
!> Compute grid rotation angle for non-latlon grids.
!!
!! @note The original gridrot subroutine was specific to polar
!! stereographic grids.  We need to compute it for Lambert Conformal
!! grids. So we need lat1,lat2.  This follows the ncl_ncarg source
!! code: ncl_ncarg-6.6.2/ni/src/ncl/GetGrids.c
!!
!! @param [in] lov  orientation angle
!! @param [in] latin1  first tangent latitude
!! @param [in] latin2  second tangent latitude
!! @param [in] lon     longitude
!! @param [inout] rot  rotation angle
!! @author Larissa Reames
subroutine gridrot(lov,latin1,latin2,lon,rot)

  use model_grid, only                : i_input,j_input
  implicit none


  real(esmf_kind_r4), intent(in)      :: lov,latin1,latin2
  real(esmf_kind_r4), intent(inout)   :: rot(i_input,j_input)
  real(esmf_kind_r8), intent(in)      :: lon(i_input,j_input)

  real(esmf_kind_r4)                  :: trot(i_input,j_input), tlon(i_input,j_input)
  real(esmf_kind_r4)                  :: dtor = 3.14159265359/180.0_esmf_kind_r4
  real(esmf_kind_r4)                  :: an
  !trot_tmp = real(lon,esmf_kind_r4)-lov
  !trot = trot_tmp
  !where(trot_tmp > 180.0) trot = trot-360.0_esmf_kind_r4
  !where(trot_tmp < -180.0) trot = trot-360.0_esmf_kind_r4

  if ( (latin1 - latin2) .lt. 0.000001 ) then
        an = sin(latin1*dtor)
  else
        an = real(log( cos(latin1*dtor) / cos(latin2*dtor) ) / &
             log( tan(dtor*(90.0-latin1)/2.) / tan(dtor*(90.0-latin2)/2.)), kind=esmf_kind_r4)
  end if

  tlon = real((mod(lon - lov + 180. + 3600., 360.) - 180.), kind=esmf_kind_r4)
  trot = an * tlon

  rot = trot * dtor

end subroutine gridrot

!> Calculate rotation angle for rotated latlon grids.
!! Needed to convert to earth-relative winds.
!!
!! @param [in] latgrid  grid latitudes
!! @param [in] longrid  grid longitudes
!! @param [in] cenlat   center latitude
!! @param [in] cenlon   center longitude
!! @param [out] alpha   grid rotation angle
!! @author Larissa Reames
subroutine calcalpha_rotlatlon(latgrid,longrid,cenlat,cenlon,alpha)

  use model_grid, only                : i_input,j_input
  implicit none

  real(esmf_kind_r8), intent(in)      :: latgrid(i_input,j_input), &
                                         longrid(i_input,j_input)
  real(esmf_kind_r4), intent(in)      :: cenlat, cenlon
  real(esmf_kind_r4), intent(out)     :: alpha(i_input,j_input)

  ! Variables local to subroutine
  real(esmf_kind_r8)             :: D2R,lon0_r,lat0_r,sphi0,cphi0
  real(esmf_kind_r8), DIMENSION(i_input,j_input) :: tlat,tlon,tph,sinalpha

  D2R = acos(-1.0_esmf_kind_r8) /  180.0_esmf_kind_r8
  if (cenlon .lt. 0) then
      lon0_r = (cenlon + 360.0)*D2R
  else
      lon0_r = cenlon*D2R
  end if
  lat0_r=cenlat*D2R
  sphi0=sin(lat0_r)
  cphi0=cos(lat0_r)

  ! deal with input lat/lon
  tlat = latgrid * D2R
  tlon = longrid * D2R

  ! Calculate alpha (rotation angle)
  tlon = -tlon + lon0_r
  tph  = asin(cphi0*sin(tlat) - sphi0*cos(tlat)*cos(tlon))
  sinalpha = sphi0 * sin(tlon) / cos(tph)
  alpha = real((-asin(sinalpha)/D2R), kind=esmf_kind_r4)
  ! returns alpha in degrees
end subroutine calcalpha_rotlatlon

!> Handle GRIB2 read error based on the user selected
!! method in the varmap file.
!!
!! @param [in] vname  grib2 variable name
!! @param [in] lev    grib2 variable level
!! @param [in] method  how missing data is handled
!! @param [in] value   fill value for missing data
!! @param [in] varnum  grib2 variable number
!! @param [inout] iret  return status code
!! @param [inout] var   4-byte array of corrected data
!! @param [inout] var8  8-byte array of corrected data
!! @param [inout] var3d 3-d array of corrected data
!! @author Larissa Reames
subroutine handle_grib_error(vname,lev,method,value,varnum, iret,var,var8,var3d)

  use, intrinsic :: ieee_arithmetic

  implicit none
  
  real(esmf_kind_r4), intent(in)    :: value
  real(esmf_kind_r4), intent(inout), optional :: var(:,:)
  real(esmf_kind_r8), intent(inout), optional :: var8(:,:)
  real(esmf_kind_r8), intent(inout), optional  :: var3d(:,:,:)
  
  character(len=20), intent(in)     :: vname, lev, method
  character(len=200)                :: err_msg
  
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
    print*, "WARNING: ,", trim(vname), " NOT AVAILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO FILL VALUE OF ", value
    if(present(var)) var(:,:) = value
    if(present(var8)) var8(:,:) = value
    if(present(var3d)) var3d(:,:,:) = value
  elseif (trim(method) == "set_to_NaN") then
    print*, "WARNING: ,", trim(vname), " NOT AVAILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO NaNs"
    if(present(var)) var(:,:) = ieee_value(var,IEEE_QUIET_NAN)
    if(present(var8)) var8(:,:) = ieee_value(var8,IEEE_QUIET_NAN)
    if(present(var3d)) var3d(:,:,:) = ieee_value(var3d,IEEE_QUIET_NAN)
  elseif (trim(method) == "stop") then
    err_msg="READING " // trim(vname) // " at level " //lev// ". TO MAKE THIS NON-" // &
            "FATAL, CHANGE STOP TO SKIP FOR THIS VARIABLE IN YOUR VARMAP FILE."
    call error_handler(err_msg, iret)
  elseif (trim(method) == "intrp") then
    print*, "WARNING: ,"//trim(vname)//" NOT AVAILABLE AT LEVEL "//trim(lev)// &
          ". WILL INTERPOLATE INTERSPERSED MISSING LEVELS AND/OR FILL MISSING"//&
          " LEVELS AT EDGES."
  else
    err_msg="ERROR USING MISSING_VAR_METHOD. PLEASE SET VALUES IN" // &
            " VARMAP TABLE TO ONE OF: set_to_fill, set_to_NaN, intrp, skip, or stop."
    call error_handler(err_msg, 1)
  endif

end subroutine handle_grib_error

!> Read soil temperature and soil moisture fields from a GRIB2 file.
!!
!! @param [in] vname         variable name in varmap table
!! @param [in] vname_file    variable name in grib2 file
!! @param [in] lugb          logical unit number for surface grib2 file
!! @param [in] pdt_num       product definition template number.
!! @param [inout] dummy3d    array of soil data
!! @author George Gayno NCEP/EMC   
 subroutine read_grib_soil(vname, vname_file, lugb, pdt_num, dummy3d)
  
 use grib_mod

 implicit none
  
 character(len=20), intent(in)           :: vname,vname_file
  
 integer, intent(in)                     :: lugb, pdt_num

 real(esmf_kind_r8), intent(inout)       :: dummy3d(:,:,:)
  
 character(len=50)                       :: slevs(lsoil_input)
 character(len=50)                       :: method

 integer                                 :: varnum, i, j, k, rc, rc2
 integer                                 :: jdisc, jgdtn, jpdtn, lugi
 integer                                 :: jids(200), jgdt(200), jpdt(200)
 integer                                 :: iscale1, iscale2

 logical                                 :: unpack

 real(esmf_kind_r4), allocatable         :: dummy2d(:,:)
 real(esmf_kind_r4)                      :: value

 type(gribfield)                         :: gfld

 allocate(dummy2d(i_input,j_input))

 if(lsoil_input == 4) then
    slevs = (/character(24)::':0-0.1 m below ground:', ':0.1-0.4 m below ground:', &
                             ':0.4-1 m below ground:', ':1-2 m below ground:'/)
 elseif(lsoil_input == 9) then
    slevs = (/character(26)::':0-0 m below ground',':0.01-0.01 m below ground:',':0.04-0.04 m below ground:', &
        ':0.1-0.1 m below ground:',':0.3-0.3 m below ground:',':0.6-0.6 m below ground:', &
        ':1-1 m below ground:',':1.6-1.6 m below ground:',':3-3 m below ground:'/)
 else
    rc = -1
    call error_handler("reading soil levels. File must have 4 or 9 soil levels.", rc)
 endif
 
 call get_var_cond(vname,this_miss_var_method=method,this_miss_var_value=value, &
                         loc=varnum)

 lugi = 0        ! unit number for index file
 jdisc   = 2     ! search for discipline - land products
 j = 0           ! search at beginning of file.
 jpdt    = -9999  ! array of values in product definition template 4.n
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jgdtn   = -1     ! search for any grid definition number.
 jpdtn   = pdt_num  ! Search for the product definition template number.
 jpdt(1) = 0        ! Section 4/Octet 10 - parameter category - veg/biomass
 if (trim(vname) == 'soilt') jpdt(2) = 2    ! Section 4/Octet 11 - parameter number - soil temp
 if (trim(vname) == 'soilw') jpdt(2) = 192  ! Section 4/Octet 11 - parameter number - total soilm
 if (trim(vname) == 'soill') then
   jpdt(1) = 3    ! Section 4/Octet 10 - soil products
   jpdt(2) = 192  ! Section 4/Octet 11 - parameter number - liquid soilm
 endif
 jpdt(10) = 106 ! Section 4/Octet 23 - depth below ground 
 jpdt(13) = 106 ! Section 4/Octet 29 - depth below ground
 unpack=.true.

 do i = 1,lsoil_input

   call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, rc2)

   if (rc2 /= 0) then  ! record not found.
     call handle_grib_error(vname_file, slevs(i),method,value,varnum,rc,var=dummy2d)
     if (rc==1 .and. trim(vname) /= "soill") then 
       ! missing_var_method == skip or no entry in varmap table
       call error_handler("READING IN "//trim(vname)//". SET A FILL "// &
                      "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",rc)
     elseif (rc==1) then
       dummy3d(:,:,:) = 0.0_esmf_kind_r8
       return
     endif
   endif

   if (rc2 == 0) then ! record found. 
     iscale1 = 10 ** gfld%ipdtmpl(11)
     iscale2 = 10 ** gfld%ipdtmpl(14)
!    print*,'getgb2 top of soil layer in m ', float(gfld%ipdtmpl(12))/float(iscale1)
!    print*,'getgb2 bot of soil layer in m ', float(gfld%ipdtmpl(15))/float(iscale2)
     dummy2d = reshape(real(gfld%fld,kind=esmf_kind_r4), (/i_input,j_input/) )
   endif 

   j = k

   dummy3d(:,:,i) = real(dummy2d,esmf_kind_r8)

 enddo

 deallocate(dummy2d)

 end subroutine read_grib_soil

!> Free up memory associated with atm data.
!!
!! @author George Gayno NCEP/EMC   
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

 do n = 1, num_tracers_input
   call ESMF_FieldDestroy(tracers_input_grid(n), rc=rc)
 enddo
 deallocate(tracers_input_grid)

 end subroutine cleanup_input_atm_data

!> Free up memory associated with nst data.
!!
!! @author George Gayno NCEP/EMC   
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

!> Free up memory associated with sfc data.
!!
!! @author George Gayno NCEP/EMC   
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
 if (.not. vgfrc_from_climo) then
   call ESMF_FieldDestroy(veg_greenness_input_grid, rc=rc)
 endif
 if (.not. minmax_vgfrc_from_climo) then
   call ESMF_FieldDestroy(min_veg_greenness_input_grid, rc=rc)
   call ESMF_FieldDestroy(max_veg_greenness_input_grid, rc=rc)
 endif
 if (.not. lai_from_climo) then
   call ESMF_FieldDestroy(lai_input_grid, rc=rc)
 endif

 end subroutine cleanup_input_sfc_data

!> Sort an array of values.
!!
!! @param a      the sorted array
!! @param first  the first value of sorted array
!! @param last   the last value of sorted array
!! @author Jili Dong NOAA/EMC
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

!> Check for and replace certain values in soil temperature.
!> At open water points (landmask=0) use skin temperature as
!> a filler value. At land points (landmask=1) with excessive
!> soil temperature, replace soil temperature with skin temperature. 
!> In GEFSv12.0 data there are some erroneous missing values at
!> land points which this corrects. At sea ice points (landmask=2),
!> store a default ice column temperature because grib2 files do not 
!> have ice column temperature which FV3 expects at these points.
!!
!! @param soilt    [inout] 3-dimensional soil temperature arrray
!! @param landmask [in]    landmask of the input grid
!! @param skint    [in]    2-dimensional skin temperature array
!! @author Larissa Reames CIMMS/NSSL

subroutine check_soilt(soilt, landmask, skint)
  implicit none
  real(esmf_kind_r8), intent(inout) ::  soilt(i_input,j_input,lsoil_input)
  real(esmf_kind_r8), intent(in)    ::  skint(i_input,j_input)
  integer(esmf_kind_i4), intent(in)    ::  landmask(i_input,j_input)
  
  integer                           :: i, j, k

  do k=1,lsoil_input
    do j = 1, j_input
      do i = 1, i_input
        if (landmask(i,j) == 0_esmf_kind_i4 ) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 1_esmf_kind_i4 .and. soilt(i,j,k) > 350.0_esmf_kind_r8) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 2_esmf_kind_i4 ) then 
          soilt(i,j,k) = ICET_DEFAULT
        endif
      enddo
    enddo
  enddo
end subroutine check_soilt

!> When using GEFS data, some points on the target grid have 
!> unreasonable canpy moisture content, so zero out any 
!> locations with unrealistic canopy moisture values (>0.5).
!!
!! @param cnwat [input] 2-dimensional canopy moisture content
!! @author Larissa Reames CIMMS/NSSL

subroutine check_cnwat(cnwat)
  implicit none 
  real(esmf_kind_r8), intent(inout) :: cnwat(i_input,j_input)
  
  real(esmf_kind_r8)                :: max_cnwat = 0.5
  
  integer :: i, j

  do i = 1,i_input
    do j = 1,j_input
      if (cnwat(i,j) .gt. max_cnwat) cnwat(i,j) = 0.0_esmf_kind_r8
    enddo
  enddo
end subroutine check_cnwat




!> Pressure to presure vertical interpolation for tracers with linear or lnP
!> interpolation. Input tracers on pres levels are interpolated 
!> to the target output pressure levels. The matching levels of input and 
!> output will keep the same. Extrapolation is also allowed but needs
!> caution. The routine is mostly for GFSV16 combined grib2 input when spfh has
!> missing levels in low and mid troposphere from U/T/HGT/DZDT. 
!!
!! @param [in] ppin  1d input pres levs
!! @param [in] xxin  1d input tracer
!! @param [in] npin  number of input levs 
!! @param [in] ppout 1d target pres levs 
!! @param [out] xxout 1d interpolated tracer
!! @param [in] npout number of target levs 
!! @param [in] linlog interp method.1:linear;not 1:log;neg:extrp allowed
!! @param [in] xmsg  fill values of missing levels (-999.0)
!! @param [out] ier  error status. non 0: failed interpolation
!! @author Jili Dong NCEP/EMC  
!! @date 2021/07/30

SUBROUTINE DINT2P(PPIN,XXIN,NPIN,PPOUT,XXOUT,NPOUT   &
                      ,LINLOG,XMSG,IER)
      IMPLICIT NONE

! NCL code for pressure level interpolation
!
! This code was designed for one simple task. It has since
! been mangled and abused for assorted reasons. For example,
! early gfortran compilers had some issues with automatic arrays.
! Hence, the C-Wrapper was used to create 'work' arrays which
! were then passed to this code.  The original focused (non-NCL) 
! task was to handle PPIN & PPOUT that had the same 'monotonicity.' 
! Extra code was added to handle the more general case. 
! Blah-Blah:  Punch line: it is embarrassingly convoluted!!!
!
!                                                ! input types
      INTEGER NPIN,NPOUT,LINLOG,IER
      real*8 PPIN(NPIN),XXIN(NPIN),PPOUT(NPOUT),XMSG
                                                ! output
      real*8 XXOUT(NPOUT)
                                                ! work
      real*8 PIN(NPIN),XIN(NPIN),P(NPIN),X(NPIN)
      real*8 POUT(NPOUT),XOUT(NPOUT)

! local
      INTEGER NP,NL,NLMAX,NLSAVE,NP1,NO1,N1,N2,LOGLIN,   &
             NLSTRT
      real*8 SLOPE,PA,PB,PC

      LOGLIN = ABS(LINLOG)

! error check: enough points: pressures consistency?

      IER = 0
      IF (NPOUT.GT.0) THEN
          DO NP = 1,NPOUT
              XXOUT(NP) = XMSG
          END DO
      END IF
! Jili Dong input levels have to be the same as output levels:
! we only interpolate for levels with missing variables
!      IF (.not. all(PPIN .eq. PPOUT)) IER = IER+1

      IF (NPIN.LT.2 .OR. NPOUT.LT.1) IER = IER + 1

      IF (IER.NE.0) THEN
!          PRINT *,'INT2P: error exit: ier=',IER
          RETURN
      END IF

! should *input arrays* be reordered? want p(1) > p(2) > p(3) etc
! so that it will match order for which code was originally designed
! copy to 'work'  arrays

      NP1 = 0
      NO1 = 0
      IF (PPIN(1).LT.PPIN(2)) THEN
          NP1 = NPIN + 1
      END IF
      IF (PPOUT(1).LT.PPOUT(2)) THEN
          NO1 = NPOUT + 1
      END IF

      DO NP = 1,NPIN
          PIN(NP) = PPIN(ABS(NP1-NP))
          XIN(NP) = XXIN(ABS(NP1-NP))
      END DO

      DO NP = 1,NPOUT
          POUT(NP) = PPOUT(ABS(NO1-NP))
      END DO

! eliminate XIN levels with missing data. 
! .   This can happen with observational data.

      NL = 0
      DO NP = 1,NPIN
          IF (XIN(NP).NE.XMSG .AND. PIN(NP).NE.XMSG) THEN
              NL = NL + 1
              P(NL) = PIN(NP)
              X(NL) = XIN(NP)
          END IF
      END DO
      NLMAX = NL

                                                ! all missing data
      IF (NLMAX.LT.2) THEN
          IER = IER + 1000
          PRINT *,'INT2P: ier=',IER
          RETURN
      END IF

! ===============> pressure in decreasing order <================
! perform the interpolation  [pin(1)>pin(2)>...>pin(npin)]
!                                                      ( p ,x)
! ------------------------- p(nl+1), x(nl+1)   example (200,5)
! .
! ------------------------- pout(np), xout(np)         (250,?)
! .
! ------------------------- p(nl)  , x(nl)             (300,10)


! exact p-level matches
      NLSTRT = 1
      NLSAVE = 1
      DO NP = 1,NPOUT
          XOUT(NP) = XMSG
          DO NL = NLSTRT,NLMAX
              IF (POUT(NP).EQ.P(NL)) THEN
                  XOUT(NP) = X(NL)
                  NLSAVE = NL + 1
                  GO TO 10
              END IF
          END DO
   10     NLSTRT = NLSAVE
      END DO

      IF (LOGLIN.EQ.1) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      SLOPE = (X(NL)-X(NL+1))/ (P(NL)-P(NL+1))
                      XOUT(NP) = X(NL+1) + SLOPE* (POUT(NP)-P(NL+1))
                  END IF
              END DO
          END DO
      ELSE
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      PA = LOG(P(NL))
                      PB = LOG(POUT(NP))
! special case: In case someome inadvertently enter p=0.
                      if (p(nl+1).gt.0.d0) then
                          PC = LOG(P(NL+1))
                      else
                          PC = LOG(1.E-4)
                      end if

                      SLOPE = (X(NL)-X(NL+1))/ (PA-PC)
                      XOUT(NP) = X(NL+1) + SLOPE* (PB-PC)
                  END IF
              END DO
          END DO
      END IF

! extrapolate?
! . use the 'last' valid slope for extrapolating

      IF (LINLOG.LT.0) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX
                  IF (POUT(NP).GT.P(1)) THEN
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(2)-X(1))/ (P(2)-P(1))
                          XOUT(NP) = X(1) + SLOPE* (POUT(NP)-P(1))
                      ELSE
                          PA = LOG(P(2))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(1))
                          SLOPE = (X(2)-X(1))/ (PA-PC)
                          XOUT(NP) = X(1) + SLOPE* (PB-PC)
                      END IF
                  ELSE IF (POUT(NP).LT.P(NLMAX)) THEN
                      N1 = NLMAX
                      N2 = NLMAX - 1
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(N1)-X(N2))/ (P(N1)-P(N2))
                          XOUT(NP) = X(N1) + SLOPE* (POUT(NP)-P(N1))
                      ELSE
                          PA = LOG(P(N1))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(N2))
                          SLOPE = (X(N1)-X(N2))/ (PA-PC)
                          !XOUT(NP) = X(N1) + SLOPE* (PB-PC) !bug fixed below
                          XOUT(NP) = X(N1) + SLOPE* (PB-PA)
                      END IF
                  END IF
              END DO
          END DO
      END IF

! place results in the return array;
! .   possibly .... reverse to original order

      if (NO1.GT.0) THEN
          DO NP = 1,NPOUT
             n1 = ABS(NO1-NP)
             PPOUT(NP) = POUT(n1)
             XXOUT(NP) = XOUT(n1)
          END DO
      ELSE
          DO NP = 1,NPOUT
             PPOUT(NP) = POUT(NP)
             XXOUT(NP) = XOUT(NP)
          END DO
      END IF


      RETURN
      END SUBROUTINE DINT2P


 end module input_data
