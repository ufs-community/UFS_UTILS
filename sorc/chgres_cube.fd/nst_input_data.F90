module nst_input_data
!> @file
!! @brief Read  NST surface data from NEMSIO and NetCDF files.
!! @author George Gayno NCEP/EMC

!> Read nst data on the input grid.
!! Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!! 'history' files, fv3 gaussian history files, and spectral gfs
!! gaussian nemsio files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author George Gayno NCEP/EMC
 use esmf
 use netcdf
#ifdef CHGRES_ALL
 use nemsio_module
#endif

 use program_setup, only          : data_dir_input_grid, &
                                    sfc_files_input_grid, &
                                    nst_files_input_grid, &
                                    input_type
 
 use model_grid, only             : input_grid,        &
                                    i_input, j_input,  &
                                    ip1_input, jp1_input,  &
                                    num_tiles_input_grid
 
 use sfc_input_data, only         : lsoil_input, &
                                    read_fv3_grid_data_netcdf, &
                                    landsea_mask_input_grid

 use utilities, only              : error_handler
 implicit none 

! Fields associated with the nst model.

 type(esmf_field), public :: c_d_input_grid   !< Coefficient 2 to calculate d(tz)/d(ts)
 type(esmf_field), public :: c_0_input_grid   !< Coefficient 1 to calculate d(tz)/d(ts)
 type(esmf_field), public :: d_conv_input_grid   !< Thickness of free convectionlayer
 type(esmf_field), public :: dt_cool_input_grid   !< Sub-layer cooling amount
 type(esmf_field), public :: ifd_input_grid   !< Model mode index. 0-diurnalmodel not
                                                     !< started; 1-diurnal model
                                                     !started.
 type(esmf_field), public :: qrain_input_grid   !< Sensible heat flux due torainfall
 type(esmf_field), public :: tref_input_grid  !< Reference temperature
 type(esmf_field), public :: w_d_input_grid   !< Coefficient 4 to calculated(tz)/d(ts)
 type(esmf_field), public :: w_0_input_grid   !< Coefficient 3 to calculated(tz)/d(ts)
 type(esmf_field), public :: xs_input_grid   !< Salinity content in diurnalthermocline layer
 type(esmf_field), public :: xt_input_grid   !< Heat content in diurnalthermocline layer
 type(esmf_field), public :: xu_input_grid   !< u-current content in diurnalthermocline layer
 type(esmf_field), public :: xv_input_grid   !< v-current content in diurnalthermocline layer
 type(esmf_field), public :: xz_input_grid   !< Diurnal thermocline layerthickness
 type(esmf_field), public :: xtts_input_grid   !< d(xt)/d(ts)
 type(esmf_field), public :: xzts_input_grid   !< d(xz)/d(ts)
 type(esmf_field), public :: z_c_input_grid   !< Sub-layer cooling thickness
 type(esmf_field), public :: zm_input_grid   !< Oceanic mixed layer depth
 
 public :: read_input_nst_data
 public :: cleanup_input_nst_data
                                   
 contains                                    
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

#ifdef CHGRES_ALL
 if (trim(input_type) == "gaussian_nemsio" .or. trim(input_type) == "gfs_gaussian_nemsio") then

   call read_input_nst_nemsio_file(localpet)

!---------------------------------------------------------------------------
! Read nst data from these netcdf formatted fv3 files: tiled history,
! tiled warm restart, and gaussian history.
!---------------------------------------------------------------------------

 else

   call read_input_nst_netcdf_file(localpet)

 endif
#else

 call read_input_nst_netcdf_file(localpet)

#endif

 end subroutine read_input_nst_data
 
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

#ifdef CHGRES_ALL

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

#endif
 
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
 
 end module nst_input_data
