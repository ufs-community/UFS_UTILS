module sfc_input_data
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
 use esmf
 use netcdf
 use nemsio_module

 use program_setup, only          : data_dir_input_grid, &
                                    sfc_files_input_grid, &
                                    grib2_file_input_grid, &
                                    convert_nst, &
                                    orog_dir_input_grid, &
                                    orog_files_input_grid, &
                                    input_type, &
                                    get_var_cond, &
                                    geogrid_file_input_grid, &
                                    external_model, &
                                    vgfrc_from_climo, &
                                    minmax_vgfrc_from_climo, &
                                    lai_from_climo,&
                                    read_from_input
                                    
 use model_grid, only             : input_grid,        &
                                    i_input, j_input,  &
                                    ip1_input, jp1_input,  &
                                    num_tiles_input_grid
 use atm_input_data, only         : terrain_input_grid

 use utilities, only              : error_handler, &
                                    netcdf_err, &
                                    handle_grib_error, &
                                    to_upper, &
                                    check_soilt, &
                                    check_cnwat
                                    
! Fields associated with the land-surface model.

 integer, public                 :: veg_type_landice_input = 15 !< NOAH land ice option
                                                                !< defined at this veg type.
                                                                !< Default is igbp.
 real                            :: ICET_DEFAULT = 265.0    !< Default value of soil and skin
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
 
 public :: read_input_sfc_data
 public :: cleanup_input_sfc_data
 public :: init_sfc_esmf_fields
 
 contains
 
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

   use mpi_f08
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
         call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var= dummy2d)
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

!      print*,'lai ', maxval(gfld%fld),minval(gfld%fld)
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
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc,var8=dummy2d_8)
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
        call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8=dummy2d_8)
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
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8=dummy2d_8)
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
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8=dummy2d_8)
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
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8=dummy2d_8)
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
       call check_cnwat(dummy2d_8,i_input,j_input)
     else
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8=dummy2d_8)
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
       call handle_grib_error(vname, slev ,method,value,varnum,read_from_input,rc, var8= dummy2d_8)
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
     call check_soilt(dummy3d,slmsk_save,tsk_save,ICET_DEFAULT,i_input,j_input,lsoil_input)
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
     call handle_grib_error(vname_file, slevs(i),method,value,varnum,read_from_input,rc,var=dummy2d)
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
 
 end module sfc_input_data
