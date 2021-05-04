!> @file
!! @brief Read atmospheric and surface data from GRIB2, NEMSIO and NetCDF files.
!! @author George Gayno, Kyle Gerheiser NCEP/EMC

!> Read atmospheric data on the input grid.
!! Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!! 'history' files, fv3 gaussian history files, spectral gfs
!! gaussian nemsio files, and spectral gfs sigio/sfcio files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author George Gayno, Kyle Gerheiser NCEP/EMC
module atm_input_data_mod
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
       longitude_input_grid, &
       inv_file

  use input_data_utils_mod, only: handle_grib_error, terrain_input_grid
  
  implicit none


  private
  ! Fields associated with the atmospheric model.

  type(esmf_field), public              :: dzdt_input_grid       !< vert velocity
  type(esmf_field)                      :: dpres_input_grid      !< pressure thickness
  type(esmf_field), public              :: pres_input_grid       !< 3-d pressure
  type(esmf_field), public              :: ps_input_grid         !< surface pressure
  type(esmf_field), public              :: temp_input_grid       !< temperature
  type(esmf_field)                      :: u_input_grid          !< u/v wind at grid
  type(esmf_field)                      :: v_input_grid          !< box center
  type(esmf_field), public              :: wind_input_grid       !< 3-component wind
  type(esmf_field), allocatable, public :: tracers_input_grid(:) !< tracers

  integer, public                 :: lev_input      !< number of atmospheric layers
  integer, public                 :: levp1_input    !< number of atmos layer interfaces

  character(len=50), private, allocatable :: slevs(:) !< The atmospheric levels in the GRIB2 input file.

  public :: read_input_atm_data
  public :: cleanup_input_atm_data

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

  !> Read input grid atmospheric fv3gfs grib2 files.
  !!
  !! @param[in] localpet  ESMF local persistent execution thread 
  !! @author George Gayno NCEP/EMC   
  subroutine read_input_atm_grib2_file(localpet)

    use wgrib2api

    use grib2_util, only                   : rh2spfh, convert_omega
    use input_data_utils_mod, only: quicksort

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
         hasspfh=.true., &
         isnative=.false.

    real(esmf_kind_r8), allocatable       :: rlevs(:)
    real(esmf_kind_r4), allocatable       :: dummy2d(:,:)
    real(esmf_kind_r8), allocatable       :: dummy3d(:,:,:), dummy2d_8(:,:),&
         u_tmp_3d(:,:,:), v_tmp_3d(:,:,:)
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
    iret = grb2_inq(the_file,inv_file,":var0_2","_0_0:",":10 hybrid level:")

    if (iret <= 0) then
       lvl_str = "mb:" 
       lvl_str_space = " mb:"
       lvl_str_space_len = 4
       isnative = .false.
       iret = grb2_inq(the_file,inv_file,":UGRD:",lvl_str_space)
       lev_input=iret
       if (localpet == 0) print*,"- DATA IS ON ", lev_input, " ISOBARIC LEVELS."
    else
       lvl_str = " level:"
       lvl_str_space = " hybrid "
       lvl_str_space_len = 7
       isnative = .true.
       iret = grb2_inq(the_file,inv_file,":UGRD:",lvl_str_space, " level:")
       if (iret < 0) call error_handler("READING VERTICAL LEVEL TYPE.", iret)
       lev_input=iret
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
       if (.not. isnative) rlevs(i) = rlevs(i) * 100.0
       if (localpet==0) print*, "- LEVEL = ", slevs(i)
    enddo

    ! Jili Dong add sort to re-order isobaric levels.

    call quicksort(rlevs,1,lev_input)

    if (.not. isnative) then
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
    endif

    if (localpet == 0) print*,"- FIND SPFH OR RH IN FILE"
    iret = grb2_inq(the_file,inv_file,trim(trac_names_grib_1(1)),trac_names_grib_2(1),lvl_str_space)

    if (iret <= 0) then
       iret = grb2_inq(the_file,inv_file, ':var0_2','_1_1:',lvl_str_space)
       if (iret <= 0) call error_handler("READING ATMOSPHERIC WATER VAPOR VARIABLE.", iret)
       hasspfh = .false.
       trac_names_grib_2(1)='_1_1:'
       if (localpet == 0) print*,"- FILE CONTAINS RH."
    else
       if (localpet == 0) print*,"- FILE CONTAINS SPFH."
    endif

    if (localpet == 0) print*,"- FIND ICMR, SCLIWC, OR CICE IN FILE"
    iret = grb2_inq(the_file,inv_file,trac_names_grib_1(4),trac_names_grib_2(4),lvl_str_space)

    if (iret <= 0) then
       vname = trac_names_vmap(4)
       print*, "vname = ", vname
       call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
            this_field_var_name=tmpstr,loc=varnum)
       iret = grb2_inq(the_file,inv_file, ':var0_2','_1_84:',lvl_str_space)
       if (iret <= 0) then
          iret = grb2_inq(the_file,inv_file, ':var0_2','_6_0:',lvl_str_space)
          if (iret <= 0 ) then 
             call handle_grib_error(vname, slevs(1),method,value,varnum,rc,var=dummy2d)
          else
             trac_names_grib_2(4) = '_6_0'
             if (localpet == 0) print*,"- FILE CONTAINS CICE."
          endif
       else
          trac_names_grib_2(4)='_1_84:'
          if (localpet == 0) print*,"- FILE CONTAINS SCLIWC."
       endif
    else
       if (localpet == 0) print*,"- FILE CONTAINS ICMR."
    endif

    if (localpet == 0) print*,"- FIND CLWMR or SCLLWC IN FILE"
    iret = grb2_inq(the_file,inv_file,trac_names_grib_1(5),trac_names_grib_2(5),lvl_str_space)

    if (iret <= 0) then
       vname = trac_names_vmap(5)
       print*, "vname = ", vname
       call get_var_cond(vname,this_miss_var_method=method, this_miss_var_value=value, &
            this_field_var_name=tmpstr,loc=varnum)
       iret = grb2_inq(the_file,inv_file, ':var0_2','_1_83:',lvl_str_space)
       if (iret <= 0) then 
          call handle_grib_error(vname, slevs(1),method,value,varnum,rc,var=dummy2d)
       elseif (iret <=0 .and. rc .ne. 1) then
          call error_handler("READING CLOUD WATER VARIABLE.", iret)
       else
          trac_names_grib_2(4)='_1_83:'
          if (localpet == 0) print*,"- FILE CONTAINS SCLLWC."
       endif
    else
       if (localpet == 0) print*,"- FILE CONTAINS CLWMR."
    endif

    do n = 1, num_tracers

       vname = tracers_input(n)

       i = maxloc(merge(1.,0.,trac_names_vmap == vname),dim=1)

       tracers_input_grib_1(n) = trac_names_grib_1(i)
       tracers_input_grib_2(n) = trac_names_grib_2(i)
       tracers_input_vmap(n)=trac_names_vmap(i)
       tracers(n)=tracers_default(i)

    enddo

    if (localpet==0) print*, "- NUMBER OF TRACERS TO BE PROCESSED = ", num_tracers

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

    !----------------------------------------------------------------------------------
    ! This program expects field levels from bottom to top. Fields in non-native 
    ! files read in from top to bottom. We will flip indices later. Fields on 
    ! native vertical coordinates read from bottom to top so those need no adjustments.
    !----------------------------------------------------------------------------------

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

    call read_winds(the_file,inv_file,u_tmp_3d,v_tmp_3d, localpet)

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
       vname = ":var0_2"
       vname2 = "_2_9:"
       do vlev = 1, lev_input
          iret = grb2_inq(the_file,inv_file,vname,vname2,slevs(vlev),data2=dummy2d)
          if (iret <= 0 ) then
             print*,"DZDT not available at level ", trim(slevs(vlev)), " so checking for VVEL"
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

    else
       ! For native files, read in pressure field directly from file but don't flip levels
       if (localpet == 0) then
          print*,"- READ PRESSURE."
          vname = ":PRES:"
          do vlev = 1, lev_input
             iret = grb2_inq(the_file,inv_file,vname,slevs(vlev),data2=dummy2d)
             if (iret<=0) then
                call error_handler("READING IN PRESSURE AT LEVEL "//trim(slevs(vlev)),iret)
             endif
             dummy3d(:,:,vlev) = real(dummy2d,esmf_kind_r8)
             print*,'pres check after read ',vlev, dummy3d(1,1,vlev)
          enddo
       endif

       if (localpet == 0) print*,"- CALL FieldScatter FOR INPUT GRID PRESSURE."
       call ESMF_FieldScatter(pres_input_grid, dummy3d, rootpet=0, rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldScatter", rc)
    endif
    deallocate(dummy3d) 

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


  !> Read winds from a grib2 file.  Rotate winds
  !! to be earth relative if necessary.
  !!
  !! @param [in] file  grib2 file to be read
  !! @param [in] inv   grib2 inventory file
  !! @param [inout] u  u-component wind
  !! @param [inout] v  v-component wind
  !! @param[in] localpet  ESMF local persistent execution thread
  !! @author Larissa Reames
  subroutine read_winds(file,inv,u,v,localpet)

    use wgrib2api
    use netcdf
    use program_setup, only      : get_var_cond, fix_dir_input_grid
    use model_grid, only         : input_grid_type
    implicit none

    character(len=250), intent(in)          :: file
    character(len=10), intent(in)            :: inv
    integer, intent(in)                     :: localpet
    real(esmf_kind_r8), intent(inout), allocatable :: u(:,:,:),v(:,:,:)

    real(esmf_kind_r4), dimension(i_input,j_input)  :: alpha
    real(esmf_kind_r8), dimension(i_input,j_input)  :: lon, lat
    real(esmf_kind_r4), allocatable                 :: u_tmp(:,:),v_tmp(:,:)
    real(esmf_kind_r4), dimension(i_input,j_input)  :: ws,wd
    real(esmf_kind_r4)                      :: value_u, value_v,lov,latin1,latin2
    real(esmf_kind_r8)                      :: d2r

    integer                                 :: varnum_u, varnum_v, vlev, & !ncid, id_var, &
         error, iret, istr

    character(len=20)                       :: vname
    character(len=50)                       :: method_u, method_v
    character(len=250)                      :: file_coord
    character(len=10000)                    :: temp_msg

    d2r=acos(-1.0_esmf_kind_r8) / 180.0_esmf_kind_r8
    if (localpet==0) then
       allocate(u(i_input,j_input,lev_input))
       allocate(v(i_input,j_input,lev_input))
    else
       allocate(u(0,0,0))
       allocate(v(0,0,0))
    endif

    file_coord = trim(fix_dir_input_grid)//"/latlon_grid3.32769.nc"

    vname = "u"
    call get_var_cond(vname,this_miss_var_method=method_u, this_miss_var_value=value_u, &
         loc=varnum_u)
    vname = "v"
    call get_var_cond(vname,this_miss_var_method=method_v, this_miss_var_value=value_v, &
         loc=varnum_v)

    if (trim(input_grid_type)=="rotated_latlon") then
       print*,"- CALL FieldGather FOR INPUT GRID LONGITUDE"
       call ESMF_FieldGather(longitude_input_grid, lon, rootPet=0, tile=1, rc=error)
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", error)
       print*,"- CALL FieldGather FOR INPUT GRID LATITUDE"
       call ESMF_FieldGather(latitude_input_grid, lat, rootPet=0, tile=1, rc=error)
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", error)

       if (localpet==0) then
          print*,"- CALCULATE ROTATION ANGLE FOR ROTATED_LATLON INPUT GRID"
          error = grb2_inq(file, inv,grid_desc=temp_msg)
          !1:0:grid_template=32769:winds(grid):
          !   I am not an Arakawa E-grid.
          !   I am rotated but have no rotation angle.
          !   I am staggered. What am I?
          !   (953 x 834) units 1e-06 input WE:SN output WE:SN res 56
          !   lat0 -10.590603 lat-center 54.000000 dlat 121.813000
          !   lon0 220.914154 lon-center 254.000000 dlon 121.813000 #points=794802

          istr = index(temp_msg, "lat-center ") + len("lat_center ")
          read(temp_msg(istr:istr+9),"(F8.5)") latin1
          istr = index(temp_msg, "lon-center ") + len("lon-center ")
          read(temp_msg(istr:istr+10),"(F9.6)") lov

          print*, "- CALL CALCALPHA_ROTLATLON with center lat,lon = ",latin1,lov
          call calcalpha_rotlatlon(lat,lon,latin1,lov,alpha)
          print*, " alpha min/max = ",MINVAL(alpha),MAXVAL(alpha)
       endif
    elseif (trim(input_grid_type) == "lambert") then
       !# NG this has been edited to correctly calculate gridrot for Lambert grids
       !  Previously was incorrectly using polar-stereographic formation
       print*,"- CALL FieldGather FOR INPUT GRID LONGITUDE"
       call ESMF_FieldGather(longitude_input_grid, lon, rootPet=0, tile=1, rc=error)
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", error)

       if (localpet==0) then
          error = grb2_inq(file, inv,grid_desc=temp_msg)
          !1:0:grid_template=30:winds(grid):
          !   Lambert Conformal: (1799 x 1059) input WE:SN output WE:SN res 8
          !   Lat1 21.138123 Lon1 237.280472 LoV 262.500000
          !   LatD 38.500000 Latin1 38.500000 Latin2 38.500000
          !   LatSP 0.000000 LonSP 0.000000
          !   North Pole (1799 x 1059) Dx 3000.000000 m Dy 3000.000000 m mode 8

          istr = index(temp_msg, "LoV ") + len("LoV ")
          read(temp_msg(istr:istr+10),"(F9.6)") lov
          istr = index(temp_msg, "Latin1 ") + len("Latin1 ")
          read(temp_msg(istr:istr+9),"(F8.5)") latin1
          istr = index(temp_msg, "Latin2 ") + len("Latin2 ")
          read(temp_msg(istr:istr+9),"(F8.5)") latin2

          print*, "- CALL GRIDROT for LC grid with lov,latin1/2 = ",lov,latin1,latin2
          call gridrot(lov,latin1,latin2,lon,alpha)
          print*, " alpha min/max = ",MINVAL(alpha),MAXVAL(alpha)
       endif
    endif

    if (localpet==0) then
       do vlev = 1, lev_input

          vname = ":UGRD:"
          iret = grb2_inq(file,inv,vname,slevs(vlev),data2=u_tmp)
          if (iret <= 0) then
             call handle_grib_error(vname, slevs(vlev),method_u,value_u,varnum_u,iret,var=u_tmp)
             if (iret==1) then ! missing_var_method == skip
                call error_handler("READING IN U AT LEVEL "//trim(slevs(vlev))//". SET A FILL "// &
                     "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
             endif
          endif

          vname = ":VGRD:"
          iret = grb2_inq(file,inv,vname,slevs(vlev),data2=v_tmp)
          if (iret <= 0) then
             call handle_grib_error(vname, slevs(vlev),method_v,value_v,varnum_v,iret,var=v_tmp)
             if (iret==1) then ! missing_var_method == skip
                call error_handler("READING IN V AT LEVEL "//trim(slevs(vlev))//". SET A FILL "// &
                     "VALUE IN THE VARMAP TABLE IF THIS ERROR IS NOT DESIRABLE.",iret)
             endif
          endif

          if (trim(input_grid_type) == "latlon") then
             if (external_model == 'UKMET') then
                u(:,:,vlev) = u_tmp
                v(:,:,vlev) = (v_tmp(:,2:jp1_input) + v_tmp(:,1:j_input))/2
             else
                u(:,:,vlev) = u_tmp
                v(:,:,vlev) = v_tmp
             endif
          else if (trim(input_grid_type) == "rotated_latlon") then
             ws = sqrt(u_tmp**2 + v_tmp**2)
             wd = atan2(-u_tmp,-v_tmp) / d2r ! calculate grid-relative wind direction
             wd = wd + alpha + 180.0 ! Rotate from grid- to earth-relative direction
             wd = 270.0 - wd ! Convert from meteorological (true N) to mathematical direction
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
       an = log( cos(latin1*dtor) / cos(latin2*dtor) ) / &
            log( tan(dtor*(90.0-latin1)/2.) / tan(dtor*(90.0-latin2)/2.))
    end if

    tlon = mod(lon - lov + 180. + 3600., 360.) - 180.
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
    alpha = -asin(sinalpha)/D2R
    ! returns alpha in degrees
  end subroutine calcalpha_rotlatlon


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

    do n = 1, num_tracers
       call ESMF_FieldDestroy(tracers_input_grid(n), rc=rc)
    enddo
    deallocate(tracers_input_grid)

  end subroutine cleanup_input_atm_data


end module atm_input_data_mod

