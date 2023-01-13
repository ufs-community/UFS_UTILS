!> @file
!! @brief Read atmospheric data from GRIB2, NEMSIO and NetCDF files.
!! @author George Gayno NCEP/EMC

!> Read atmospheric data on the input grid.
!! Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!! 'history' files, fv3 gaussian history files, spectral gfs
!! gaussian nemsio files, and spectral gfs sigio/sfcio files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author George Gayno NCEP/EMC

module atm_input_data
 use esmf
 use netcdf
 use nemsio_module

 use program_setup, only          : data_dir_input_grid, &
                                    atm_files_input_grid, &
                                    grib2_file_input_grid, &
                                    atm_core_files_input_grid, &
                                    atm_tracer_files_input_grid, &
                                    tracers_input, num_tracers_input, &
                                    tracers, &
                                    get_var_cond, &
                                    external_model, &
                                    read_from_input, &
                                    input_type
 use model_grid, only             : input_grid,        &
                                    i_input, j_input,  &
                                    ip1_input, jp1_input,  &
                                    num_tiles_input_grid, &
                                    latitude_input_grid, &
                                    longitude_input_grid
 use utilities, only              : error_handler, &
                                    netcdf_err, &
                                    handle_grib_error, &
                                    quicksort, &
                                    dint2p                                    
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
!type(esmf_field), public              :: wind_input_grid       !< 3-component wind
 type(esmf_field), public              :: xwind_input_grid       !< 3-component wind
 type(esmf_field), public              :: ywind_input_grid       !< 3-component wind
 type(esmf_field), public              :: zwind_input_grid       !< 3-component wind
 type(esmf_field), allocatable, public :: tracers_input_grid(:) !< tracers

 integer, public                 :: lev_input      !< number of atmospheric layers
 integer, public                 :: levp1_input    !< number of atmos layer interfaces

 character(len=50), private, allocatable :: slevs(:) !< The atmospheric levels in the GRIB2 input file.

 public :: read_input_atm_data
 public :: cleanup_input_atm_data
 public :: convert_winds
 
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

!print*,"- CALL FieldCreate FOR INPUT GRID 3-D WIND."
!wind_input_grid = ESMF_FieldCreate(input_grid, &
!                                  typekind=ESMF_TYPEKIND_R8, &
!                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
!                                  ungriddedLBound=(/1,1/), &
!                                  ungriddedUBound=(/lev_input,3/), rc=rc)
!if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!   call error_handler("IN FieldCreate", rc)

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

 print*,"- CALL FieldCreate FOR INPUT GRID xwind."
 xwind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID ywind."
 ywind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR INPUT GRID zwind."
 zwind_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_input/), rc=rc)
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

          call handle_grib_error(vname, slevs(vlev),method,value,varnum,read_from_input,iret,var=dummy2d)
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
         call handle_grib_error(vname, slevs(vlev),method,value,varnum,read_from_input,iret,var8=dum2d_1)
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
        call handle_grib_error(vname, slevs(vlev),method_u,value_u,varnum_u,read_from_input,iret,var=u_tmp)
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
        call handle_grib_error(vname, slevs(vlev),method_v,value_v,varnum_v,read_from_input,iret,var=v_tmp)
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

 integer                         :: clb(3), cub(3)
 integer                         :: i, j, k, rc

 real(esmf_kind_r8)              :: latrad, lonrad
!real(esmf_kind_r8), pointer     :: windptr(:,:,:,:)
 real(esmf_kind_r8), pointer     :: xptr(:,:,:)
 real(esmf_kind_r8), pointer     :: yptr(:,:,:)
 real(esmf_kind_r8), pointer     :: zptr(:,:,:)
 real(esmf_kind_r8), pointer     :: uptr(:,:,:)
 real(esmf_kind_r8), pointer     :: vptr(:,:,:)
 real(esmf_kind_r8), pointer     :: latptr(:,:)
 real(esmf_kind_r8), pointer     :: lonptr(:,:)

!print*,"- CALL FieldGet FOR 3-D WIND."
!call ESMF_FieldGet(wind_input_grid, &
!                   computationalLBound=clb, &
!                   computationalUBound=cub, &
!                   farrayPtr=windptr, rc=rc)
!if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!   call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR x."
 call ESMF_FieldGet(xwind_input_grid, &
                    computationalLBound=clb, &
                    computationalUBound=cub, &
                    farrayPtr=xptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR y."
 call ESMF_FieldGet(ywind_input_grid, &
                    farrayPtr=yptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", rc)

 print*,"- CALL FieldGet FOR z."
 call ESMF_FieldGet(zwind_input_grid, &
                    farrayPtr=zptr, rc=rc)
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
!      windptr(i,j,k,1) = uptr(i,j,k) * cos(lonrad) - vptr(i,j,k) * sin(latrad) * sin(lonrad)
       xptr(i,j,k) = uptr(i,j,k) * cos(lonrad) - vptr(i,j,k) * sin(latrad) * sin(lonrad)
!      windptr(i,j,k,2) = uptr(i,j,k) * sin(lonrad) + vptr(i,j,k) * sin(latrad) * cos(lonrad)
       yptr(i,j,k) = uptr(i,j,k) * sin(lonrad) + vptr(i,j,k) * sin(latrad) * cos(lonrad)
!      windptr(i,j,k,3) = vptr(i,j,k) * cos(latrad)
       zptr(i,j,k) = vptr(i,j,k) * cos(latrad)
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
  real(esmf_kind_r4)                  :: dtor = 3.14159265359_esmf_kind_r4/180.0_esmf_kind_r4
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
!call ESMF_FieldDestroy(wind_input_grid, rc=rc)
 call ESMF_FieldDestroy(xwind_input_grid, rc=rc)
 call ESMF_FieldDestroy(ywind_input_grid, rc=rc)
 call ESMF_FieldDestroy(zwind_input_grid, rc=rc)
 call ESMF_FieldDestroy(ps_input_grid, rc=rc)

 do n = 1, num_tracers_input
   call ESMF_FieldDestroy(tracers_input_grid(n), rc=rc)
 enddo
 deallocate(tracers_input_grid)

 end subroutine cleanup_input_atm_data
 
end module atm_input_data
