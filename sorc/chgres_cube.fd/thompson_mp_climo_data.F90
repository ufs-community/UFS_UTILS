!> @file
!! @brief Process Thompson climatological MP data.
!! @author George Gayno NOAA/EMC

!> Module to read the Thompson climatological MP data file
!! and set up the associated esmf field and grid objects.
!!
!! @author George Gayno NOAA/EMC
 module thompson_mp_climo_data

 use esmf
 use netcdf
 use program_setup, only      : cycle_mon, cycle_day, cycle_hour, &
                                thomp_mp_climo_file

 implicit none

 private

 integer                    :: i_thomp_mp_climo
                               !< i-dimension of Thompson climo data
 integer                    :: j_thomp_mp_climo
                               !< j-dimension of Thompson climo data
 integer, public            :: lev_thomp_mp_climo
                               !< number of vert lvls of Thompson climo data

 type(esmf_grid)            :: thomp_mp_climo_grid
                               !< esmf grid object for Thompson data grid

 type(esmf_field), public   :: qnifa_climo_input_grid
                               !< number concentration of ice friendly
                               !! nuclei.
 type(esmf_field), public   :: qnwfa_climo_input_grid
                               !< number concentration of water friendly
                               !! nuclei.
 type(esmf_field), public   :: thomp_pres_climo_input_grid
                               !< 3-d pressure of the Thompson climo
                               !! data points

 public                     :: read_thomp_mp_climo_data
 public                     :: cleanup_thomp_mp_climo_input_data

 contains

!> Read Thompson climatological MP data file and time interpolate data
!! to current cycle time. 
!!
!! @author George Gayno NOAA/EMC
 subroutine read_thomp_mp_climo_data

 implicit none

 integer            :: error, ncid, rc, clb(2), cub(2)
 integer            :: i, j, localpet, npets, id_var
 integer            :: jda(8), jdow, jdoy, jday, id_dim
 integer            :: mm, mmm, mmp, mon1, mon2

 real(esmf_kind_r8), allocatable :: dummy3d(:,:,:)
 real(esmf_kind_r8), allocatable :: dummy3d_mon1(:,:,:)
 real(esmf_kind_r8), allocatable :: dummy3d_mon2(:,:,:)
 real(esmf_kind_r8), pointer     :: lat_ptr(:,:), lon_ptr(:,:)
 real(esmf_kind_r8), allocatable :: lons(:), lats(:)
 real                            :: rjday, dayhf(13), wei1m, wei2m

 type(esmf_vm)                   :: vm

 type(esmf_polekind_flag)        :: polekindflag(2)

 data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,   &
             196.5,227.5,258.0,288.5,319.0,349.5,380.5/

!-----------------------------------------------------------------------------------
! Open the file and read the grid dimensions and latitude/longitude.
!-----------------------------------------------------------------------------------

 print*,"- READ THOMP_MP_CLIMO_FILE: ", trim(thomp_mp_climo_file)
 error=nf90_open(trim(thomp_mp_climo_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(thomp_mp_climo_file) )

 error=nf90_inq_dimid(ncid, 'lat', id_dim)
 call netcdf_err(error, 'reading lat id')
 error=nf90_inquire_dimension(ncid,id_dim,len=j_thomp_mp_climo)
 call netcdf_err(error, 'reading lat')

 error=nf90_inq_dimid(ncid, 'lon', id_dim)
 call netcdf_err(error, 'reading lon id')
 error=nf90_inquire_dimension(ncid,id_dim,len=i_thomp_mp_climo)
 call netcdf_err(error, 'reading lon')

 error=nf90_inq_dimid(ncid, 'plev', id_dim)
 call netcdf_err(error, 'reading plev id')
 error=nf90_inquire_dimension(ncid,id_dim,len=lev_thomp_mp_climo)
 call netcdf_err(error, 'reading plev')

 allocate(lons(i_thomp_mp_climo))
 allocate(lats(j_thomp_mp_climo))
 error=nf90_inq_varid(ncid, 'lon', id_var)
 call netcdf_err(error, 'reading lon field id' )
 error=nf90_get_var(ncid, id_var, lons)
 call netcdf_err(error, 'reading grid longitude' )
 error=nf90_inq_varid(ncid, 'lat', id_var)
 call netcdf_err(error, 'reading lat field id' )
 error=nf90_get_var(ncid, id_var, lats)
 call netcdf_err(error, 'reading grid latitude' )

!-----------------------------------------------------------------------------------
! Now that we have the grid information, create the esmf grid object.
!-----------------------------------------------------------------------------------

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", rc)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", rc)

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 print*,"- CALL GridCreate1PeriDim FOR THOMP MP CLIMO GRID."
 thomp_mp_climo_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_thomp_mp_climo,j_thomp_mp_climo/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN GridCreate1PeriDim", rc)

 print*,"- CALL GridAddCoord FOR THOMP MP CLIMO GRID."
 call ESMF_GridAddCoord(thomp_mp_climo_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

!-----------------------------------------------------------------------------------
! Set the grid object lat/lon.
!-----------------------------------------------------------------------------------

 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 nullify(lon_ptr)
 call ESMF_GridGetCoord(thomp_mp_climo_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 nullify(lat_ptr)
 call ESMF_GridGetCoord(thomp_mp_climo_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_ptr, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", rc)

 do i = clb(1), cub(1)
     lon_ptr(i,:) = lons(i)
 enddo

 do j = clb(2), cub(2)
     lat_ptr(:,j) = lats(j)
 enddo

!-----------------------------------------------------------------------------------
! Create esmf fields for the two tracers and 3-d pressure.
!-----------------------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR QNIFA INPUT CLIMO."
 qnifa_climo_input_grid = ESMF_FieldCreate(thomp_mp_climo_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR QNWFA INPUT CLIMO."
 qnwfa_climo_input_grid = ESMF_FieldCreate(thomp_mp_climo_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR THOMP PRESS CLIMO."
 thomp_pres_climo_input_grid = ESMF_FieldCreate(thomp_mp_climo_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/lev_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

!-----------------------------------------------------------------------------------
! Data are monthly and valid at the 15th of the month.  Compute time interpolation
! weights for the current cycle.
!-----------------------------------------------------------------------------------

 jda=0
 jda(1) = 2007
 if (cycle_mon == 2 .and. cycle_day == 29) then  ! leap year
   jda(2) = 3
   jda(3) = 1
 else
   jda(2) = cycle_mon
   jda(3) = cycle_day
 endif

 jda(5) = cycle_hour
 jdow = 0
 jdoy = 0
 jday = 0
 call w3doxdat(jda,jdow,jdoy,jday)
 rjday = float(jdoy) + float(jda(5)) / 24.
 if(rjday < dayhf(1)) rjday = rjday + 365.

 do mm=1,12
   mmm = mm
   mmp = mm + 1
   if(rjday >= dayhf(mmm) .and. rjday < dayhf(mmp)) then
     mon1 = mmm
     mon2 = mmp
     exit
   endif
 enddo

 wei1m = (dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
 wei2m = 1.0 - wei1m

 if (mon2==13) mon2=1

 print*,"- BOUNDING MONTHS AND INTERPOLATION WEIGHTS: ", mon1, wei1m, mon2, wei2m

!-----------------------------------------------------------------------------------
! Read tracers and 3-d pressure for each bounding month.  Then linearly
! interpolate in time.
!-----------------------------------------------------------------------------------

 if (localpet == 0) then
   allocate(dummy3d(i_thomp_mp_climo, j_thomp_mp_climo, lev_thomp_mp_climo))
   dummy3d = 0.0
   allocate(dummy3d_mon1(i_thomp_mp_climo, j_thomp_mp_climo, lev_thomp_mp_climo))
   dummy3d_mon1 = 0.0
   allocate(dummy3d_mon2(i_thomp_mp_climo, j_thomp_mp_climo, lev_thomp_mp_climo))
   dummy3d_mon2 = 0.0
 else
   allocate(dummy3d(0,0,0))
   allocate(dummy3d_mon1(0,0,0))
   allocate(dummy3d_mon2(0,0,0))
 endif

 if (localpet == 0) then
   print*,"- READ QNIFA FOR BOUNDING MONTH 1"
   error=nf90_inq_varid(ncid, 'nifa', id_var)
   call netcdf_err(error, 'reading nifa field id' )
   error=nf90_get_var(ncid, id_var, dummy3d_mon1, start=(/1,1,1,mon1/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading nifa month1 field' )
   print*,"- READ QNIFA FOR BOUNDING MONTH 2"
   error=nf90_get_var(ncid, id_var, dummy3d_mon2, start=(/1,1,1,mon2/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading nifa month2 field' )
   dummy3d(:,:,:) = wei1m * dummy3d_mon1 + wei2m * dummy3d_mon2
 endif

 print*,"- CALL FieldScatter FOR qnifa input climo."
 call ESMF_FieldScatter(qnifa_climo_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ QNWFA FOR BOUNDING MONTH 1"
   error=nf90_inq_varid(ncid, 'nwfa', id_var)
   call netcdf_err(error, 'reading nwfa field id' )
   error=nf90_get_var(ncid, id_var, dummy3d_mon1, start=(/1,1,1,mon1/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading nwfa month1 field' )
   print*,"- READ QNWFA FOR BOUNDING MONTH 2"
   error=nf90_get_var(ncid, id_var, dummy3d_mon2, start=(/1,1,1,mon2/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading nwfa month2 field' )
   dummy3d(:,:,:) = wei1m * dummy3d_mon1 + wei2m * dummy3d_mon2
 endif

 print*,"- CALL FieldScatter FOR qnwfa input climo."
 call ESMF_FieldScatter(qnwfa_climo_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   print*,"- READ PRESSURE FOR BOUNDING MONTH 1"
   error=nf90_inq_varid(ncid, 'prs', id_var)
   call netcdf_err(error, 'reading prs field id' )
   error=nf90_get_var(ncid, id_var, dummy3d_mon1, start=(/1,1,1,mon1/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading prs month1 field' )
   print*,"- READ PRESSURE FOR BOUNDING MONTH 2"
   error=nf90_get_var(ncid, id_var, dummy3d_mon2, start=(/1,1,1,mon2/), &
           count=(/i_thomp_mp_climo,j_thomp_mp_climo,lev_thomp_mp_climo,1/) )
   call netcdf_err(error, 'reading prs month2 field' )
   dummy3d(:,:,:) = wei1m * dummy3d_mon1 + wei2m * dummy3d_mon2
 endif

 print*,"- CALL FieldScatter FOR thomp press."
 call ESMF_FieldScatter(thomp_pres_climo_input_grid, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 error=nf90_close(ncid)

 deallocate(lons, lats, dummy3d, dummy3d_mon1, dummy3d_mon2)

 end subroutine read_thomp_mp_climo_data

!> Free up memory associated with this module.
!!
!! @author George Gayno NOAA/EMC
 subroutine cleanup_thomp_mp_climo_input_data

 implicit none

 integer                    :: rc

 call ESMF_GridDestroy(thomp_mp_climo_grid, rc=rc)
 call ESMF_FieldDestroy(thomp_pres_climo_input_grid, rc=rc)
 call ESMF_FieldDestroy(qnifa_climo_input_grid, rc=rc)
 call ESMF_FieldDestroy(qnwfa_climo_input_grid, rc=rc)

 end subroutine cleanup_thomp_mp_climo_input_data

 end module thompson_mp_climo_data
