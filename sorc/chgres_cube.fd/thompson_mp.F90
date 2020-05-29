 module thompson_mp

 use esmf
 use netcdf

 implicit none

 private

 integer :: i_thomp_mp_climo = 288
 integer :: j_thomp_mp_climo = 181
 integer, public :: k_thomp_mp_climo = 30

 type(esmf_grid),  public               :: thomp_mp_climo_grid

 type(esmf_field), public :: qnifa, thomp_press

 public :: read_thomp_mp_data

 contains


 subroutine read_thomp_mp_data

 implicit none

 character(len=150) :: thomp_mp_climo_file
 character(len=3)   :: month(12)
 character(len=14)  :: record
 character(len=2)   :: level

 integer            :: error, ncid, rc, clb(2), cub(2)
 integer            :: i, j, k, localpet, npets, id_var

 real(esmf_kind_r8), allocatable :: dummy3d(:,:,:)

 real(esmf_kind_r8), pointer :: lat_ptr(:,:), lon_ptr(:,:)
 real(esmf_kind_r8) :: lon, lat

 type(esmf_vm)                :: vm

 type(esmf_polekind_flag)         :: polekindflag(2)

 data month /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
             'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /

! First create esmf grid object for the climo grid.

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", rc)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", rc)

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 print*,"- CALL GridCreate1PeriDim FOR THOMP_MP CLIMO GRID."
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

 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(thomp_mp_climo_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", rc)

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
   lon = -180.0 + ( float(i-1) * (360.0/288.0) )
   do j = clb(2), cub(2)
     lon_ptr(i,j) = lon
   enddo
 enddo

 do j = clb(2), cub(2)
   lat = -90.0 + ( float(j-1) * 1.0 )
   do i = clb(1), cub(1)
     lat_ptr(i,j) = lat
   enddo
 enddo

 if(localpet == 1) print*,'lon check ',lon_ptr(:,clb(2))
 if(localpet == 0) print*,'lat check ',lat_ptr(clb(1),:)

 print*,"- CALL FieldCreate FOR QNIFA."
 qnifa = ESMF_FieldCreate(thomp_mp_climo_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/k_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)

 print*,"- CALL FieldCreate FOR THOMP PRESS"
 thomp_press = ESMF_FieldCreate(thomp_mp_climo_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   ungriddedLBound=(/1/), &
                                   ungriddedUBound=(/k_thomp_mp_climo/), rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", rc)


 thomp_mp_climo_file = "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/chgres_thomp_mp/QNWFA_QNIFA_SIGMA_MONTHLY.dat.nc"


 print*,"- READ THOMP_MP_CLIMO_FILE: ", trim(thomp_mp_climo_file)
 error=nf90_open(trim(thomp_mp_climo_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(thomp_mp_climo_file) )
 
 
 if (localpet == 0) then
   allocate(dummy3d(i_thomp_mp_climo, j_thomp_mp_climo, k_thomp_mp_climo))
   dummy3d = 0.0
 else
   allocate(dummy3d(0,0,0))
 endif

 if (localpet == 0) then
   do k = 1, k_thomp_mp_climo
     write(level, "(I2.2)") k
     if (k < 10) then
     record = "QNIFA_" // month(1) // "__" // level
     else
     record = "QNIFA_" // month(1) // "__0" // level
     endif
     print*,'record is ',record
     error=nf90_inq_varid(ncid, trim(record), id_var)
     call netcdf_err(error, 'reading  field id' )
     error=nf90_get_var(ncid, id_var, dummy3d(:,:,k))
     call netcdf_err(error, 'reading field' )
     print*,'maxval ', record,  maxval(dummy3d(:,:,k)),minval(dummy3d(:,:,k))
   enddo
 endif

 print*,"- CALL FieldScatter FOR qnifa."
 call ESMF_FieldScatter(qnifa, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)

 if (localpet == 0) then
   do k = 1, k_thomp_mp_climo
     write(level, "(I2.2)") k
     if (k < 10) then
     record = "P_WIF_" // month(1) // "__" // level
     else
     record = "P_WIF_" // month(1) // "__0" // level
     endif
     print*,'record is ',record
     error=nf90_inq_varid(ncid, trim(record), id_var)
     call netcdf_err(error, 'reading  field id' )
     error=nf90_get_var(ncid, id_var, dummy3d(:,:,k))
     call netcdf_err(error, 'reading field' )
     print*,'maxval ', record,  maxval(dummy3d(:,:,k)),minval(dummy3d(:,:,k))
   enddo
 endif

 print*,"- CALL FieldScatter FOR thomp press."
 call ESMF_FieldScatter(thomp_press, dummy3d, rootpet=0, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldScatter", rc)
 error=nf90_close(ncid)

 end subroutine read_thomp_mp_data





 end module thompson_mp
