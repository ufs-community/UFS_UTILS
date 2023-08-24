!> @file
!! @brief Contains utility routines that read and write data.
!! 
!! @author Xu Li, Hang Lei, George Gayno NOAA/EMC

!> This module contains routines that read and write data.
!! @author Xu Li, Hang Lei, George Gayno NOAA/EMC
MODULE READ_WRITE_DATA

 USE NETCDF

 PRIVATE

! DATA STRUCTURE TO HOLD ALL NSST RECORDS.

 TYPE, PUBLIC :: NSST_DATA
   REAL, ALLOCATABLE :: C_0(:)
   REAL, ALLOCATABLE :: C_D(:)
   REAL, ALLOCATABLE :: D_CONV(:)
   REAL, ALLOCATABLE :: DT_COOL(:)
   REAL, ALLOCATABLE :: IFD(:)
   REAL, ALLOCATABLE :: QRAIN(:)
   REAL, ALLOCATABLE :: TREF(:)
   REAL, ALLOCATABLE :: TFINC(:)
   REAL, ALLOCATABLE :: W_0(:)
   REAL, ALLOCATABLE :: W_D(:)
   REAL, ALLOCATABLE :: XS(:)
   REAL, ALLOCATABLE :: XT(:)
   REAL, ALLOCATABLE :: XTTS(:)
   REAL, ALLOCATABLE :: XU(:)
   REAL, ALLOCATABLE :: XV(:)
   REAL, ALLOCATABLE :: XZ(:)
   REAL, ALLOCATABLE :: XZTS(:)
   REAL, ALLOCATABLE :: Z_C(:)
   REAL, ALLOCATABLE :: ZM(:)
 END TYPE NSST_DATA

 INTEGER, PUBLIC              :: IDIM_GAUS !< 'i' dimension of GSI gaussian
                                           !! grid.
 INTEGER, PUBLIC              :: JDIM_GAUS !< 'j' dimension of GSI gaussian
                                           !! grid.
 INTEGER, ALLOCATABLE, PUBLIC :: SLMSK_GAUS(:,:) !< GSI land mask on the
                                                 !! gaussian grid.

 INTEGER, ALLOCATABLE, PUBLIC :: SOILSNOW_GAUS(:,:) !< GSI soil / snow mask for land on the
                                                    !! gaussian grid. 
                                                    !! 1 - soil, 2 - snow, 0 - not land

 REAL, ALLOCATABLE, PUBLIC    :: DTREF_GAUS(:,:) !< GSI foundation temperature
                                                 !! increment on the gaussian grid.

 REAL, ALLOCATABLE, PUBLIC    :: STC_INC_GAUS(:,:,:) !< GSI soil temperature increments
                                                     !! on the gaussian grid.

 REAL, ALLOCATABLE, PUBLIC    :: SLC_INC_GAUS(:,:,:) !< GSI soil moisture increments
                                                     !! on the gaussian grid.

 PUBLIC :: READ_DATA
 PUBLIC :: READ_GSI_DATA
 PUBLIC :: READ_LAT_LON_OROG
 PUBLIC :: WRITE_DATA
 public :: read_tf_clim_grb,get_tf_clm_dim
 public :: read_salclm_gfs_nc,get_dim_nc

 CONTAINS

   !> Update surface records - and nsst records if selected -
   !! on a single cubed-sphere tile to a pre-existing model 
   !! restart file (in netcdf).
   !! 
   !! @note The model restart files contain an additional snow field -
   !! snow cover (snocvr). That field is required for bit identical
   !! reproducability. If that record does not exist, the model will
   !! compute it as an initialization step. Because this program does not
   !! contain the snow cover algorithm, it will let the model compute it.
   !!
   !! @param[in] idim 'i' dimension of a tile.
   !! @param[in] jdim 'j' dimension of a tile.
   !! @param[in] lensfc Total number of points on a tile.
   !! @param[in] lsoil Number of soil layers.
   !! @param[in] do_nsst When true, nsst fields were processed.
   !! @param[in] nsst Data structure containing nsst fields.
   !! @param[in] slifcs Land-sea mask.
   !! @param[in] tsffcs Skin temperature.
   !! @param[in] vegfcs Vegetation greenness.
   !! @param[in] swefcs Snow water equivalent
   !! @param[in] tg3fcs Soil substrate temperature.
   !! @param[in] zorfcs Roughness length.
   !! @param[in] albfcs Snow-free albedo.
   !! @param[in] alffcs Fractional coverage for strong/weak zenith angle
   !! dependent albedo.
   !! @param[in] cnpfcs Plant canopy moisture content.
   !! @param[in] f10m log((z0+10)/z0). See model routine sfc_diff.f for
   !! details.
   !! @param[in] t2m Two-meter air temperature.
   !! @param[in] q2m Two-meter specific humidity.
   !! @param[in] vetfcs Vegetation type.
   !! @param[in] sotfcs Soil type.
   !! @param[in] ustar Friction velocity.
   !! @param[in] fmm log((z0+z1)/z0). See model routine sfc_diff.f for
   !! details.
   !! @param[in] fhh log(ztmax+z1)/ztmax).  See model routine sfc_diff.f for
   !! details.
   !! @param[in] sicfcs Sea ice concentraton.
   !! @param[in] sihfcs Sea ice depth.
   !! @param[in] sitfcs Sea ice temperature.
   !! @param[in] tprcp Precipitation.
   !! @param[in] srflag Snow/rain flag.
   !! @param[in] swdfcs Physical snow depth.
   !! @param[in] vmnfcs Minimum vegetation greenness.
   !! @param[in] vmxfcs Maximum vegetation greenness.
   !! @param[in] slpfcs Slope type.
   !! @param[in] absfcs Maximum snow albedo.
   !! @param[in] slcfcs Liquid portion of volumetric soil moisture.
   !! @param[in] smcfcs Total volumetric soil moisture.
   !! @param[in] stcfcs Soil temperature.
   !!
   !! @author George Gayno NOAA/EMC

 subroutine write_data(lensfc,idim,jdim,lsoil, &
                       do_nsst,nsst,slifcs,tsffcs,vegfcs,swefcs, &
                       tg3fcs,zorfcs,albfcs,alffcs, &
                       cnpfcs,f10m,t2m,q2m,vetfcs, &
                       sotfcs,ustar,fmm,fhh,sicfcs, &
                       sihfcs,sitfcs,tprcp,srflag,  &
                       swdfcs,vmnfcs,vmxfcs,slpfcs, &
                       absfcs,slcfcs,smcfcs,stcfcs)

 use mpi

 implicit none

 integer, intent(in)              :: lensfc, lsoil
 integer, intent(in)              :: idim, jdim

 logical, intent(in)              :: do_nsst

 real, intent(in), optional       :: slifcs(lensfc),tsffcs(lensfc)
 real, intent(in), optional       :: swefcs(lensfc),tg3fcs(lensfc)
 real, intent(in), optional       :: zorfcs(lensfc),albfcs(lensfc,4)
 real, intent(in), optional       :: alffcs(lensfc,2),cnpfcs(lensfc)
 real, intent(in), optional       :: f10m(lensfc),t2m(lensfc)
 real, intent(in), optional       :: q2m(lensfc),vegfcs(lensfc)
 real, intent(in), optional       :: vetfcs(lensfc),sotfcs(lensfc)
 real, intent(in), optional       :: ustar(lensfc),fmm(lensfc)
 real, intent(in), optional       :: fhh(lensfc), sicfcs(lensfc)
 real, intent(in), optional       :: sihfcs(lensfc), sitfcs(lensfc)
 real, intent(in), optional       :: tprcp(lensfc), srflag(lensfc)
 real, intent(in), optional       :: swdfcs(lensfc), vmnfcs(lensfc)
 real, intent(in), optional       :: vmxfcs(lensfc), slpfcs(lensfc)
 real, intent(in), optional       :: absfcs(lensfc), slcfcs(lensfc,lsoil)
 real, intent(in), optional       :: smcfcs(lensfc,lsoil), stcfcs(lensfc,lsoil)

 type(nsst_data), intent(in)      :: nsst

 integer :: dim_x, dim_y, dim_time, dims_3d(3)

 real :: dum2d(idim,jdim), dum3d(idim,jdim,lsoil)
 
 character(len=50) :: fnbgso
 character(len=3)  :: rankch

 integer           :: myrank, error, ncid, id_var

 call mpi_comm_rank(mpi_comm_world, myrank, error)

 write(rankch, '(i3.3)') (myrank+1)

 fnbgso = "./fnbgso." // rankch

 print*
 print*,"update OUTPUT SFC DATA TO: ",trim(fnbgso)

 ERROR=NF90_OPEN(TRIM(fnbgso),NF90_WRITE,NCID)
 CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(fnbgso) )

 if(present(slifcs)) then
   error=nf90_inq_varid(ncid, "slmsk", id_var)
   call netcdf_err(error, 'reading slmsk id' )
   dum2d = reshape(slifcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing slmsk record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(tsffcs)) then
   error=nf90_inq_varid(ncid, "tsea", id_var)
   call netcdf_err(error, 'reading tsea id' )
   dum2d = reshape(tsffcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing tsea record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(swefcs)) then
   error=nf90_inq_varid(ncid, "sheleg", id_var)
   call netcdf_err(error, 'reading sheleg id' )
   dum2d = reshape(swefcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing sheleg record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(tg3fcs)) then
   error=nf90_inq_varid(ncid, "tg3", id_var)
   call netcdf_err(error, 'reading tg3 id' )
   dum2d = reshape(tg3fcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing tg3 record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(zorfcs)) then
   error=nf90_inq_varid(ncid, "zorl", id_var)
   call netcdf_err(error, 'reading zorl id' )
   dum2d = reshape(zorfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing zorl record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(albfcs)) then
   error=nf90_inq_varid(ncid, "alvsf", id_var)
   call netcdf_err(error, 'reading alvsf id' )
   dum2d = reshape(albfcs(:,1), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing alvsf record' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "alvwf", id_var)
   call netcdf_err(error, 'reading alvwf id' )
   dum2d = reshape(albfcs(:,2), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing alvwf record' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "alnsf", id_var)
   call netcdf_err(error, 'reading alnsf id' )
   dum2d = reshape(albfcs(:,3), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing alnsf record' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "alnwf", id_var)
   call netcdf_err(error, 'reading alnwf id' )
   dum2d = reshape(albfcs(:,4), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing alnwf record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(alffcs)) then
   error=nf90_inq_varid(ncid, "facsf", id_var)
   call netcdf_err(error, 'reading facsf id' )
   dum2d = reshape(alffcs(:,1), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing facsf record' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "facwf", id_var)
   call netcdf_err(error, 'reading facwf id' )
   dum2d = reshape(alffcs(:,2), (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing facwf record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(vegfcs)) then
   error=nf90_inq_varid(ncid, "vfrac", id_var)
   call netcdf_err(error, 'reading vfrac id' )
   dum2d = reshape(vegfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing vegfcs record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(cnpfcs)) then
   error=nf90_inq_varid(ncid, "canopy", id_var)
   call netcdf_err(error, 'reading canopy id' )
   dum2d = reshape(cnpfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing canopy record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(f10m)) then
   error=nf90_inq_varid(ncid, "f10m", id_var)
   call netcdf_err(error, 'reading f10m id' )
   dum2d = reshape(f10m, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing f10m record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(t2m)) then
   error=nf90_inq_varid(ncid, "t2m", id_var)
   call netcdf_err(error, 'reading t2m id' )
   dum2d = reshape(t2m, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing t2m record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(q2m)) then
   error=nf90_inq_varid(ncid, "q2m", id_var)
   call netcdf_err(error, 'reading q2m id' )
   dum2d = reshape(q2m, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing q2m record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(vetfcs)) then
   error=nf90_inq_varid(ncid, "vtype", id_var)
   call netcdf_err(error, 'reading vtype id' )
   dum2d = reshape(vetfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing vtype record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(sotfcs)) then
   error=nf90_inq_varid(ncid, "stype", id_var)
   call netcdf_err(error, 'reading stype id' )
   dum2d = reshape(sotfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing stype record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(ustar)) then
   error=nf90_inq_varid(ncid, "uustar", id_var)
   call netcdf_err(error, 'reading uustar id' )
   dum2d = reshape(ustar, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing uustar record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(fmm)) then
   error=nf90_inq_varid(ncid, "ffmm", id_var)
   call netcdf_err(error, 'reading ffmm id' )
   dum2d = reshape(fmm, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing ffmm record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(fhh)) then
   error=nf90_inq_varid(ncid, "ffhh", id_var)
   call netcdf_err(error, 'reading ffhh id' )
   dum2d = reshape(fhh, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing ffhh record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(sicfcs)) then
   error=nf90_inq_varid(ncid, "fice", id_var)
   call netcdf_err(error, 'reading fice id' )
   dum2d = reshape(sicfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing fice record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(sihfcs)) then
   error=nf90_inq_varid(ncid, "hice", id_var)
   call netcdf_err(error, 'reading hice id' )
   dum2d = reshape(sihfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing hice record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(sitfcs)) then
   error=nf90_inq_varid(ncid, "tisfc", id_var)
   call netcdf_err(error, 'reading tisfc id' )
   dum2d = reshape(sitfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing tisfc record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(tprcp)) then
   error=nf90_inq_varid(ncid, "tprcp", id_var)
   call netcdf_err(error, 'reading tprcp id' )
   dum2d = reshape(tprcp, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing tprcp record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(srflag)) then
   error=nf90_inq_varid(ncid, "srflag", id_var)
   call netcdf_err(error, 'reading srflag id' )
   dum2d = reshape(srflag, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing srflag record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(swdfcs)) then
   error=nf90_inq_varid(ncid, "snwdph", id_var)
   call netcdf_err(error, 'reading snwdph id' )
   dum2d = reshape(swdfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing snwdph record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(vmnfcs)) then
   error=nf90_inq_varid(ncid, "shdmin", id_var)
   call netcdf_err(error, 'reading shdmin id' )
   dum2d = reshape(vmnfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing shdmin record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(vmxfcs)) then
   error=nf90_inq_varid(ncid, "shdmax", id_var)
   call netcdf_err(error, 'reading shdmax id' )
   dum2d = reshape(vmxfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing shdmax record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(slpfcs)) then
   error=nf90_inq_varid(ncid, "slope", id_var)
   call netcdf_err(error, 'reading slope id' )
   dum2d = reshape(slpfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing slope record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(absfcs)) then
   error=nf90_inq_varid(ncid, "snoalb", id_var)
   call netcdf_err(error, 'reading snoalb id' )
   dum2d = reshape(absfcs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'writing snoalb record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(slcfcs)) then
   error=nf90_inq_varid(ncid, "slc", id_var)
   call netcdf_err(error, 'reading slc id' )
   dum3d = reshape(slcfcs, (/idim,jdim,lsoil/))
   error = nf90_put_var( ncid, id_var, dum3d)
   call netcdf_err(error, 'writing slc record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(smcfcs)) then
   error=nf90_inq_varid(ncid, "smc", id_var)
   call netcdf_err(error, 'reading smc id' )
   dum3d = reshape(smcfcs, (/idim,jdim,lsoil/))
   error = nf90_put_var( ncid, id_var, dum3d)
   call netcdf_err(error, 'writing smc record' )
   call remove_checksum(ncid, id_var)
 endif

 if(present(stcfcs)) then
   error=nf90_inq_varid(ncid, "stc", id_var)
   call netcdf_err(error, 'reading stc id' )
   dum3d = reshape(stcfcs, (/idim,jdim,lsoil/))
   error = nf90_put_var( ncid, id_var, dum3d)
   call netcdf_err(error, 'writing stc record' )
   call remove_checksum(ncid, id_var)
 endif

 if(do_nsst) then

   error=nf90_inq_varid(ncid, "tref", id_var)
   call netcdf_err(error, 'reading tref id' )
   dum2d = reshape(nsst%tref, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING TREF RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "z_c", id_var)
   call netcdf_err(error, 'reading z_c id' )
   dum2d = reshape(nsst%z_c, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING Z_C RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "c_0", id_var)
   call netcdf_err(error, 'reading c_0 id' )
   dum2d = reshape(nsst%c_0, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING C_0 RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "c_d", id_var)
   call netcdf_err(error, 'reading c_d id' )
   dum2d = reshape(nsst%c_d, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING C_D RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "w_0", id_var)
   call netcdf_err(error, 'reading w_0 id' )
   dum2d = reshape(nsst%w_0, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING W_0 RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "w_d", id_var)
   call netcdf_err(error, 'reading w_d id' )
   dum2d = reshape(nsst%w_d, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING W_D RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xt", id_var)
   call netcdf_err(error, 'reading xt id' )
   dum2d = reshape(nsst%xt, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XT RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xs", id_var)
   call netcdf_err(error, 'reading xs id' )
   dum2d = reshape(nsst%xs, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XS RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xu", id_var)
   call netcdf_err(error, 'reading xu id' )
   dum2d = reshape(nsst%xu, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XU RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xv", id_var)
   call netcdf_err(error, 'reading xv id' )
   dum2d = reshape(nsst%xv, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XV RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xz", id_var)
   call netcdf_err(error, 'reading xz id' )
   dum2d = reshape(nsst%xz, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XZ RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "zm", id_var)
   call netcdf_err(error, 'reading zm id' )
   dum2d = reshape(nsst%zm, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING ZM RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xtts", id_var)
   call netcdf_err(error, 'reading xtts id' )
   dum2d = reshape(nsst%xtts, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XTTS RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "xzts", id_var)
   call netcdf_err(error, 'reading xzts id' )
   dum2d = reshape(nsst%xzts, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING XZTS RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "d_conv", id_var)
   call netcdf_err(error, 'reading d_conv id' )
   dum2d = reshape(nsst%d_conv, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING D_CONV RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "ifd", id_var)
   call netcdf_err(error, 'reading idf id' )
   dum2d = reshape(nsst%ifd, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING IFD RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "dt_cool", id_var)
   call netcdf_err(error, 'reading dt_cool id' )
   dum2d = reshape(nsst%dt_cool, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING DT_COOL RECORD' )
   call remove_checksum(ncid, id_var)

   error=nf90_inq_varid(ncid, "qrain", id_var)
   call netcdf_err(error, 'reading qrain id' )
   dum2d = reshape(nsst%qrain, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING QRAIN RECORD' )
   call remove_checksum(ncid, id_var)

! Some files don't include 'tfinc', which is just diagnostic. 
! If missing, then add it to the restart file.
   error=nf90_inq_varid(ncid, "tfinc", id_var)
   if (error /= 0) then
     error=nf90_inq_dimid(ncid, "xaxis_1", dim_x)
     call netcdf_err(error, 'finding xaxis_1' )
     error=nf90_inq_dimid(ncid, "yaxis_1", dim_y)
     call netcdf_err(error, 'finding yaxis_1' )
     error=nf90_inq_dimid(ncid, "Time", dim_time)
     call netcdf_err(error, 'finding Time' )
     dims_3d(1) = dim_x
     dims_3d(2) = dim_y
     dims_3d(3) = dim_time
     error=nf90_redef(ncid)
     error = nf90_def_var(ncid, 'tfinc', NF90_DOUBLE, dims_3d, id_var)
     call netcdf_err(error, 'DEFINING tfinc' )
     error = nf90_put_att(ncid, id_var, "long_name", "tfinc")
     call netcdf_err(error, 'DEFINING tfinc LONG NAME' )
     error = nf90_put_att(ncid, id_var, "units", "none")
     call netcdf_err(error, 'DEFINING tfinc UNITS' )
     error=nf90_enddef(ncid)
   endif
   dum2d = reshape(nsst%tfinc, (/idim,jdim/))
   error = nf90_put_var( ncid, id_var, dum2d)
   call netcdf_err(error, 'WRITING TFINC RECORD' )

 endif

 error = nf90_close(ncid)

 end subroutine write_data

!> Remove the checksum attribute from a netcdf record.
!!
!! @param[in] ncid netcdf file id
!! @param[in] id_var netcdf variable id.
!!
!! @author George Gayno NCEP/EMC
 subroutine remove_checksum(ncid, id_var)

 implicit none

 integer, intent(in)       :: ncid, id_var

 integer                   :: error

 error=nf90_inquire_attribute(ncid, id_var, 'checksum')

 if (error == 0) then ! attribute was found

   error = nf90_redef(ncid)
   call netcdf_err(error, 'entering define mode' )

   error=nf90_del_att(ncid, id_var, 'checksum')
   call netcdf_err(error, 'deleting checksum' )

   error= nf90_enddef(ncid)
   call netcdf_err(error, 'ending define mode' )

 endif

 end subroutine remove_checksum

 !> Read latitude and longitude for the cubed-sphere tile from the
 !! 'grid' file.  Read the filtered and unfiltered orography and
 !! optionally the land fraction from the 'orography' file.
 !!
 !! @param[in] IDIM 'i' dimension of cubed-sphere tile.
 !! @param[in] JDIM 'j' dimension of cubed-sphere tile.
 !! @param[in] IJDIM Total number of points on the cubed-sphere tile.
 !! @param[out] RLA Latitude on the cubed-sphere tile.
 !! @param[out] RLO Longitude on the cubed-sphere tile.
 !! @param[out] OROG Filtered orography.
 !! @param[out] OROG_UF Unfiltered orography.
 !! @param[out] TILE_NUM Cubed-sphere tile number
 !! @param[out] LANDFRAC Land fraction.
 !! @author George Gayno NOAA/EMC
 SUBROUTINE READ_LAT_LON_OROG(RLA,RLO,OROG,OROG_UF,&
                              TILE_NUM,IDIM,JDIM,IJDIM,LANDFRAC)

 USE MPI
 USE MACHINE

 IMPLICIT NONE

 INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM

 CHARACTER(LEN=5), INTENT(OUT) :: TILE_NUM

 REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
 REAL, INTENT(OUT)      :: OROG(IJDIM),OROG_UF(IJDIM)
 REAL(KIND=KIND_IO8), INTENT(OUT), OPTIONAL  :: LANDFRAC(IJDIM)

 CHARACTER(LEN=50)      :: FNOROG, FNGRID
 CHARACTER(LEN=3)       :: RANKCH

 INTEGER                :: ERROR, NCID, NCID_OROG
 INTEGER                :: I, II, J, JJ, MYRANK
 INTEGER                :: ID_DIM, ID_VAR, NX, NY

 REAL, ALLOCATABLE         :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
 REAL(KIND=4), ALLOCATABLE :: DUMMY4(:,:)

 CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)

 WRITE(RANKCH, '(I3.3)') (MYRANK+1)

 FNGRID = "./fngrid." // RANKCH

 PRINT*
 PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(FNGRID)

 ERROR=NF90_OPEN(TRIM(FNGRID),NF90_NOWRITE,NCID)
 CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNGRID) )

 ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
 CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )

 ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
 CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )

 IF ((NX/2) /= IDIM .OR. (NY/2) /= JDIM) THEN
   PRINT*,'FATAL ERROR: DIMENSIONS IN FILE: ',(NX/2),(NY/2)
   PRINT*,'DO NOT MATCH GRID DIMENSIONS: ',IDIM,JDIM
   CALL MPI_ABORT(MPI_COMM_WORLD, 130, ERROR)
 ENDIF

 ALLOCATE(GEOLON(NX+1,NY+1))
 ALLOCATE(GEOLAT(NX+1,NY+1))

 ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR)
 CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
 CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )

 ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR)
 CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
 CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )

 ALLOCATE(DUMMY(IDIM,JDIM))

 DO J = 1, JDIM
   DO I = 1, IDIM
     II = 2*I
     JJ = 2*J
     DUMMY(I,J) = GEOLON(II,JJ)
   ENDDO
 ENDDO

 RLO = RESHAPE(DUMMY, (/IJDIM/))

 DEALLOCATE(GEOLON)

 DO J = 1, JDIM
   DO I = 1, IDIM
     II = 2*I
     JJ = 2*J
     DUMMY(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 RLA = RESHAPE(DUMMY, (/IJDIM/))

 DEALLOCATE(GEOLAT, DUMMY)

 ERROR=NF90_INQ_VARID(NCID, 'tile', ID_VAR)
 CALL NETCDF_ERR(ERROR, 'ERROR READING TILE ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TILE_NUM)
 CALL NETCDF_ERR(ERROR, 'ERROR READING TILE RECORD' )

 ERROR = NF90_CLOSE(NCID)

 FNOROG = "./fnorog." // RANKCH

 PRINT*
 PRINT*, "READ FV3 OROG INFO FROM: "//TRIM(FNOROG)

 ERROR=NF90_OPEN(TRIM(FNOROG),NF90_NOWRITE,NCID_OROG)
 CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNOROG) )

 ALLOCATE(DUMMY4(IDIM,JDIM))

 ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_raw', ID_VAR)
 CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw ID' )
 ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
 CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw RECORD' )
 OROG_UF = RESHAPE(DUMMY4, (/IJDIM/))

 ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_filt', ID_VAR)
 CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt ID' )
 ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
 CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt RECORD' )
 OROG = RESHAPE(DUMMY4, (/IJDIM/))

 IF(PRESENT(LANDFRAC))THEN
   ERROR=NF90_INQ_VARID(NCID_OROG, 'land_frac', ID_VAR)
   CALL NETCDF_ERR(ERROR, 'ERROR READING land_frac ID' )
   ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
   CALL NETCDF_ERR(ERROR, 'ERROR READING land_frac RECORD' )
   LANDFRAC = RESHAPE(DUMMY4, (/IJDIM/))
 ENDIF

 DEALLOCATE(DUMMY4)

 ERROR = NF90_CLOSE(NCID_OROG)

 END SUBROUTINE READ_LAT_LON_OROG

 !> If a NetCDF call returns an error, print out a user-supplied
 !! message and the NetCDF library message.  Then stop processing.
 !!
 !! @param[in] ERR NetCDF error code.
 !! @param[in] STRING User-defined error message.
 !! @author George Gayno NOAA/EMC
 SUBROUTINE NETCDF_ERR( ERR, STRING )

 USE MPI

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: ERR
 CHARACTER(LEN=*), INTENT(IN) :: STRING
 CHARACTER(LEN=80) :: ERRMSG
 INTEGER :: IRET

 IF( ERR == NF90_NOERR )RETURN
 ERRMSG = NF90_STRERROR(ERR)
 PRINT*,''
 PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
 PRINT*,'STOP.'
 CALL MPI_ABORT(MPI_COMM_WORLD, 999, IRET)

 RETURN
 END SUBROUTINE NETCDF_ERR

 !> Read increment file from the GSI containing either the foundation 
 !! temperature increments and mask, or the soil increments. 
 !!
 !! The data is in NetCDF and on a gaussian grid. The grid contains two
 !! extra rows for each pole. The interpolation from gaussian to
 !! native grid assumes no pole points, so these are removed.
 !!
 !! @param[in] GSI_FILE Path/name of the GSI file to be read.
 !! @param[in] FILE_TYPE file-type to be read in, 'NST' or 'LND'.
 !! @param[in] LSOIL Number of model soil levels.
 !! 
 !! @author George Gayno NOAA/EMC
 SUBROUTINE READ_GSI_DATA(GSI_FILE, FILE_TYPE, LSOIL) 

 IMPLICIT NONE

 CHARACTER(LEN=*), INTENT(IN)     :: GSI_FILE
 CHARACTER(LEN=3), INTENT(IN)     :: FILE_TYPE 
 INTEGER, INTENT(IN), OPTIONAL    :: LSOIL

 INTEGER                          :: ERROR, ID_DIM, NCID
 INTEGER                          :: ID_VAR, J

 INTEGER(KIND=1), ALLOCATABLE     :: IDUMMY(:,:)

 REAL(KIND=8), ALLOCATABLE        :: DUMMY(:,:)

 CHARACTER(LEN=1)                 :: K_CH
 CHARACTER(LEN=10)                :: INCVAR
 CHARACTER(LEN=80)                :: err_msg
 INTEGER                          :: K

 PRINT*
 PRINT*, "READ INPUT GSI DATA FROM: "//TRIM(GSI_FILE)

 ERROR=NF90_OPEN(TRIM(GSI_FILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(GSI_FILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'latitude', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'READING latitude ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM_GAUS)
 CALL NETCDF_ERR(ERROR, 'READING latitude length' )
 JDIM_GAUS = JDIM_GAUS - 2  ! WILL IGNORE POLE POINTS
 
 ERROR=NF90_INQ_DIMID(NCID, 'longitude', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'READING longitude ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM_GAUS)
 CALL NETCDF_ERR(ERROR, 'READING longitude length' )

 IF (FILE_TYPE=='NST') then
     ALLOCATE(DUMMY(IDIM_GAUS,JDIM_GAUS+2))
     ALLOCATE(DTREF_GAUS(IDIM_GAUS,JDIM_GAUS))

     ERROR=NF90_INQ_VARID(NCID, "dtf", ID_VAR)
     CALL NETCDF_ERR(ERROR, 'READING dtf ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
     CALL NETCDF_ERR(ERROR, 'READING dtf' )

     ALLOCATE(IDUMMY(IDIM_GAUS,JDIM_GAUS+2))
     ALLOCATE(SLMSK_GAUS(IDIM_GAUS,JDIM_GAUS))

     ERROR=NF90_INQ_VARID(NCID, "msk", ID_VAR)
     CALL NETCDF_ERR(ERROR, 'READING msk ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, IDUMMY)
     CALL NETCDF_ERR(ERROR, 'READING msk' )

    ! REMOVE POLE POINTS.

     DO J = 1, JDIM_GAUS
       SLMSK_GAUS(:,J) = IDUMMY(:,J+1)
       DTREF_GAUS(:,J) = DUMMY(:,J+1)
     ENDDO

 ELSEIF (FILE_TYPE=='LND') then 

     ALLOCATE(DUMMY(IDIM_GAUS,JDIM_GAUS+2))
     ALLOCATE(STC_INC_GAUS(LSOIL,IDIM_GAUS,JDIM_GAUS))
     ALLOCATE(SLC_INC_GAUS(LSOIL,IDIM_GAUS,JDIM_GAUS))

     ! read in soil temperature increments in each layer
     DO K = 1, LSOIL
         WRITE(K_CH, '(I1)') K

         INCVAR = "soilt"//K_CH//"_inc"
         ERROR=NF90_INQ_VARID(NCID, INCVAR, ID_VAR)
         err_msg = "reading "//INCVAR//" ID" 
         CALL NETCDF_ERR(ERROR, trim(err_msg))
         ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
         err_msg = "reading "//INCVAR//" data"
         CALL NETCDF_ERR(ERROR, err_msg) 

         DO J = 1, JDIM_GAUS
           STC_INC_GAUS(K,:,J) = DUMMY(:,J+1)
         ENDDO

         INCVAR = "slc"//K_CH//"_inc"
         ERROR=NF90_INQ_VARID(NCID, INCVAR, ID_VAR)
         err_msg = "reading "//INCVAR//" ID"
         CALL NETCDF_ERR(ERROR, trim(err_msg))
         ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
         err_msg = "reading "//INCVAR//" data"
         CALL NETCDF_ERR(ERROR, err_msg)

         DO J = 1, JDIM_GAUS
           SLC_INC_GAUS(K,:,J) = DUMMY(:,J+1)
         ENDDO

     ENDDO

     ALLOCATE(IDUMMY(IDIM_GAUS,JDIM_GAUS+2))
     ALLOCATE(SOILSNOW_GAUS(IDIM_GAUS,JDIM_GAUS))

     ERROR=NF90_INQ_VARID(NCID, "soilsnow_mask", ID_VAR)
     CALL NETCDF_ERR(ERROR, 'READING soilsnow_mask ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, IDUMMY)
     CALL NETCDF_ERR(ERROR, 'READING soilsnow_mask' )

    ! REMOVE POLE POINTS.

     DO J = 1, JDIM_GAUS
       SOILSNOW_GAUS(:,J) = IDUMMY(:,J+1)
     ENDDO


 ELSE 
    print *, 'WARNING: FILE_TYPE', FILE_TYPE, 'not recognised.', &      
                ', no increments read in'  
 ENDIF

 IF(ALLOCATED(DUMMY)) DEALLOCATE(DUMMY)
 IF(ALLOCATED(IDUMMY)) DEALLOCATE(IDUMMY)

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_GSI_DATA

 !> Read the first guess surface records and nsst records (if
 !! selected) for a single cubed-sphere tile.
 !!
 !! @param[in] LSOIL Number of soil layers.
 !! @param[in] LENSFC Total number of points on a tile.
 !! @param[in] DO_NSST When true, nsst fields are read.
 !! @param[in] INC_FILE When true, read from an increment file.
 !!                     False reads from a restart file.
 !! @param[out] IS_NOAHMP When true, process for the Noah-MP LSM.
 !! @param[out] TSFFCS Skin Temperature.
 !! @param[out] SMCFCS Total volumetric soil moisture.
 !! @param[out] SWEFCS Snow water equivalent.
 !! @param[out] STCFCS Soil temperature.
 !! @param[out] TG3FCS Soil substrate temperature.
 !! @param[out] ZORFCS Roughness length.
 !! @param[out] CVFCS Cloud cover.
 !! @param[out] CVBFCS Cloud base.
 !! @param[out] CVTFCS Cloud top.
 !! @param[out] ALBFCS Snow-free albedo.
 !! @param[out] SLIFCS Land-sea mask including ice flag.
 !! @param[out] VEGFCS Vegetation greenness.
 !! @param[out] CNPFCS Plant canopy moisture content.
 !! @param[out] F10M log((z0+10)/z0). See model routine sfc_diff.f for details.
 !! @param[out] VETFCS Vegetation type.
 !! @param[out] SOTFCS Soil type.
 !! @param[out] ALFFCS Fractional coverage for strong/weak zenith angle
 !! dependent albedo.
 !! @param[out] USTAR Friction velocity.
 !! @param[out] FMM log((z0+z1)/z0). See model routine sfc_diff.f for details.
 !! @param[out] FHH log((ztmax+z1)/ztmax). See model routine sfc_diff.f for
 !! details.
 !! @param[out] SIHFCS Sea ice depth.
 !! @param[out] SICFCS Sea ice concentration.
 !! @param[out] SITFCS Sea ice temperature.
 !! @param[out] TPRCP Precipitation.
 !! @param[out] SRFLAG Snow/rain flag.
 !! @param[out] SNDFCS Snow depth.
 !! @param[out] VMNFCS Minimum vegetation greenness.
 !! @param[out] VMXFCS Maximum vegetation greenness.
 !! @param[out] SLCFCS Liquid portion of volumetric soil moisture.
 !! @param[out] SLPFCS Slope type.
 !! @param[out] ABSFCS Maximum snow albedo.
 !! @param[out] T2M Two-meter air temperature.
 !! @param[out] Q2M Two-meter specific humidity.
 !! @param[out] SLMASK Land-sea mask without ice flag.
 !! @param[out] ZSOIL Soil layer thickness.
 !! @param[out] NSST Data structure containing nsst fields.
 !! @author George Gayno NOAA/EMC
 SUBROUTINE READ_DATA(LSOIL,LENSFC,DO_NSST,INC_FILE,IS_NOAHMP, &
                      TSFFCS,SMCFCS,SWEFCS,STCFCS, &
                      TG3FCS,ZORFCS, &
                      CVFCS,CVBFCS,CVTFCS,ALBFCS, &
                      VEGFCS,SLIFCS,CNPFCS,F10M, &
                      VETFCS,SOTFCS,ALFFCS, &
                      USTAR,FMM,FHH, &
                      SIHFCS,SICFCS,SITFCS, &
                      TPRCP,SRFLAG,SNDFCS,  &
                      VMNFCS,VMXFCS,SLCFCS, &
                      SLPFCS,ABSFCS,T2M,Q2M,SLMASK, &
                      ZSOIL,NSST)
 USE MPI

 IMPLICIT NONE

 INTEGER, INTENT(IN)       :: LSOIL, LENSFC
 LOGICAL, INTENT(IN)       :: DO_NSST, INC_FILE

 LOGICAL, OPTIONAL, INTENT(OUT)      :: IS_NOAHMP

 REAL, OPTIONAL, INTENT(OUT)         :: CVFCS(LENSFC), CVBFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: CVTFCS(LENSFC), ALBFCS(LENSFC,4)
 REAL, OPTIONAL, INTENT(OUT)         :: SLIFCS(LENSFC), CNPFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: VEGFCS(LENSFC), F10M(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: VETFCS(LENSFC), SOTFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: TSFFCS(LENSFC), SWEFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: TG3FCS(LENSFC), ZORFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: ALFFCS(LENSFC,2), USTAR(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: FMM(LENSFC), FHH(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: SIHFCS(LENSFC), SICFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: SITFCS(LENSFC), TPRCP(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: SRFLAG(LENSFC), SNDFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: VMNFCS(LENSFC), VMXFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: SLPFCS(LENSFC), ABSFCS(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: T2M(LENSFC), Q2M(LENSFC), SLMASK(LENSFC)
 REAL, OPTIONAL, INTENT(OUT)         :: SLCFCS(LENSFC,LSOIL)
 REAL, OPTIONAL, INTENT(OUT)         :: SMCFCS(LENSFC,LSOIL)
 REAL, OPTIONAL, INTENT(OUT)         :: STCFCS(LENSFC,LSOIL)
 REAL(KIND=4), OPTIONAL, INTENT(OUT) :: ZSOIL(LSOIL)

 TYPE(NSST_DATA), OPTIONAL           :: NSST ! intent(out) will crash 
                                             ! because subtypes are allocated in main.

 CHARACTER(LEN=50)         :: FNBGSI
 CHARACTER(LEN=3)          :: RANKCH

 INTEGER                   :: ERROR, ERROR2, NCID, MYRANK
 INTEGER                   :: IDIM, JDIM, ID_DIM
 INTEGER                   :: ID_VAR, IERR

 REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:), DUMMY3D(:,:,:)

 CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)

 WRITE(RANKCH, '(I3.3)') (MYRANK+1)
 
 IF (INC_FILE) THEN
        FNBGSI = "./xainc." // RANKCH
 ELSE
        FNBGSI = "./fnbgsi." // RANKCH
 ENDIF

 PRINT*
 PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

 ERROR=NF90_OPEN(TRIM(FNBGSI),NF90_NOWRITE,NCID)
 CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNBGSI) )

 ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
 CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

 ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
 CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
 CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

 IF ((IDIM*JDIM) /= LENSFC) THEN
   PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
   CALL MPI_ABORT(MPI_COMM_WORLD, 88, IERR)
 ENDIF

! Check for records that indicate the restart file is
! for the Noah-MP land surface model.

 IF(PRESENT(IS_NOAHMP))THEN
   ERROR=NF90_INQ_VARID(NCID, "canliqxy", ID_VAR)
   ERROR2=NF90_INQ_VARID(NCID, "tsnoxy", ID_VAR)
   IS_NOAHMP=.FALSE.
   IF(ERROR == 0 .AND. ERROR2 == 0) THEN
     IS_NOAHMP=.TRUE.
     print*,"- WILL PROCESS FOR NOAH-MP LSM."
   ENDIF
 ENDIF

 ALLOCATE(DUMMY(IDIM,JDIM))

 IF (PRESENT(TSFFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "tsea", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING tsea ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING tsea' )
 TSFFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SWEFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING sheleg' )
 SWEFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(TG3FCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "tg3", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING tg3 ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING tg3' )
 TG3FCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(ZORFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "zorl", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING zorl ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING zorl' )
 ZORFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(ALBFCS)) THEN

 ERROR=NF90_INQ_VARID(NCID, "alvsf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING alvsf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING alvsf' )
 ALBFCS(:,1) = RESHAPE(DUMMY, (/LENSFC/))

 ERROR=NF90_INQ_VARID(NCID, "alvwf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING alvwf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING alvwf' )
 ALBFCS(:,2) = RESHAPE(DUMMY, (/LENSFC/))

 ERROR=NF90_INQ_VARID(NCID, "alnsf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING alnsf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING alnsf' )
 ALBFCS(:,3) = RESHAPE(DUMMY, (/LENSFC/))

 ERROR=NF90_INQ_VARID(NCID, "alnwf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING alnwf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING alnwf' )
 ALBFCS(:,4) = RESHAPE(DUMMY, (/LENSFC/))
  
 ENDIF

 IF (PRESENT(SLIFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING slmsk' )
 SLIFCS = RESHAPE(DUMMY, (/LENSFC/))
 SLMASK = SLIFCS
 WHERE (SLMASK > 1.5) SLMASK=0.0  ! remove sea ice
 ENDIF

 IF (PRESENT(CNPFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "canopy", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING canopy ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING canopy' )
 CNPFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(VEGFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "vfrac", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING vfrac ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING vfrac' )
 VEGFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(F10M)) THEN
 ERROR=NF90_INQ_VARID(NCID, "f10m", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING f10m ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING f10m' )
 F10M = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(VETFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING vtype' )
 VETFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SOTFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "stype", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING stype ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING stype' )
 SOTFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(ALFFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "facsf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING facsf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING facsf' )
 ALFFCS(:,1) = RESHAPE(DUMMY, (/LENSFC/))
  
 ERROR=NF90_INQ_VARID(NCID, "facwf", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING facwf ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING facwf' )
 ALFFCS(:,2) = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(USTAR)) THEN
 ERROR=NF90_INQ_VARID(NCID, "uustar", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING uustar ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING uustar' )
 USTAR = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(FMM)) THEN
 ERROR=NF90_INQ_VARID(NCID, "ffmm", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING ffmm ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING ffmm' )
 FMM = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(FHH)) THEN
 ERROR=NF90_INQ_VARID(NCID, "ffhh", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING ffhh ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING ffhh' )
 FHH = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SIHFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "hice", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING hice ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING hice' )
 SIHFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SICFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "fice", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING fice ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING fice' )
 SICFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SITFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "tisfc", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING tisfc ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING tisfc' )
 SITFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(TPRCP)) THEN
 ERROR=NF90_INQ_VARID(NCID, "tprcp", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING tprcp ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING tprcp' )
 TPRCP = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SRFLAG)) THEN
 ERROR=NF90_INQ_VARID(NCID, "srflag", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING srflag ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING srflag' )
 SRFLAG = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SNDFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING snwdph' )
 SNDFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(VMNFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "shdmin", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING shdmin ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING shdmin' )
 VMNFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(VMXFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "shdmax", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING shdmax ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING shdmax' )
 VMXFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(SLPFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "slope", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING slope ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING slope' )
 SLPFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(ABSFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "snoalb", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING snoalb ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING snoalb' )
 ABSFCS = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(T2M)) THEN
 ERROR=NF90_INQ_VARID(NCID, "t2m", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING t2m ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING t2m' )
 T2M = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF

 IF (PRESENT(Q2M)) THEN
 ERROR=NF90_INQ_VARID(NCID, "q2m", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING q2m ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
 CALL NETCDF_ERR(ERROR, 'READING q2m' )
 Q2M = RESHAPE(DUMMY, (/LENSFC/))
 ENDIF
  
 NSST_READ : IF(DO_NSST) THEN

   PRINT*
   PRINT*,"WILL READ NSST RECORDS."

   ERROR=NF90_INQ_VARID(NCID, "c_0", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING c_0 ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING c_0' )
   NSST%C_0 = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "c_d", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING c_d ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING c_d' )
   NSST%C_D = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "d_conv", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING d_conv ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING d_conv' )
   NSST%D_CONV = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "dt_cool", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING dt_cool ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING dt_cool' )
   NSST%DT_COOL = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "ifd", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING ifd ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING ifd' )
   NSST%IFD = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "qrain", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING qrain ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING qrain' )
   NSST%QRAIN = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "tref", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING tref ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING tref' )
   NSST%TREF = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "w_0", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING w_0 ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING w_0' )
   NSST%W_0 = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "w_d", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING w_d ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING w_d' )
   NSST%W_D = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xs", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xs ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xs' )
   NSST%XS = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xt", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xt ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xt' )
   NSST%XT = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xtts", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xtts ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xtts' )
   NSST%XTTS = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xu", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xu ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xu' )
   NSST%XU = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xv", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xv ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xv' )
   NSST%XV = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xz", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xz ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xz' )
   NSST%XZ = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "xzts", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING xzts ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING xzts' )
   NSST%XZTS = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "z_c", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING z_c ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING z_c' )
   NSST%Z_C = RESHAPE(DUMMY, (/LENSFC/))

   ERROR=NF90_INQ_VARID(NCID, "zm", ID_VAR)
   CALL NETCDF_ERR(ERROR, 'READING zm ID' )
   ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
   CALL NETCDF_ERR(ERROR, 'READING zm' )
   NSST%ZM = RESHAPE(DUMMY, (/LENSFC/))

 END IF NSST_READ

 DEALLOCATE(DUMMY)

 ALLOCATE(DUMMY3D(IDIM,JDIM,LSOIL))

 IF (PRESENT(SMCFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "smc", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING smc ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy3d)
 CALL NETCDF_ERR(ERROR, 'READING smc' )
 SMCFCS = RESHAPE(DUMMY3D, (/LENSFC,LSOIL/))
 ENDIF

 IF (PRESENT(SLCFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "slc", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING slc ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy3d)
 CALL NETCDF_ERR(ERROR, 'READING slc' )
 SLCFCS = RESHAPE(DUMMY3D, (/LENSFC,LSOIL/))
 ENDIF

 IF (PRESENT(STCFCS)) THEN
 ERROR=NF90_INQ_VARID(NCID, "stc", ID_VAR)
 CALL NETCDF_ERR(ERROR, 'READING stc ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy3d)
 CALL NETCDF_ERR(ERROR, 'READING stc' )
 STCFCS = RESHAPE(DUMMY3D, (/LENSFC,LSOIL/))
 ENDIF

 DEALLOCATE(DUMMY3D)

! cloud fields not in warm restart files.  set to zero?

 IF (PRESENT(CVFCS)) CVFCS = 0.0
 IF (PRESENT(CVTFCS)) CVTFCS = 0.0
 IF (PRESENT(CVBFCS)) CVBFCS = 0.0

! soil layer thicknesses not in warm restart files.  hardwire
! for now.
 
 IF (PRESENT(ZSOIL)) THEN
 ZSOIL(1) = -0.1
 ZSOIL(2) = -0.4
 ZSOIL(3) = -1.0
 ZSOIL(4) = -2.0
 ENDIF

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_DATA
 
 !> Read a GRIB1 sst climatological analysis file.
 !!
 !! Read the sst analysis and save it as an expanded and
 !! transposed array.
 !!
 !! @note The data is stored from north to south, but this 
 !! routine flips the poles. 
 !!   
 !! @param[in] file_sst File name of the sst file.
 !! @param[in] mlat_sst 'j' dimension of the sst data.
 !! @param[in] mlon_sst 'i' dimension of the sst data.
 !! @param[in] mon The month of the year.
 !! @param[out] sst The sst analysis data.
 !! @param[out] rlats_sst The latitudes of the sst data points.
 !! @param[out] rlons_sst The longitudes of the sst data points.
 !! @author Xu Li NOAA/EMC @date 2019-03-13
subroutine read_tf_clim_grb(file_sst,sst,rlats_sst,rlons_sst,mlat_sst,mlon_sst,mon)

  use mpi

  implicit none

! declare passed variables and arrays
  character(*)                       , intent(in   ) :: file_sst
  integer                            , intent(in   ) :: mon,mlat_sst,mlon_sst
  real, dimension(mlat_sst)          , intent(  out) :: rlats_sst
  real, dimension(mlon_sst)          , intent(  out) :: rlons_sst
  real, dimension(mlon_sst,mlat_sst) , intent(  out) :: sst

! declare local parameters
  integer,parameter:: lu_sst = 21   ! fortran unit number of grib sst file
  real, parameter :: deg2rad = 3.141593/180.0

! declare local variables and arrays
  logical(1), allocatable, dimension(:)    ::  lb

  integer :: nlat_sst !< Latitudinal dimension of the sst data.
  integer :: nlon_sst !< Longitudinal dimension of the sst data.
  integer :: iret,ni,nj
  integer :: mscan,kb1,ierr
  integer :: jincdir,i,iincdir,kb2,kb3,kf,kg,k,j,jf
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

  real :: xsst0 !< Latitude of the origin.
  real :: ysst0 !< Longitude of the origin.
  real :: dres  !< Latitude/longitude increment.
  real, allocatable, dimension(:) :: f
  
!************+******************************************************************************
!
! open sst analysis file (grib)
  write(*,*) ' sstclm : ',file_sst
  call baopenr(lu_sst,trim(file_sst),iret)
  if (iret /= 0 ) then
     write(6,*)'FATAL ERROR in read_tf_clm_grb: error opening sst file.'
     CALL MPI_ABORT(MPI_COMM_WORLD, 111, ierr)
  endif

! define sst variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=11        ! sst variable
  jpds(6)=1         ! surface
  jpds(9) = mon
  call getgbh(lu_sst,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  nlat_sst = kgds(3)   ! number points on longitude circle (360)
  nlon_sst = kgds(2)   ! number points on latitude circle (720)

! write(*,*) 'nlat_sst, nlon_sst, mon : ',nlat_sst, nlon_sst, mon
! write(*,*) 'mlat_sst, mlon_sst : ',mlat_sst, mlon_sst

! allocate arrays
  allocate(lb(nlat_sst*nlon_sst))
  allocate(f(nlat_sst*nlon_sst))
  jf=nlat_sst*nlon_sst

! read in the analysis
  call getgb(lu_sst,0,jf,j,jpds,jgds,kf,k,kpds,kgds,lb,f,iret)
  if (iret /= 0) then
     write(6,*)'FATAL ERROR in read_tf_clm_grb: error reading sst analysis data record.'
     deallocate(lb,f)
     CALL MPI_ABORT(MPI_COMM_WORLD, 111, ierr)
  endif

  if ( (nlat_sst /= mlat_sst) .or. (nlon_sst /= mlon_sst) ) then
     write(6,*)'FATAL ERROR in read_rtg_org:  inconsistent dimensions.  mlat_sst,mlon_sst=',&
          mlat_sst,mlon_sst,' -versus- nlat_sst,nlon_sst=',nlat_sst,nlon_sst
     deallocate(lb,f)
     CALL MPI_ABORT(MPI_COMM_WORLD, 111, ierr)
  endif

!
! get xlats and xlons
!
  dres = 180.0/real(nlat_sst)
  ysst0 = 0.5*dres-90.0
  xsst0 = 0.5*dres

! get lat_sst & lon_sst
  do j = 1, nlat_sst
     rlats_sst(j) = ysst0 + real(j-1)*dres
  enddo

  do i = 1, nlon_sst 
     rlons_sst(i) = (xsst0 + real(i-1)*dres)
  enddo

! load dimensions and grid specs.  check for unusual values
  ni=kgds(2)                ! 720 for 0.5 x 0.5
  nj=kgds(3)                ! 360 for 0.5 x 0.5 resolution

  mscan=kgds(11)
  kb1=ibits(mscan,7,1)   ! i scan direction
  kb2=ibits(mscan,6,1)   ! j scan direction
  kb3=ibits(mscan,5,1)   ! (i,j) or (j,i)

! get i and j scanning directions from kb1 and kb2.
! 0 yields +1, 1 yields -1. +1 is west to east, -1 is east to west.
  iincdir = 1-kb1*2

! 0 yields -1, 1 yields +1. +1 is south to north, -1 is north to south.
  jincdir = kb2*2 - 1

! write(*,*) 'read_tf_clim_grb,iincdir,jincdir : ',iincdir,jincdir
  do k=1,kf

!    kb3 from scan mode indicates if i points are consecutive
!    or if j points are consecutive
     if(kb3==0)then     !  (i,j)
        i=(ni+1)*kb1+(mod(k-1,ni)+1)*iincdir
        j=(nj+1)*(1-kb2)+jincdir*((k-1)/ni+1)
     else                !  (j,i)
        j=(nj+1)*(1-kb2)+(mod(k-1,nj)+1)*jincdir
        i=(ni+1)*kb1+iincdir*((k-1)/nj+1)
     endif
     sst (i,j) = f(k)
  end do
     
  deallocate(lb,f)

  call baclose(lu_sst,iret)
  if (iret /= 0 ) then
     write(6,*)'FATAL ERROR in read_tf_clm_grb: error closing sst file.'
     CALL MPI_ABORT(MPI_COMM_WORLD, 121, ierr)
  endif
  
end subroutine read_tf_clim_grb

!> Get the i/j dimensions of RTG SST climatology file.
!! The file is GRIB1.
!!
!! @param[in] file_sst File name of the sst file.
!! @param[in] mlat_sst The 'j' dimension of the data.
!! @param[in] mlon_sst The 'i' dimension of the data.
!! @author Xu Li NOAA/EMC @date 2019-03-13
subroutine get_tf_clm_dim(file_sst,mlat_sst,mlon_sst)
  use mpi

  implicit none

! declare passed variables and arrays
  character(*)                                   , intent(in ) :: file_sst
  integer                                        , intent(out) :: mlat_sst,mlon_sst

! declare local parameters
  integer,parameter:: lu_sst = 21   ! fortran unit number of grib sst file

  integer :: iret
  integer :: kf,kg,k,j,ierr
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

!************+******************************************************************************
!
! open sst analysis file (grib)
  call baopenr(lu_sst,trim(file_sst),iret)
  if (iret /= 0 ) then
     write(6,*)'FATAL ERROR in get_tf_clm_dim: error opening sst file.'
     CALL MPI_ABORT(MPI_COMM_WORLD, 111, ierr)
  endif

! define sst variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=11        ! sst variable
  jpds(6)=1         ! surface
  jpds(9) = 1
  call getgbh(lu_sst,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  mlat_sst = kgds(3)   ! number points on longitude circle (360)
  mlon_sst = kgds(2)   ! number points on latitude circle (720)

  write(*,*) 'mlat_sst, mlon_sst : ',mlat_sst, mlon_sst

  call baclose(lu_sst,iret)
  if (iret /= 0 ) then
     write(6,*)'FATAL ERROR in get_tf_clm_dim: error closing sst file.'
     CALL MPI_ABORT(MPI_COMM_WORLD, 121, ierr)
  endif
end subroutine get_tf_clm_dim

!> Read the woa05 salinity monthly climatology file.
!! The file is NetCDF.
!!
!! @param[in] filename The name of the climatology file.
!! @param[in] nlat The 'j' dimension of the data in the file.
!! @param[in] nlon The 'i' dimension of the data in the file.
!! @param[in] itime The monthly record to read.
!! @param[out] xlats The latitude of the data points.
!! @param[out] xlons The longitude of the data points.
!! @param[out] sal The salinity.
!! @author Xu Li NOAA/EMC
subroutine read_salclm_gfs_nc(filename,sal,xlats,xlons,nlat,nlon,itime)
  use netcdf
  implicit none
  
  ! This is the name of the data file we will read.
  character (len=*),             intent(in)  :: filename
  integer,                       intent(in)  :: nlat,nlon
  integer,                       intent(in)  :: itime
  real,    dimension(nlat),      intent(out) :: xlats
  real,    dimension(nlon),      intent(out) :: xlons
  real,    dimension(nlon,nlat), intent(out) :: sal
! Local variables
  integer :: ncid

  integer, parameter :: ndims = 3
  character (len = *), parameter :: lat_name = "latitude"
  character (len = *), parameter :: lon_name = "longitude"
  character (len = *), parameter :: t_name = "time"
  character (len = *), parameter :: sal_name="sal"
  integer :: time_varid,lon_varid, lat_varid, sal_varid

  ! The start and count arrays will tell the netCDF library where to read our data.
  integer, dimension(ndims) :: start, count

  character (len = *), parameter :: units = "units"
  character (len = *), parameter :: sal_units = "psu" 
                                  ! PSU (Practical SalinitUnit). 1 PSU = 1g/kg
  character (len = *), parameter :: time_units = "months"
  character (len = *), parameter :: lat_units = "degrees_north"
  character (len = *), parameter :: lon_units = "degrees_east"

! Open the file. 
  call nc_check( nf90_open(filename, nf90_nowrite, ncid) )

! Get the varids of time, latitude, longitude & depth coordinate variables.
  call nc_check( nf90_inq_varid(ncid, t_name,   time_varid) )
  call nc_check( nf90_inq_varid(ncid, lat_name, lat_varid) )
  call nc_check( nf90_inq_varid(ncid, lon_name, lon_varid) )

! Read the time, latitude and longitude data.
! call nc_check( nf90_get_var(ncid, time_varid, ntime) )
  call nc_check( nf90_get_var(ncid, lat_varid,  xlats) )
  call nc_check( nf90_get_var(ncid, lon_varid,  xlons) )

! Get the varids of the sal netCDF variables.
  call nc_check( nf90_inq_varid(ncid, sal_name,sal_varid) )

! Read 1 record of nlat*nlon values, starting at the beginning 
! of the record (the (1, 1, 1, rec) element in the netCDF file).
  start = (/ 1, 1, itime /)
  count = (/ nlon, nlat, 1 /)

!  write(*,*) 'read_salclm_gfs_nc itime : ',itime
! Read the sal data from the file, one record at a time.
  call nc_check( nf90_get_var(ncid, sal_varid, sal, start, count) )

! Close the file. This frees up any internal netCDF resources
! associated with the file.
  call nc_check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee! 
!  print *,"*** SUCCESS reading file ", filename, "!"

end subroutine read_salclm_gfs_nc

!> Get the i/j dimensions of the data from a NetCDF file.
!!
!! @param[in] filename Name of the file to be read.
!! @param[out] nlat 'j' dimension of the data in the file.
!! @param[out] nlon 'i' dimension of the data in the file.
!! @author Xu Li NOAA/EMC
subroutine get_dim_nc(filename,nlat,nlon)
  use netcdf
  implicit none
  
  character (len=*), intent(in)  :: filename
  integer,           intent(out) :: nlat,nlon
! Local variables
  character (len = *), parameter :: lat_name = "latitude"
  character (len = *), parameter :: lon_name = "longitude"
  integer :: ncid
  integer :: LatDimID,LonDimID

! Open the file. 
  call nc_check( nf90_open(filename, nf90_nowrite, ncid) )

! Get dimensions 
  call nc_check( nf90_inq_dimid(ncid,lat_name,LatDimID) )
  call nc_check( nf90_inq_dimid(ncid,lon_name,LonDimID) )
  call nc_check( nf90_inquire_dimension(ncid,LatDimID,len=nlat) )
  call nc_check( nf90_inquire_dimension(ncid,LonDimID,len=nlon) )

!  write(*,'(a,1x,a6,2I8)') 'get_dim_nc, file, nlat, nlon : ',filename,nlat,nlon

! Close the file. This frees up any internal netCDF resources
! associated with the file.
  call nc_check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee! 
!  print *,"*** SUCCESS get dimensions from nc file ", filename, "!"

end subroutine get_dim_nc

!> Check the NetCDF status code. If there is an error,
!! print the library error message and stop processing.
!!
!! @param[in] status NetCDF status code.
!! @author Xu Li NOAA/EMC
subroutine nc_check(status)

  use mpi
  use netcdf

  integer, intent ( in) :: status
  integer :: ierr

  if(status /= nf90_noerr) then
    print *, 'FATAL ERROR:'
    print *, trim(nf90_strerror(status))
    CALL MPI_ABORT(MPI_COMM_WORLD, 122, ierr)
  end if
end subroutine nc_check

 END MODULE READ_WRITE_DATA
