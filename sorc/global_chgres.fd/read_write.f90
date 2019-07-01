 SUBROUTINE READ_FV3GFS_ATMS_DATA_NEMSIO(GFILEI, GFSDATAI, GFSHEADI, &
                                         VCOORD, LEVSP1, NVCOORD)

 USE NEMSIO_MODULE
 USE NEMSIO_GFS

 IMPLICIT NONE

 TYPE(NEMSIO_GFILE)  :: GFILEI
 TYPE(NEMSIO_DBTA)   :: GFSDATAI
 TYPE(NEMSIO_HEAD)   :: GFSHEADI

 INTEGER, INTENT(IN) :: LEVSP1, NVCOORD
 
 REAL, INTENT(IN)    :: VCOORD(LEVSP1, NVCOORD)

 INTEGER             :: I, J, L, IRET
 INTEGER             :: LONB, LATB, LEVSI

 REAL, ALLOCATABLE   :: TMP(:), P_INTERFACE(:)

 print*,''
 print*,'READ FV3GFS ATMOSPHERIC NEMSIO FILE'

 LONB = GFSHEADI%DIMX
 LATB = GFSHEADI%DIMY
 LEVSI = GFSHEADI%DIMZ

 ALLOCATE(TMP(LONB*LATB))

 PRINT*,'READ HGT'
 CALL NEMSIO_READRECV(GFILEI, 'hgt', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 GFSDATAI%ZS = RESHAPE(TMP,(/LONB,LATB/))

 PRINT*,'READ U WINDS'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'ugrd', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%U(:,:,L) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ V WINDS'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'vgrd', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%V(:,:,L) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ T'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'tmp', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%T(:,:,L) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ Q'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'spfh', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,1) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ O3'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'o3mr', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,2) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ CLWMR'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'clwmr', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,3) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ RWMR'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'rwmr', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,4) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ ICMR'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'icmr', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,5) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ SNMR'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'snmr', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,6) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ GRLE'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'grle', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%Q(:,:,L,7) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ DZDT'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'dzdt', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%W(:,:,L) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

 PRINT*,'READ DPRES'
 DO L = 1, LEVSI
   CALL NEMSIO_READRECV(GFILEI, 'dpres', 'mid layer', L, TMP, IRET=IRET)
   IF (IRET /= 0) GOTO 99
   GFSDATAI%DP(:,:,L) = RESHAPE(TMP,(/LONB,LATB/))
 ENDDO

! COMPUTE SURFACE PRESSURE AND MID-LAYER PRESSURE FROM DELTA P.
! DO NOT USE THE SURFACE PRESSURE IN THE FILE.  AFTER INTERPOLATION
! FROM THE MODEL GRID TO THE GAUSSIAN GRID IN THE WRITE COMPONENT,
! THE SURFACE PRESSURE IS NO LONGER CONSISTENT WITH DELTA P.

 ALLOCATE(P_INTERFACE(LEVSI+1))

 DO J = 1, LATB
 DO I = 1, LONB
   P_INTERFACE(LEVSI+1) = VCOORD(LEVSI+1,1) ! MODEL TOP PRESSURE
   DO L = LEVSI, 1, -1
     P_INTERFACE(L) = P_INTERFACE(L+1) + GFSDATAI%DP(I,J,L)
   ENDDO
   GFSDATAI%PS(I,J) = P_INTERFACE(1)  ! SURFACE PRESSURE
   DO L = 1, LEVSI
     GFSDATAI%P(I,J,L) = (P_INTERFACE(L) + P_INTERFACE(L+1)) * 0.5
   ENDDO
 ENDDO
 ENDDO

 DEALLOCATE(P_INTERFACE)

 DEALLOCATE(TMP)

 RETURN

 99 CONTINUE
 PRINT*,'FATAL ERROR READING FV3GFS ATMOSPHERIC NEMSIO FILE.'
 PRINT*,'IRET IS: ', IRET
 CALL ERREXIT(22)

 END SUBROUTINE READ_FV3GFS_ATMS_DATA_NEMSIO

 SUBROUTINE WRITE_FV3_ATMS_HEADER_NETCDF(LEVS_P1, NTRACM, NVCOORD, VCOORD)

 use netcdf

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: LEVS_P1
 INTEGER, INTENT(IN) :: NTRACM
 INTEGER, INTENT(IN) :: NVCOORD

 REAL, INTENT(IN)    :: VCOORD(LEVS_P1, NVCOORD)

 CHARACTER(LEN=13)   :: OUTFILE

 INTEGER             :: ERROR, NCID
 INTEGER             :: DIM_NVCOORD, DIM_LEVSP
 INTEGER             :: ID_NTRAC, ID_VCOORD
 INTEGER             :: FSIZE=65536, INITAL = 0
 INTEGER             :: HEADER_BUFFER_VAL = 16384

 REAL(KIND=8)        :: TMP(LEVS_P1,NVCOORD)

 OUTFILE = "./gfs_ctrl.nc"

 ERROR = NF90_CREATE(OUTFILE, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                     NCID, INITIALSIZE=INITAL, CHUNKSIZE=FSIZE)
 CALL NETCDF_ERROR(ERROR, 'Creating file '//TRIM(OUTFILE) )

 ERROR = NF90_DEF_DIM(NCID, 'nvcoord', NVCOORD, DIM_NVCOORD)
 CALL NETCDF_ERROR(ERROR, 'define dimension nvcoord for file='//TRIM(OUTFILE) )

 ERROR = NF90_DEF_DIM(NCID, 'levsp', LEVS_P1, DIM_LEVSP)
 CALL NETCDF_ERROR(ERROR, 'define dimension levsp for file='//TRIM(OUTFILE) )

 ERROR = NF90_DEF_VAR(NCID, 'ntrac', NF90_INT, ID_NTRAC)
 CALL NETCDF_ERROR(ERROR, 'define var ntrac for file='//TRIM(OUTFILE) )

 ERROR = NF90_DEF_VAR(NCID, 'vcoord', NF90_DOUBLE, (/DIM_LEVSP, DIM_NVCOORD/), ID_VCOORD)
 CALL NETCDF_ERROR(ERROR, 'define var vcoord for file='//TRIM(OUTFILE) )   

 ERROR = NF90_ENDDEF(NCID, HEADER_BUFFER_VAL,4,0,4)
 CALL NETCDF_ERROR(ERROR, 'end meta define for file='//TRIM(OUTFILE) )

 ERROR = NF90_PUT_VAR( NCID, ID_NTRAC, NTRACM)
 CALL NETCDF_ERROR(ERROR, 'write var ntrac for file='//TRIM(OUTFILE) )

 TMP(1:LEVS_P1,:) = VCOORD(LEVS_P1:1:-1,:)
 ERROR = NF90_PUT_VAR( NCID, ID_VCOORD, TMP)
 CALL NETCDF_ERROR(ERROR, 'write var vcoord for file='//TRIM(OUTFILE) )

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE WRITE_FV3_ATMS_HEADER_NETCDF

 subroutine netcdf_error( err, string )
 use netcdf
 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 call errexit(999)

 return
 end subroutine netcdf_error

 subroutine write_fv3_sfc_data_netcdf(lonb, latb, lsoil, sfcoutput, f10m, &
                           t2m, q2m, uustar, ffmm, ffhh, tprcp, srflag, tile, &
                           num_nsst_fields, nsst_output)

 use netcdf

 use surface_chgres, only        : sfc1d

 implicit none

 integer, intent(in)            :: latb, lonb, lsoil, tile
 integer, intent(in)            :: num_nsst_fields
 character(len=128)             :: outfile

 integer                        :: fsize=65536, inital = 0
 integer                        :: header_buffer_val = 16384
 integer                        :: dim_lon, dim_lat, dim_lsoil
 integer                        :: error, ncid, i
 integer                        :: id_lon, id_lat, id_lsoil
 integer                        :: id_geolon, id_geolat, id_slmsk
 integer                        :: id_tsea, id_sheleg, id_tg3
 integer                        :: id_zorl, id_alvsf, id_alvwf
 integer                        :: id_alnsf, id_alnwf, id_vfrac
 integer                        :: id_canopy, id_f10m, id_t2m
 integer                        :: id_q2m, id_vtype, id_stype
 integer                        :: id_facsf, id_facwf, id_uustar
 integer                        :: id_ffmm, id_ffhh, id_hice
 integer                        :: id_fice, id_tisfc, id_tprcp
 integer                        :: id_srflag, id_snwdph, id_shdmin
 integer                        :: id_shdmax, id_slope, id_snoalb
 integer                        :: id_stc, id_smc, id_slc
 integer                        :: id_tref, id_z_c, id_c_0
 integer                        :: id_c_d, id_w_0, id_w_d
 integer                        :: id_xt, id_xs, id_xu, id_xv
 integer                        :: id_xz, id_zm, id_xtts, id_xzts
 integer                        :: id_d_conv, id_ifd, id_dt_cool
 integer                        :: id_qrain
 
 logical                        :: write_nsst

 real, intent(in)               :: f10m(lonb,latb)
 real, intent(in)               :: q2m(lonb,latb)
 real, intent(in)               :: t2m(lonb,latb)
 real, intent(in)               :: uustar(lonb,latb)
 real, intent(in)               :: ffmm(lonb,latb)
 real, intent(in)               :: ffhh(lonb,latb)
 real, intent(in)               :: tprcp(lonb,latb)
 real, intent(in)               :: srflag(lonb,latb)
 real, intent(in), optional     :: nsst_output(lonb*latb,num_nsst_fields)
 real(kind=4)                   :: lsoil_data(lsoil)
 real(kind=4), allocatable      :: dum2d(:,:), dum3d(:,:,:)

 type(sfc1d)                    :: sfcoutput

 write_nsst = .false.
 if (present(nsst_output)) write_nsst = .true.

 if (write_nsst) then
   print*,'- WRITE FV3 SURFACE AND NSST DATA TO NETCDF FILE'
 else
   print*,'- WRITE FV3 SURFACE DATA TO NETCDF FILE'
 endif

 if (tile < 10) then
   write(outfile, '(A, I1, A)'), 'out.sfc.tile', tile, '.nc'
 else
   write(outfile, '(A, I2, A)'), 'out.sfc.tile', tile, '.nc'
 endif

!--- open the file
 error = nf90_create(outfile, ior(nf90_netcdf4,nf90_classic_model), &
                     ncid, initialsize=inital, chunksize=fsize)
 call netcdf_error(error, 'CREATING FILE='//trim(outfile) )

!--- define dimension
 error = nf90_def_dim(ncid, 'lon', lonb, dim_lon)
 call netcdf_error(error, 'DEFINING LON DIMENSION' )
 error = nf90_def_dim(ncid, 'lat', latb, dim_lat)
 call netcdf_error(error, 'DEFINING LAT DIMENSION' )
 error = nf90_def_dim(ncid, 'lsoil', lsoil, dim_lsoil)
 call netcdf_error(error, 'DEFINING LSOIL DIMENSION' )

 !--- define field
 error = nf90_def_var(ncid, 'lon', NF90_FLOAT, (/dim_lon/), id_lon)
 call netcdf_error(error, 'DEFINING LON FIELD' )
 error = nf90_put_att(ncid, id_lon, "cartesian_axis", "X")
 call netcdf_error(error, 'WRITING LON FIELD' )
 error = nf90_def_var(ncid, 'lat', NF90_FLOAT, (/dim_lat/), id_lat)
 call netcdf_error(error, 'DEFINING LAT FIELD' )
 error = nf90_put_att(ncid, id_lat, "cartesian_axis", "Y")
 call netcdf_error(error, 'WRITING LAT FIELD' )
 error = nf90_def_var(ncid, 'lsoil', NF90_FLOAT, (/dim_lsoil/), id_lsoil)
 call netcdf_error(error, 'DEFINING LSOIL FIELD' )
 error = nf90_put_att(ncid, id_lsoil, "cartesian_axis", "Z")
 call netcdf_error(error, 'WRITING LSOIL FIELD' )
 error = nf90_def_var(ncid, 'geolon', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolon)
 call netcdf_error(error, 'DEFINING GEOLON' )
 error = nf90_def_var(ncid, 'geolat', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolat)
 call netcdf_error(error, 'DEFINING GEOLAT' )
 error = nf90_def_var(ncid, 'slmsk', NF90_FLOAT, (/dim_lon,dim_lat/), id_slmsk)
 call netcdf_error(error, 'DEFINING SLMSK' )
 error = nf90_def_var(ncid, 'tsea', NF90_FLOAT, (/dim_lon,dim_lat/), id_tsea)
 call netcdf_error(error, 'DEFINING TSEA' )
 error = nf90_def_var(ncid, 'sheleg', NF90_FLOAT, (/dim_lon,dim_lat/), id_sheleg)
 call netcdf_error(error, 'DEFINING SHELEG' )
 error = nf90_def_var(ncid, 'tg3', NF90_FLOAT, (/dim_lon,dim_lat/), id_tg3)
 call netcdf_error(error, 'DEFINING TG3' )
 error = nf90_def_var(ncid, 'zorl', NF90_FLOAT, (/dim_lon,dim_lat/), id_zorl)
 call netcdf_error(error, 'DEFINING ZORL' )
 error = nf90_def_var(ncid, 'alvsf', NF90_FLOAT, (/dim_lon,dim_lat/), id_alvsf)
 call netcdf_error(error, 'DEFINING ALVSF' )
 error = nf90_def_var(ncid, 'alvwf', NF90_FLOAT, (/dim_lon,dim_lat/), id_alvwf)
 call netcdf_error(error, 'DEFINING ALVWF' )
 error = nf90_def_var(ncid, 'alnsf', NF90_FLOAT, (/dim_lon,dim_lat/), id_alnsf)
 call netcdf_error(error, 'DEFINING ALNSF' )
 error = nf90_def_var(ncid, 'alnwf', NF90_FLOAT, (/dim_lon,dim_lat/), id_alnwf)
 call netcdf_error(error, 'DEFINING ALNWF' )
 error = nf90_def_var(ncid, 'vfrac', NF90_FLOAT, (/dim_lon,dim_lat/), id_vfrac)
 call netcdf_error(error, 'DEFINING VFRAC' )
 error = nf90_def_var(ncid, 'canopy', NF90_FLOAT, (/dim_lon,dim_lat/), id_canopy)
 call netcdf_error(error, 'DEFINING CANOPY' )
 error = nf90_def_var(ncid, 'f10m', NF90_FLOAT, (/dim_lon,dim_lat/), id_f10m)
 call netcdf_error(error, 'DEFINING F10M' )
 error = nf90_def_var(ncid, 't2m', NF90_FLOAT, (/dim_lon,dim_lat/), id_t2m)
 call netcdf_error(error, 'DEFINING T2M' )
 error = nf90_def_var(ncid, 'q2m', NF90_FLOAT, (/dim_lon,dim_lat/), id_q2m)
 call netcdf_error(error, 'DEFINING Q2M' )
 error = nf90_def_var(ncid, 'vtype', NF90_FLOAT, (/dim_lon,dim_lat/), id_vtype)
 call netcdf_error(error, 'DEFINING VTYPE' )
 error = nf90_def_var(ncid, 'stype', NF90_FLOAT, (/dim_lon,dim_lat/), id_stype)
 call netcdf_error(error, 'DEFINING STYPE' )
 error = nf90_def_var(ncid, 'facsf', NF90_FLOAT, (/dim_lon,dim_lat/), id_facsf)
 call netcdf_error(error, 'DEFINING FACSF' )
 error = nf90_def_var(ncid, 'facwf', NF90_FLOAT, (/dim_lon,dim_lat/), id_facwf)
 call netcdf_error(error, 'DEFINING FACWF' )
 error = nf90_def_var(ncid, 'uustar', NF90_FLOAT, (/dim_lon,dim_lat/), id_uustar)
 call netcdf_error(error, 'DEFINING UUSTAR' )
 error = nf90_def_var(ncid, 'ffmm', NF90_FLOAT, (/dim_lon,dim_lat/), id_ffmm)
 call netcdf_error(error, 'DEFINING FFMM' )
 error = nf90_def_var(ncid, 'ffhh', NF90_FLOAT, (/dim_lon,dim_lat/), id_ffhh)
 call netcdf_error(error, 'DEFINING FFHH' )
 error = nf90_def_var(ncid, 'hice', NF90_FLOAT, (/dim_lon,dim_lat/), id_hice)
 call netcdf_error(error, 'DEFINING HICE' )
 error = nf90_def_var(ncid, 'fice', NF90_FLOAT, (/dim_lon,dim_lat/), id_fice)
 call netcdf_error(error, 'DEFINING FICE' )
 error = nf90_def_var(ncid, 'tisfc', NF90_FLOAT, (/dim_lon,dim_lat/), id_tisfc)
 call netcdf_error(error, 'DEFINING TISFC' )
 error = nf90_def_var(ncid, 'tprcp', NF90_FLOAT, (/dim_lon,dim_lat/), id_tprcp)
 call netcdf_error(error, 'DEFINING TPRCP' )
 error = nf90_def_var(ncid, 'srflag', NF90_FLOAT, (/dim_lon,dim_lat/), id_srflag)
 call netcdf_error(error, 'DEFINING SRFLAG' )
 error = nf90_def_var(ncid, 'snwdph', NF90_FLOAT, (/dim_lon,dim_lat/), id_snwdph)
 call netcdf_error(error, 'DEFINING SNWDPH' )
 error = nf90_def_var(ncid, 'shdmin', NF90_FLOAT, (/dim_lon,dim_lat/), id_shdmin)
 call netcdf_error(error, 'DEFINING SHDMIN' )
 error = nf90_def_var(ncid, 'shdmax', NF90_FLOAT, (/dim_lon,dim_lat/), id_shdmax)
 call netcdf_error(error, 'DEFINING SHDMAX' )
 error = nf90_def_var(ncid, 'slope', NF90_FLOAT, (/dim_lon,dim_lat/), id_slope)
 call netcdf_error(error, 'DEFINING SLOPE' )
 error = nf90_def_var(ncid, 'snoalb', NF90_FLOAT, (/dim_lon,dim_lat/), id_snoalb)
 call netcdf_error(error, 'DEFINING SNOALB' )
 error = nf90_def_var(ncid, 'stc', NF90_FLOAT, (/dim_lon,dim_lat,dim_lsoil/), id_stc)
 call netcdf_error(error, 'DEFINING STC' )
 error = nf90_def_var(ncid, 'smc', NF90_FLOAT, (/dim_lon,dim_lat,dim_lsoil/), id_smc)
 call netcdf_error(error, 'DEFINING SMC' )
 error = nf90_def_var(ncid, 'slc', NF90_FLOAT, (/dim_lon,dim_lat,dim_lsoil/), id_slc)
 call netcdf_error(error, 'DEFINING SLC' )
 if (write_nsst) then
   error = nf90_def_var(ncid, 'tref', NF90_FLOAT, (/dim_lon,dim_lat/), id_tref)
   call netcdf_error(error, 'DEFINING TREF' )
   error = nf90_def_var(ncid, 'z_c', NF90_FLOAT, (/dim_lon,dim_lat/), id_z_c)
   call netcdf_error(error, 'DEFINING Z_C' )
   error = nf90_def_var(ncid, 'c_0', NF90_FLOAT, (/dim_lon,dim_lat/), id_c_0)
   call netcdf_error(error, 'DEFINING C_0' )
   error = nf90_def_var(ncid, 'c_d', NF90_FLOAT, (/dim_lon,dim_lat/), id_c_d)
   call netcdf_error(error, 'DEFINING C_D' )
   error = nf90_def_var(ncid, 'w_0', NF90_FLOAT, (/dim_lon,dim_lat/), id_w_0)
   call netcdf_error(error, 'DEFINING W_0' )
   error = nf90_def_var(ncid, 'w_d', NF90_FLOAT, (/dim_lon,dim_lat/), id_w_d)
   call netcdf_error(error, 'DEFINING W_D' )
   error = nf90_def_var(ncid, 'xt', NF90_FLOAT, (/dim_lon,dim_lat/), id_xt)
   call netcdf_error(error, 'DEFINING XT' )
   error = nf90_def_var(ncid, 'xs', NF90_FLOAT, (/dim_lon,dim_lat/), id_xs)
   call netcdf_error(error, 'DEFINING XS' )
   error = nf90_def_var(ncid, 'xu', NF90_FLOAT, (/dim_lon,dim_lat/), id_xu)
   call netcdf_error(error, 'DEFINING XU' )
   error = nf90_def_var(ncid, 'xv', NF90_FLOAT, (/dim_lon,dim_lat/), id_xv)
   call netcdf_error(error, 'DEFINING XV' )
   error = nf90_def_var(ncid, 'xz', NF90_FLOAT, (/dim_lon,dim_lat/), id_xz)
   call netcdf_error(error, 'DEFINING XZ' )
   error = nf90_def_var(ncid, 'zm', NF90_FLOAT, (/dim_lon,dim_lat/), id_zm)
   call netcdf_error(error, 'DEFINING ZM' )
   error = nf90_def_var(ncid, 'xtts', NF90_FLOAT, (/dim_lon,dim_lat/), id_xtts)
   call netcdf_error(error, 'DEFINING XTTS' )
   error = nf90_def_var(ncid, 'xzts', NF90_FLOAT, (/dim_lon,dim_lat/), id_xzts)
   call netcdf_error(error, 'DEFINING XZTS' )
   error = nf90_def_var(ncid, 'd_conv', NF90_FLOAT, (/dim_lon,dim_lat/), id_d_conv)
   call netcdf_error(error, 'DEFINING D_CONV' )
   error = nf90_def_var(ncid, 'ifd', NF90_FLOAT, (/dim_lon,dim_lat/), id_ifd)
   call netcdf_error(error, 'DEFINING IFD' )
   error = nf90_def_var(ncid, 'dt_cool', NF90_FLOAT, (/dim_lon,dim_lat/), id_dt_cool)
   call netcdf_error(error, 'DEFINING DT_COOL' )
   error = nf90_def_var(ncid, 'qrain', NF90_FLOAT, (/dim_lon,dim_lat/), id_qrain)
   call netcdf_error(error, 'DEFINING QRAIN' )
 endif

 error = nf90_enddef(ncid, header_buffer_val, 4, 0, 4)
 call netcdf_error(error, 'DEFINING HEADER' )

 allocate(dum2d(lonb,latb))

 dum2d = reshape(sfcoutput%lons, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_lon, dum2d(:,1))
 call netcdf_error(error, 'WRITING LON HEADER RECORD' )

 dum2d = reshape(sfcoutput%lats, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_lat, dum2d(1,:))
 call netcdf_error(error, 'WRITING LAT HEADER RECORD' )

 do i = 1, lsoil
   lsoil_data(i) = float(i)
 enddo
 error = nf90_put_var( ncid, id_lsoil, lsoil_data)
 call netcdf_error(error, 'WRITING LSOIL HEADER' )

 dum2d = reshape(sfcoutput%lons, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_geolon, dum2d)
 call netcdf_error(error, 'WRITING GEOLON RECORD' )

 dum2d = reshape(sfcoutput%lats, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_geolat, dum2d)
 call netcdf_error(error, 'WRITING GEOLAT RECORD' )

 dum2d = reshape(sfcoutput%lsmask, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_slmsk, dum2d)
 call netcdf_error(error, 'WRITING SLMSK RECORD' )

 dum2d = reshape(sfcoutput%skin_temp, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_tsea, dum2d)
 call netcdf_error(error, 'WRITING TSEA RECORD' )

 dum2d = reshape(sfcoutput%snow_liq_equiv, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_sheleg, dum2d)
 call netcdf_error(error, 'WRITING SHELEG RECORD' )

 dum2d = reshape(sfcoutput%substrate_temp, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_tg3, dum2d)
 call netcdf_error(error, 'WRITING TG3 RECORD' )

 dum2d = reshape(sfcoutput%z0, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_zorl, dum2d)
 call netcdf_error(error, 'WRITING ZORL RECORD' )

 dum2d = reshape(sfcoutput%alvsf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_alvsf, dum2d)
 call netcdf_error(error, 'WRITING ALVSF RECORD' )

 dum2d = reshape(sfcoutput%alvwf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_alvwf, dum2d)
 call netcdf_error(error, 'WRITING ALVWF RECORD' )

 dum2d = reshape(sfcoutput%alnsf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_alnsf, dum2d)
 call netcdf_error(error, 'WRITING ALNSF RECORD' )

 dum2d = reshape(sfcoutput%alnwf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_alnwf, dum2d)
 call netcdf_error(error, 'WRITING ALNWF RECORD' )

 dum2d = reshape(sfcoutput%greenfrc, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_vfrac, dum2d)
 call netcdf_error(error, 'WRITING VFRAC RECORD' )

 dum2d = reshape(sfcoutput%canopy_mc, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_canopy, dum2d)
 call netcdf_error(error, 'WRITING CANOPY RECORD' )

 dum2d = f10m
 error = nf90_put_var( ncid, id_f10m, dum2d)
 call netcdf_error(error, 'WRITING F10M RECORD' )

 dum2d = t2m
 error = nf90_put_var( ncid, id_t2m, dum2d)
 call netcdf_error(error, 'WRITING T2M RECORD' )

 dum2d = q2m
 error = nf90_put_var( ncid, id_q2m, dum2d)
 call netcdf_error(error, 'WRITING Q2M RECORD' )

 dum2d = reshape(float(sfcoutput%veg_type), (/lonb,latb/) )
 error = nf90_put_var( ncid, id_vtype, dum2d)
 call netcdf_error(error, 'WRITING VTYPE RECORD' )

 dum2d = reshape(float(sfcoutput%soil_type), (/lonb,latb/) )
 error = nf90_put_var( ncid, id_stype, dum2d)
 call netcdf_error(error, 'WRITING STYPE RECORD' )

 dum2d = reshape(sfcoutput%facsf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_facsf, dum2d)
 call netcdf_error(error, 'WRITING FACSF RECORD' )

 dum2d = reshape(sfcoutput%facwf, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_facwf, dum2d)
 call netcdf_error(error, 'WRITING FACWF RECORD' )

 dum2d = uustar
 error = nf90_put_var( ncid, id_uustar, dum2d)
 call netcdf_error(error, 'WRITING UUSTAR RECORD' )

 dum2d = ffmm
 error = nf90_put_var( ncid, id_ffmm, dum2d)
 call netcdf_error(error, 'WRITING FFMM RECORD' )

 dum2d = ffhh
 error = nf90_put_var( ncid, id_ffhh, dum2d)
 call netcdf_error(error, 'WRITING FFHH RECORD' )

 dum2d = reshape(sfcoutput%sea_ice_depth, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_hice, dum2d)
 call netcdf_error(error, 'WRITING HICE RECORD' )

 dum2d = reshape(sfcoutput%sea_ice_fract, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_fice, dum2d)
 call netcdf_error(error, 'WRITING FICE RECORD' )

 dum2d = reshape(sfcoutput%sea_ice_temp, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_tisfc, dum2d)
 call netcdf_error(error, 'WRITING TISFC RECORD' )

 dum2d = tprcp
 error = nf90_put_var( ncid, id_tprcp, dum2d)
 call netcdf_error(error, 'WRITING TPRCP RECORD' )

 dum2d = srflag
 error = nf90_put_var( ncid, id_srflag, dum2d)
 call netcdf_error(error, 'WRITING SRFLAG RECORD' )

 dum2d = reshape(sfcoutput%snow_depth, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_snwdph, dum2d)
 call netcdf_error(error, 'WRITING SNWDPH RECORD' )

 dum2d = reshape(sfcoutput%greenfrc_min, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_shdmin, dum2d)
 call netcdf_error(error, 'WRITING SHDMIN RECORD' )

 dum2d = reshape(sfcoutput%greenfrc_max, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_shdmax, dum2d)
 call netcdf_error(error, 'WRITING SHDMAX RECORD' )

 dum2d = reshape(float(sfcoutput%slope_type), (/lonb,latb/) )
 error = nf90_put_var( ncid, id_slope, dum2d)
 call netcdf_error(error, 'WRITING SLOPE RECORD' )

 dum2d = reshape(sfcoutput%mxsnow_alb, (/lonb,latb/) )
 error = nf90_put_var( ncid, id_snoalb, dum2d)
 call netcdf_error(error, 'WRITING SNOALB RECORD' )

 deallocate (dum2d)

 allocate(dum3d(lonb,latb,lsoil))

 dum3d = reshape(sfcoutput%soil_temp, (/lonb,latb,lsoil/) )
 error = nf90_put_var( ncid, id_stc, dum3d)
 call netcdf_error(error, 'WRITING STC RECORD' )

 dum3d = reshape(sfcoutput%soilm_tot, (/lonb,latb,lsoil/) )
 error = nf90_put_var( ncid, id_smc, dum3d)
 call netcdf_error(error, 'WRITING SMC RECORD' )

 dum3d = reshape(sfcoutput%soilm_liq, (/lonb,latb,lsoil/) )
 error = nf90_put_var( ncid, id_slc, dum3d)
 call netcdf_error(error, 'WRITING SLC RECORD' )

 deallocate (dum3d)

 if (write_nsst) then

   allocate(dum2d(lonb,latb))

   dum2d = reshape(nsst_output(:,17), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_tref, dum2d)
   call netcdf_error(error, 'WRITING TREF RECORD' )

   dum2d = reshape(nsst_output(:,10), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_z_c, dum2d)
   call netcdf_error(error, 'WRITING Z_C RECORD' )

   dum2d = reshape(nsst_output(:,11), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_c_0, dum2d)
   call netcdf_error(error, 'WRITING C_0 RECORD' )

   dum2d = reshape(nsst_output(:,12), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_c_d, dum2d)
   call netcdf_error(error, 'WRITING C_D RECORD' )

   dum2d = reshape(nsst_output(:,13), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_w_0, dum2d)
   call netcdf_error(error, 'WRITING W_0 RECORD' )

   dum2d = reshape(nsst_output(:,14), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_w_d, dum2d)
   call netcdf_error(error, 'WRITING W_D RECORD' )

   dum2d = reshape(nsst_output(:,1), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xt, dum2d)
   call netcdf_error(error, 'WRITING XT RECORD' )

   dum2d = reshape(nsst_output(:,2), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xs, dum2d)
   call netcdf_error(error, 'WRITING XS RECORD' )

   dum2d = reshape(nsst_output(:,3), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xu, dum2d)
   call netcdf_error(error, 'WRITING XU RECORD' )

   dum2d = reshape(nsst_output(:,4), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xv, dum2d)
   call netcdf_error(error, 'WRITING XV RECORD' )

   dum2d = reshape(nsst_output(:,5), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xz, dum2d)
   call netcdf_error(error, 'WRITING XZ RECORD' )

   dum2d = reshape(nsst_output(:,6), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_zm, dum2d)
   call netcdf_error(error, 'WRITING ZM RECORD' )

   dum2d = reshape(nsst_output(:,7), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xtts, dum2d)
   call netcdf_error(error, 'WRITING XTTS RECORD' )

   dum2d = reshape(nsst_output(:,8), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_xzts, dum2d)
   call netcdf_error(error, 'WRITING XZTS RECORD' )

   dum2d = reshape(nsst_output(:,15), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_d_conv, dum2d)
   call netcdf_error(error, 'WRITING D_CONV RECORD' )

   dum2d = reshape(nsst_output(:,16), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_ifd, dum2d)
   call netcdf_error(error, 'WRITING IFD RECORD' )

   dum2d = reshape(nsst_output(:,9), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_dt_cool, dum2d)
   call netcdf_error(error, 'WRITING DT_COOL RECORD' )

   dum2d = reshape(nsst_output(:,18), (/lonb,latb/) )
   error = nf90_put_var(ncid, id_qrain, dum2d)
   call netcdf_error(error, 'WRITING QRAIN RECORD' )

   deallocate(dum2d)

 endif

 error = nf90_close(ncid)

 end subroutine write_fv3_sfc_data_netcdf

 SUBROUTINE READ_FV3_LATLON_NETCDF(TILE_NUM, IMO, JMO, GEOLON, GEOLAT)

 use netcdf

 IMPLICIT NONE

 INTEGER, INTENT(IN)     :: TILE_NUM, IMO, JMO

 REAL, INTENT(OUT)       :: GEOLON(IMO,JMO), GEOLAT(IMO,JMO)

 CHARACTER(LEN=256)      :: TILEFILE

 INTEGER                 :: ERROR, ID_DIM, NCID, NX, NY
 INTEGER                 :: ID_VAR

 REAL, ALLOCATABLE       :: TMPVAR(:,:)

 WRITE(TILEFILE, "(A,I1)") "chgres.fv3.grd.t", TILE_NUM

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERROR(ERROR, 'OPENING FILE: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX' )

 ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY' )

 IF ((NX/2) /= IMO .OR. (NY/2) /= JMO) THEN
   PRINT*,'FATAL ERROR: DIMENSIONS IN GRID FILE WRONG.'
   CALL ERREXIT(160)
 ENDIF

 ALLOCATE(TMPVAR(NX,NY))

 ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X RECORD' )

 GEOLON(1:IMO,1:JMO)     = TMPVAR(2:NX:2,2:NY:2)

 ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y RECORD' )

 GEOLAT(1:IMO,1:JMO)     = TMPVAR(2:NX:2,2:NY:2)

 DEALLOCATE(TMPVAR)

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_FV3_LATLON_NETCDF

 SUBROUTINE READ_FV3_GRID_DIMS_NETCDF(TILE_NUM,IMO,JMO)
   
 use netcdf

 IMPLICIT NONE

 INTEGER, INTENT(IN)   :: TILE_NUM
 INTEGER, INTENT(OUT)  :: IMO, JMO

 CHARACTER(LEN=256)    :: TILEFILE

 INTEGER               :: ERROR, NCID, ID_DIM
  
 IF (TILE_NUM < 10) THEN
   WRITE(TILEFILE, "(A,I1)") "chgres.fv3.orog.t", TILE_NUM
 ELSE
   WRITE(TILEFILE, "(A,I2)") "chgres.fv3.orog.t", TILE_NUM
 ENDIF

 PRINT*,'WILL READ GRID DIMENSIONS FROM: ', TRIM(TILEFILE)

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERROR(ERROR, 'OPENING: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'READING LON ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IMO)
 CALL NETCDF_ERROR(ERROR, 'READING LON VALUE' )

 PRINT*,'I-DIRECTION GRID DIM: ',IMO

 ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'READING LAT ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JMO)
 CALL NETCDF_ERROR(ERROR, 'READING LAT VALUE' )

 PRINT*,'J-DIRECTION GRID DIM: ',JMO

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_FV3_GRID_DIMS_NETCDF

 SUBROUTINE READ_FV3_GRID_DATA_NETCDF(FIELD,TILE_NUM,IMO,JMO,SFCDATA)
   
 use netcdf

 IMPLICIT NONE

 CHARACTER(LEN=*)      :: FIELD

 INTEGER, INTENT(IN)   :: IMO, JMO, TILE_NUM

 REAL, INTENT(OUT)     :: SFCDATA(IMO,JMO)

 CHARACTER(LEN=256)    :: TILEFILE

 INTEGER               :: ERROR, NCID, LAT, LON, ID_DIM
 INTEGER               :: ID_VAR
  
 IF (TILE_NUM < 10) THEN
   WRITE(TILEFILE, "(A,I1)") "chgres.fv3.orog.t", TILE_NUM
 ELSE
   WRITE(TILEFILE, "(A,I2)") "chgres.fv3.orog.t", TILE_NUM
 ENDIF

 PRINT*,'WILL READ ',TRIM(FIELD), ' FROM: ', TRIM(TILEFILE)

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERROR(ERROR, 'OPENING: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'READING LON ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=LON)
 CALL NETCDF_ERROR(ERROR, 'READING LON VALUE' )

 PRINT*,'LON IS ',LON
 IF(LON/=IMO) THEN
   PRINT*,'FATAL ERROR: I-DIMENSIONS DO NOT MATCH ',LON,IMO
   CALL ERREXIT(101)
 ENDIF

 ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'READING LAT ID' )
 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=LAT)
 CALL NETCDF_ERROR(ERROR, 'READING LAT VALUE' )

 PRINT*,'LAT IS ',LAT
 IF(LAT/=JMO) THEN
   PRINT*,'FATAL ERROR: J-DIMENSIONS DO NOT MATCH ',LAT,JMO
   CALL ERREXIT(102)
 ENDIF

 ERROR=NF90_INQ_VARID(NCID, FIELD, ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'READING FIELD ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, SFCDATA)
 CALL NETCDF_ERROR(ERROR, 'READING FIELD' )

 ERROR = NF90_CLOSE(NCID)

 END SUBROUTINE READ_FV3_GRID_DATA_NETCDF

 SUBROUTINE WRITE_FV3_ATMS_BNDY_NETCDF(ZS,PS,T,W,U,V,Q,VCOORD,LONB,LATB,&
                                  LEVSO,NTRACM,NVCOORD,HALO,INPTYP, &
                                  MODELNAME)

!---------------------------------------------------------------------------
!
! Output data along the four halo boundaries.  The naming convention is
! based on point (1,1) being in the lower left corner of the grid:
!
!          --------------- TOP ---------------
!          |                                 |
!          |                                 |
!     LEFT |                                 | RIGHT
!          |                                 |
!          |PT(1,1)                          |
!          ------------- BOTTOM --------------
!
!---------------------------------------------------------------------------

 use netcdf

 IMPLICIT NONE

 CHARACTER(LEN=8), INTENT(IN) :: MODELNAME

 INTEGER,  INTENT(IN)  :: LONB, LATB, LEVSO, NTRACM
 INTEGER,  INTENT(IN)  :: NVCOORD, HALO, INPTYP

 REAL, INTENT(IN)      :: PS(LONB,LATB), ZS(LONB,LATB)
 REAL, INTENT(IN)      :: T(LONB,LATB,LEVSO), W(LONB,LATB,LEVSO)
 REAL, INTENT(IN)      :: U(LONB,LATB,LEVSO), V(LONB,LATB,LEVSO)
 REAL, INTENT(IN)      :: Q(LONB,LATB,LEVSO,NTRACM)
 REAL, INTENT(IN)      :: VCOORD(LEVSO+1,NVCOORD)

 CHARACTER(LEN=256)    :: OUTFILE, TILEFILE

 INTEGER               :: I, II, J, JJ, IHALO, JHALO, K
 INTEGER               :: HALO_P1, IM, JM, JM2, ID_VAR
 INTEGER               :: ID_I_BOTTOM, ID_J_BOTTOM
 INTEGER               :: ID_I_TOP, ID_J_TOP
 INTEGER               :: ID_I_RIGHT, ID_J_RIGHT
 INTEGER               :: ID_I_LEFT, ID_J_LEFT
 INTEGER               :: ID_I_W_BOTTOM, ID_J_W_BOTTOM
 INTEGER               :: ID_I_W_TOP, ID_J_W_TOP
 INTEGER               :: ID_I_W_RIGHT, ID_J_W_RIGHT
 INTEGER               :: ID_I_W_LEFT, ID_J_W_LEFT
 INTEGER               :: ID_I_S_BOTTOM, ID_J_S_BOTTOM
 INTEGER               :: ID_I_S_TOP, ID_J_S_TOP
 INTEGER               :: ID_I_S_RIGHT, ID_J_S_RIGHT
 INTEGER               :: ID_I_S_LEFT, ID_J_S_LEFT
 INTEGER               :: ID_PS_TOP, ID_PS_BOTTOM
 INTEGER               :: ID_PS_RIGHT, ID_PS_LEFT
 INTEGER               :: ID_T_TOP, ID_T_BOTTOM
 INTEGER               :: ID_T_RIGHT, ID_T_LEFT
 INTEGER               :: ID_SPHUM_TOP, ID_SPHUM_BOTTOM
 INTEGER               :: ID_SPHUM_RIGHT, ID_SPHUM_LEFT
 INTEGER               :: ID_CLWMR_TOP, ID_CLWMR_BOTTOM
 INTEGER               :: ID_CLWMR_RIGHT, ID_CLWMR_LEFT
 INTEGER               :: ID_O3MR_TOP, ID_O3MR_BOTTOM
 INTEGER               :: ID_O3MR_RIGHT, ID_O3MR_LEFT
 INTEGER               :: ID_RWMR_TOP, ID_RWMR_BOTTOM
 INTEGER               :: ID_RWMR_RIGHT, ID_RWMR_LEFT
 INTEGER               :: ID_ICMR_TOP, ID_ICMR_BOTTOM
 INTEGER               :: ID_ICMR_RIGHT, ID_ICMR_LEFT
 INTEGER               :: ID_SNMR_TOP, ID_SNMR_BOTTOM
 INTEGER               :: ID_SNMR_RIGHT, ID_SNMR_LEFT
 INTEGER               :: ID_GRLE_TOP, ID_GRLE_BOTTOM
 INTEGER               :: ID_GRLE_RIGHT, ID_GRLE_LEFT
 INTEGER               :: ID_W_TOP, ID_W_BOTTOM
 INTEGER               :: ID_W_RIGHT, ID_W_LEFT
 INTEGER               :: ID_ZH_TOP, ID_ZH_BOTTOM
 INTEGER               :: ID_ZH_RIGHT, ID_ZH_LEFT
 INTEGER               :: ID_U_S_TOP, ID_U_S_BOTTOM
 INTEGER               :: ID_U_S_RIGHT, ID_U_S_LEFT
 INTEGER               :: ID_U_W_TOP, ID_U_W_BOTTOM
 INTEGER               :: ID_U_W_RIGHT, ID_U_W_LEFT
 INTEGER               :: ID_V_S_TOP, ID_V_S_BOTTOM
 INTEGER               :: ID_V_S_RIGHT, ID_V_S_LEFT
 INTEGER               :: ID_V_W_TOP, ID_V_W_BOTTOM
 INTEGER               :: ID_V_W_RIGHT, ID_V_W_LEFT
 INTEGER               :: ERROR, ID_DIM, NX, NY, NCID, NCID2
 INTEGER               :: DIM_HALO, DIM_HALOP, DIM_LON, DIM_LAT
 INTEGER               :: DIM_LATM, DIM_LONP
 INTEGER               :: DIM_LEV, DIM_LEVP
 INTEGER               :: ISTART, IEND, JSTART, JEND
 INTEGER               :: LEVSO_P1
 INTEGER               :: INITAL=0, FSIZE=65536
 INTEGER               :: HEADER_BUFFER_VAL = 16384
 INTEGER, ALLOCATABLE  :: IDUM(:)

 REAL, ALLOCATABLE     :: GEOLAT(:,:), GEOLON(:,:)
 REAL, ALLOCATABLE     :: GEOLAT_HALO(:,:), GEOLON_HALO(:,:)
 REAL, ALLOCATABLE     :: HALO_2D(:,:), HALO_3D(:,:,:), HALO_3D2(:,:,:)
 REAL, ALLOCATABLE     :: AK(:), BK(:), ZH(:,:,:)

 REAL(KIND=4), ALLOCATABLE  :: HALO_2D_4BYTE(:,:)
 REAL(KIND=4), ALLOCATABLE  :: HALO_3D_4BYTE(:,:,:)

 print*,''
 print*,'- COMPUTE AND OUTPUT LATERAL BOUNDARY DATA.'

 HALO_P1 = HALO + 1

 LEVSO_P1 = LEVSO + 1

 ALLOCATE(AK(LEVSO_P1))
 ALLOCATE(BK(LEVSO_P1))
 ALLOCATE(ZH(LONB,LATB,LEVSO_P1))

 AK = VCOORD(:,1)
 BK = VCOORD(:,2)

 CALL COMPUTE_ZH(LONB,LATB,LEVSO,AK,BK,PS,ZS,T,Q,ZH)
    
 DEALLOCATE(AK, BK)

!----------------------------------------------------------------------------------
! Read FV3 grid file.  This routine only works for a regional domain
! and assumes that domain is tile number 7.
!----------------------------------------------------------------------------------

 TILEFILE="chgres.fv3.grd.t7"

 PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(TILEFILE)

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERROR(ERROR, 'OPENING FILE: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX' )

 ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY' )

 IF (MOD(NX,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NX IS NOT EVEN'
   CALL ERREXIT(130)
 ENDIF

 IF (MOD(NY,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NY IS NOT EVEN'
   CALL ERREXIT(131)
 ENDIF

 IM = NX/2
 JM = NY/2

 ALLOCATE(GEOLON(NX+1,NY+1))
 ALLOCATE(GEOLAT(NX+1,NY+1))

 ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X RECORD' )

 ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y RECORD' )

 ERROR = NF90_CLOSE(NCID2)

!----------------------------------------------------------------------------------
! Create output file header.
!----------------------------------------------------------------------------------

 WRITE(OUTFILE, '(A, I1, A)'), 'gfs_bndy.tile', 7, '.nc'
 ERROR = NF90_CREATE(OUTFILE, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                     NCID2, INITIALSIZE=INITAL, CHUNKSIZE=FSIZE)
 CALL NETCDF_ERROR(ERROR, 'CREATING FILE: '//TRIM(OUTFILE) )

 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'FV3GFS GAUSSIAN NEMSIO FILE')
 ELSEIF (INPTYP == 1) THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'GFS NEMSIO FILE')
 ELSEIF (INPTYP == 2) THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'GFS SIGIO FILE')
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING GLOBAL SOURCE ATTRIBUTE')

 ERROR = NF90_DEF_DIM(NCID2, 'lon', IM, DIM_LON)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LON DIMENSION')

 JM2 = JM - (2*HALO)
 ERROR = NF90_DEF_DIM(NCID2, 'lat', JM2, DIM_LAT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LAT DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'lonp', (IM+1), DIM_LONP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LONP DIMENSION')

 JM2 = (JM + 1) - (2*HALO_P1)
 ERROR = NF90_DEF_DIM(NCID2, 'latm', JM2, DIM_LATM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LATM DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'halo', HALO, DIM_HALO)
 CALL NETCDF_ERROR(ERROR, 'DEFINING HALO DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'halop', HALO_P1, DIM_HALOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING HALOP DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'lev', LEVSO, DIM_LEV)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LEV DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'levp', LEVSO_P1, DIM_LEVP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LEVP DIMENSION')

 ERROR = NF90_DEF_VAR(NCID2, 'i_bottom', NF90_INT, &
                             (/DIM_LON/), ID_I_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_BOTTOM, "long_name", "i-indices bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_bottom', NF90_INT, &
                             (/DIM_HALO/), ID_J_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_BOTTOM, "long_name", "j-indices bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_top', NF90_INT, &
                             (/DIM_LON/), ID_I_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_TOP, "long_name", "i-indices top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_top', NF90_INT, &
                             (/DIM_HALO/), ID_J_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_TOP, "long_name", "j-indices top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_right', NF90_INT, &
                             (/DIM_HALO/), ID_I_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_RIGHT, "long_name", "i-indices right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_right', NF90_INT, &
                             (/DIM_LAT/), ID_J_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_RIGHT, "long_name", "j-indices right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_left', NF90_INT, &
                             (/DIM_HALO/), ID_I_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_LEFT, "long_name", "i-indices left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_left', NF90_INT, &
                             (/DIM_LAT/), ID_J_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_LEFT, "long_name", "j-indices left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'ps_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO/), ID_PS_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_BOTTOM, "long_name", "surface pressure bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_BOTTOM, "units", "Pa")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'ps_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO/), ID_PS_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_TOP, "long_name", "surface pressure top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_TOP, "units", "Pa")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'ps_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT/), ID_PS_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_RIGHT, "long_name", "surface pressure right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_RIGHT, "units", "Pa")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'ps_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT/), ID_PS_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_LEFT, "long_name", "surface pressure left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS_LEFT, "units", "Pa")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'w_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_W_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_BOTTOM')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_BOTTOM, "long_name", "vertical velocity bottom bndy")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_BOTTOM, "long_name", "omega bottom bndy")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_BOTTOM ATTRIBUTE')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_BOTTOM, "units", "m/s")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_BOTTOM, "units", "Pa/s")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'w_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_W_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_TOP')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_TOP, "long_name", "vertical velocity top bndy")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_TOP, "long_name", "omega top bndy")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_TOP ATTRIBUTE')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_TOP, "units", "m/s")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_TOP, "units", "Pa/s")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'w_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_W_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_RIGHT')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_RIGHT, "long_name", "vertical velocity right bndy")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_RIGHT, "long_name", "omega right bndy")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_RIGHT ATTRIBUTE')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_RIGHT, "units", "m/s")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_RIGHT, "units", "Pa/s")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'w_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_W_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_LEFT')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_LEFT, "long_name", "vertical velocity left bndy")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_LEFT, "long_name", "omega left bndy")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_LEFT ATTRIBUTE')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W_LEFT, "units", "m/s")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W_LEFT, "units", "Pa/s")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'zh_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEVP/), ID_ZH_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_BOTTOM, "long_name", "height bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_BOTTOM, "units", "m")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'zh_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEVP/), ID_ZH_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_TOP, "long_name", "height top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_TOP, "units", "m")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'zh_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEVP/), ID_ZH_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_RIGHT, "long_name", "height right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_RIGHT, "units", "m")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'zh_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEVP/), ID_ZH_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_LEFT, "long_name", "height left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH_LEFT, "units", "m")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 't_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_T_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_BOTTOM, "long_name", "temperature bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_BOTTOM, "units", "kelvin")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 't_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_T_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_TOP, "long_name", "temperature top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_TOP, "units", "kelvin")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 't_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_T_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_RIGHT, "long_name", "temperature right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_RIGHT, "units", "kelvin")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 't_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_T_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_LEFT, "long_name", "temperature left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_T_LEFT, "units", "kelvin")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'sphum_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_SPHUM_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_BOTTOM, "long_name", "specific humidity bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_BOTTOM, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'sphum_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_SPHUM_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_TOP, "long_name", "specific humidity top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_TOP, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'sphum_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_SPHUM_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_RIGHT, "long_name", "specific humidity right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_RIGHT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'sphum_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_SPHUM_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_LEFT, "long_name", "specific humidity left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM_LEFT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'o3mr_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_O3MR_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_BOTTOM, "long_name", "ozone bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_BOTTOM, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'o3mr_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_O3MR_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_TOP, "long_name", "ozone top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_TOP, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'o3mr_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_O3MR_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_RIGHT, "long_name", "ozone right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_RIGHT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'o3mr_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_O3MR_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_LEFT, "long_name", "ozone left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR_LEFT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'liq_wat_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_CLWMR_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LIQ_WAT_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_BOTTOM, "long_name", "cloud liq water mixing ratio bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_BOTTOM, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'liq_wat_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_CLWMR_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LIQ_WAT_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_TOP, "long_name", "cloud liq water mixing ratio top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_TOP, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'liq_wat_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_CLWMR_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LIQ_WAT_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_RIGHT, "long_name", "cloud liq water mixing ratio right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_RIGHT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'liq_wat_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_CLWMR_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LIQ_WAT_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_LEFT, "long_name", "cloud liq water mixing ratio left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR_LEFT, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR_LEFT UNITS')

 IF (NTRACM > 3) THEN

   ERROR = NF90_DEF_VAR(NCID2, 'rainwat_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_RWMR_BOTTOM)
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_BOTTOM')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_BOTTOM, "long_name", "rain water mixing ratio bottom bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_BOTTOM ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_BOTTOM, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_BOTTOM UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'rainwat_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_RWMR_TOP)
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_TOP')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_TOP, "long_name", "rain water mixing ratio top bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_TOP ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_TOP, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_TOP UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'rainwat_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_RWMR_RIGHT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_RIGHT')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_RIGHT, "long_name", "rain water mixing ratio right bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_RIGHT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_RIGHT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_RIGHT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'rainwat_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_RWMR_LEFT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_LEFT')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_LEFT, "long_name", "rain water mixing ratio left bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_LEFT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR_LEFT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR_LEFT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'ice_wat_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_ICMR_BOTTOM)
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_BOTTOM')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_BOTTOM, "long_name", "ice water mixing ratio bottom bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_BOTTOM ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_BOTTOM, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_BOTTOM UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'ice_wat_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_ICMR_TOP)
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_TOP')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_TOP, "long_name", "ice water mixing ratio top bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_TOP ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_TOP, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_TOP UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'ice_wat_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_ICMR_RIGHT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_RIGHT')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_RIGHT, "long_name", "ice water mixing ratio right bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_RIGHT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_RIGHT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_RIGHT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'ice_wat_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_ICMR_LEFT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_LEFT')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_LEFT, "long_name", "ice water mixing ratio left bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_LEFT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR_LEFT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR_LEFT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'snowwat_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_SNMR_BOTTOM)
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_BOTTOM')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_BOTTOM, "long_name", "snow water mixing ratio bottom bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_BOTTOM ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_BOTTOM, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_BOTTOM UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'snowwat_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_SNMR_TOP)
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_TOP')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_TOP, "long_name", "snow water mixing ratio top bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_TOP ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_TOP, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_TOP UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'snowwat_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_SNMR_RIGHT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_RIGHT')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_RIGHT, "long_name", "snow water mixing ratio right bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_RIGHT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_RIGHT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_RIGHT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'snowwat_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_SNMR_LEFT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_LEFT')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_LEFT, "long_name", "snow water mixing ratio left bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_LEFT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR_LEFT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR_LEFT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'graupel_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_GRLE_BOTTOM)
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_BOTTOM')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_BOTTOM, "long_name", "graupel mixing ratio bottom bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_BOTTOM ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_BOTTOM, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_BOTTOM UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'graupel_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALO, DIM_LEV/), ID_GRLE_TOP)
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_TOP')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_TOP, "long_name", "graupel mixing ratio top bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_TOP ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_TOP, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_TOP UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'graupel_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_GRLE_RIGHT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_RIGHT')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_RIGHT, "long_name", "graupel mixing ratio right bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_RIGHT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_RIGHT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_RIGHT UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'graupel_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LAT, DIM_LEV/), ID_GRLE_LEFT)
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_LEFT')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_LEFT, "long_name", "graupel mixing ratio left bndy")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_LEFT ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE_LEFT, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE_LEFT UNITS')

 ENDIF

 ERROR = NF90_DEF_VAR(NCID2, 'i_w_bottom', NF90_INT, &
                             (/DIM_LONP/), ID_I_W_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_W_BOTTOM, "long_name", "i-indices west edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_w_bottom', NF90_INT, &
                             (/DIM_HALO/), ID_J_W_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_W_BOTTOM, "long_name", "j-indices west edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_w_top', NF90_INT, &
                             (/DIM_LONP/), ID_I_W_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_W_TOP, "long_name", "i-indices west edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_w_top', NF90_INT, &
                             (/DIM_HALO/), ID_J_W_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_W_TOP, "long_name", "j-indices west edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_w_right', NF90_INT, &
                             (/DIM_HALOP/), ID_I_W_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_W_RIGHT, "long_name", "i-indices west edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_w_right', NF90_INT, &
                             (/DIM_LAT/), ID_J_W_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_W_RIGHT, "long_name", "j-indices west edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_w_left', NF90_INT, &
                             (/DIM_HALOP/), ID_I_W_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_W_LEFT, "long_name", "i-indices west edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_W_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_w_left', NF90_INT, &
                             (/DIM_LAT/), ID_J_W_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_W_LEFT, "long_name", "j-indices west edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_W_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'u_w_bottom', NF90_FLOAT, &
                             (/DIM_LONP, DIM_HALO, DIM_LEV/), ID_U_W_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_BOTTOM, "long_name", "u-component wind west edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_BOTTOM, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_w_top', NF90_FLOAT, &
                             (/DIM_LONP, DIM_HALO, DIM_LEV/), ID_U_W_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_TOP, "long_name", "u-component wind west edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_TOP, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_w_right', NF90_FLOAT, &
                             (/DIM_HALOP, DIM_LAT, DIM_LEV/), ID_U_W_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_RIGHT, "long_name", "u-component wind west edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_RIGHT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_w_left', NF90_FLOAT, &
                             (/DIM_HALOP, DIM_LAT, DIM_LEV/), ID_U_W_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_LEFT, "long_name", "u-component wind west edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W_LEFT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_w_bottom', NF90_FLOAT, &
                             (/DIM_LONP, DIM_HALO, DIM_LEV/), ID_V_W_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_BOTTOM, "long_name", "v-component wind west edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_BOTTOM, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_w_top', NF90_FLOAT, &
                             (/DIM_LONP, DIM_HALO, DIM_LEV/), ID_V_W_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_TOP, "long_name", "v-component wind west edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_TOP, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_w_right', NF90_FLOAT, &
                             (/DIM_HALOP, DIM_LAT, DIM_LEV/), ID_V_W_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_RIGHT, "long_name", "v-component wind west edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_RIGHT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_w_left', NF90_FLOAT, &
                             (/DIM_HALOP, DIM_LAT, DIM_LEV/), ID_V_W_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_LEFT, "long_name", "v-component wind west edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W_LEFT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'i_s_bottom', NF90_INT, &
                             (/DIM_LON/), ID_I_S_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_S_BOTTOM, "long_name", "i-indices south edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_s_bottom', NF90_INT, &
                             (/DIM_HALOP/), ID_J_S_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_S_BOTTOM, "long_name", "j-indices south edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_BOTTOM ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_s_top', NF90_INT, &
                             (/DIM_LON/), ID_I_S_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_S_TOP, "long_name", "i-indices south edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_s_top', NF90_INT, &
                             (/DIM_HALOP/), ID_J_S_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_S_TOP, "long_name", "j-indices south edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_TOP ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_s_right', NF90_INT, &
                             (/DIM_HALO/), ID_I_S_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_S_RIGHT, "long_name", "i-indices south edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_s_right', NF90_INT, &
                             (/DIM_LATM/), ID_J_S_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_S_RIGHT, "long_name", "j-indices south edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_RIGHT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'i_s_left', NF90_INT, &
                             (/DIM_HALO/), ID_I_S_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_I_S_LEFT, "long_name", "i-indices south edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING I_S_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'j_s_left', NF90_INT, &
                             (/DIM_LATM/), ID_J_S_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_J_S_LEFT, "long_name", "j-indices south edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING J_S_LEFT ATTRIBUTE')

 ERROR = NF90_DEF_VAR(NCID2, 'u_s_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALOP, DIM_LEV/), ID_U_S_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_BOTTOM, "long_name", "u-component wind south edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_BOTTOM, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_s_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALOP, DIM_LEV/), ID_U_S_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_TOP, "long_name", "u-component wind south edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_TOP, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_s_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LATM, DIM_LEV/), ID_U_S_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_RIGHT, "long_name", "u-component wind south edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_RIGHT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_s_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LATM, DIM_LEV/), ID_U_S_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_LEFT, "long_name", "u-component wind south edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S_LEFT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S_LEFT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_s_bottom', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALOP, DIM_LEV/), ID_V_S_BOTTOM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_BOTTOM')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_BOTTOM, "long_name", "v-component wind south edge bottom bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_BOTTOM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_BOTTOM, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_BOTTOM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_s_top', NF90_FLOAT, &
                             (/DIM_LON, DIM_HALOP, DIM_LEV/), ID_V_S_TOP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_TOP')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_TOP, "long_name", "v-component wind south edge top bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_TOP ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_TOP, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_TOP UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_s_right', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LATM, DIM_LEV/), ID_V_S_RIGHT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_RIGHT')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_RIGHT, "long_name", "v-component wind south edge right bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_RIGHT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_RIGHT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_RIGHT UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_s_left', NF90_FLOAT, &
                             (/DIM_HALO, DIM_LATM, DIM_LEV/), ID_V_S_LEFT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_LEFT')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_LEFT, "long_name", "v-component wind south edge left bndy")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_LEFT ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S_LEFT, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S_LEFT UNITS')

 ERROR = NF90_ENDDEF(NCID2, HEADER_BUFFER_VAL, 4, 0, 4)
 CALL NETCDF_ERROR(ERROR, 'DEFINING END OF HEADER')

!----------------------------------------------------------------------------------
! "Bottom" boundary.
!----------------------------------------------------------------------------------

 ISTART = 1
 IEND   = IM
 JSTART = 1
 JEND   = HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO
 
 ALLOCATE(HALO_2D(ISTART:IEND,JSTART:JEND))
 ALLOCATE(HALO_2D_4BYTE(ISTART:IEND,JSTART:JEND))

 CALL GL2ANY(0, 1, PS, LONB, LATB, HALO_2D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 HALO_2D_4BYTE = REAL(HALO_2D,4)

 ERROR = NF90_PUT_VAR(NCID2, ID_PS_BOTTOM, HALO_2D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING PS_BOTTOM')

 DEALLOCATE(HALO_2D, HALO_2D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANY(0, LEVSO, Q(:,:,:,1), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_SPHUM_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING SPHUM_BOTTOM')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,2), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_O3MR_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING O3MR_BOTTOM')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,3), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_CLWMR_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING LIQ_WAT_BOTTOM')

 IF (NTRACM > 3) THEN

   CALL GL2ANY(0, LEVSO, Q(:,:,:,4), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_RWMR_BOTTOM, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING RWMR_BOTTOM')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,5), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_ICMR_BOTTOM, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING ICMR_BOTTOM')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,6), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_SNMR_BOTTOM, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING SNMR_BOTTOM')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,7), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_GRLE_BOTTOM, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING GRLE_BOTTOM')

 ENDIF

 CALL GL2ANY(0, LEVSO, T, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_T_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING T_BOTTOM')

 CALL GL2ANY(0, LEVSO, W, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_W_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING W_BOTTOM')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO_P1))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO_P1))

 CALL GL2ANY(0, LEVSO_P1, ZH, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO_P1
   HALO_3D_4BYTE(:,:,LEVSO_P1-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_ZH_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING ZH_BOTTOM')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = IM+1
 JSTART = 1
 JEND   = HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_W_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_W_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_W_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_W_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = (2*I)-1
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_W_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_W_BOTTOM')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_W_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_W_BOTTOM')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = IM
 JSTART = 1
 JEND   = HALO + 1

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_S_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_S_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_S_BOTTOM, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_S_BOTTOM')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = (2*J) - 1
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_S_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_S_BOTTOM')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_S_BOTTOM, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_S_BOTTOM')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

!----------------------------------------------------------------------------------
! "Top" boundary.
!----------------------------------------------------------------------------------

 ISTART = 1
 IEND   = IM
 JSTART = JM - HALO + 1
 JEND   = JM

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_TOP')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_TOP')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_2D(ISTART:IEND,JSTART:JEND))
 ALLOCATE(HALO_2D_4BYTE(ISTART:IEND,JSTART:JEND))

 CALL GL2ANY(0, 1, PS, LONB, LATB, HALO_2D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 HALO_2D_4BYTE = REAL(HALO_2D,4)

 ERROR = NF90_PUT_VAR(NCID2, ID_PS_TOP, HALO_2D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING PS_TOP')

 DEALLOCATE(HALO_2D, HALO_2D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANY(0, LEVSO, Q(:,:,:,1), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_SPHUM_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING SPHUM_TOP')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,2), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_O3MR_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING O3MR_TOP')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,3), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_CLWMR_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING CLWMR_TOP')

 IF (NTRACM > 3) THEN

   CALL GL2ANY(0, LEVSO, Q(:,:,:,4), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_RWMR_TOP, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING RWMR_TOP')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,5), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_ICMR_TOP, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING ICMR_TOP')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,6), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_SNMR_TOP, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING SNMR_TOP')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,7), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_GRLE_TOP, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING GRLE_TOP')

 ENDIF

 CALL GL2ANY(0, LEVSO, T, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_T_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING T_TOP')

 CALL GL2ANY(0, LEVSO, W, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_W_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING W_TOP')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO_P1))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO_P1))

 CALL GL2ANY(0, LEVSO_P1, ZH, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO_P1
   HALO_3D_4BYTE(:,:,LEVSO_P1-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_ZH_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING ZH_TOP')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = IM+1
 JSTART = JM - HALO + 1
 JEND   = JM

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_W_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_W_TOP')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_W_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_W_TOP')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = (2*I)-1
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_W_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_W_TOP')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_W_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_W_TOP')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = IM
 JSTART = JM - HALO + 1
 JEND   = JM + 1

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_S_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_S_TOP')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_S_TOP, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_S_TOP')

 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = (2*J) - 1
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_S_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_S_TOP')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_S_TOP, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_S_TOP')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

!----------------------------------------------------------------------------------
! "Left" boundary.
!----------------------------------------------------------------------------------

 ISTART = 1
 IEND   = HALO
 JSTART = HALO + 1
 JEND   = JM - HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_2D(ISTART:IEND,JSTART:JEND))
 ALLOCATE(HALO_2D_4BYTE(ISTART:IEND,JSTART:JEND))

 CALL GL2ANY(0, 1, PS, LONB, LATB, HALO_2D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 HALO_2D_4BYTE = REAL(HALO_2D,4)

 ERROR = NF90_PUT_VAR(NCID2, ID_PS_LEFT, HALO_2D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING PS_LEFT')

 DEALLOCATE(HALO_2D, HALO_2D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANY(0, LEVSO, Q(:,:,:,1), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_SPHUM_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING SPHUM_LEFT')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,2), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_O3MR_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING O3MR_LEFT')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,3), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_CLWMR_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING CLWMR_LEFT')

 IF (NTRACM > 3) THEN

   CALL GL2ANY(0, LEVSO, Q(:,:,:,4), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_RWMR_LEFT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING RWMR_LEFT')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,5), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_ICMR_LEFT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING ICMR_LEFT')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,6), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_SNMR_LEFT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING SNMR_LEFT')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,7), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_GRLE_LEFT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING GRLE_LEFT')

 ENDIF

 CALL GL2ANY(0, LEVSO, T, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_T_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING T_LEFT')

 CALL GL2ANY(0, LEVSO, W, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_W_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING W_LEFT')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO_P1))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO_P1))

 CALL GL2ANY(0, LEVSO_P1, ZH, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO_P1
   HALO_3D_4BYTE(:,:,LEVSO_P1-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_ZH_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING ZH_LEFT')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = HALO + 1
 JSTART = HALO + 1
 JEND   = JM - HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_W_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_W_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_W_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_W_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = (2*I)-1
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_W_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_W_LEFT')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_W_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_W_LEFT')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = 1
 IEND   = HALO
 JSTART = HALO_P1 + 1
 JEND   = JM + 1 - HALO_P1

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_S_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_S_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_S_LEFT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_S_LEFT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = (2*J) - 1
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_S_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_S_LEFT')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_S_LEFT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_S_LEFT')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

!----------------------------------------------------------------------------------
! "Right" boundary.
!----------------------------------------------------------------------------------

 ISTART = IM - HALO + 1
 IEND   = IM
 JSTART = HALO + 1
 JEND   = JM - HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_2D(ISTART:IEND,JSTART:JEND))
 ALLOCATE(HALO_2D_4BYTE(ISTART:IEND,JSTART:JEND))

 CALL GL2ANY(0, 1, PS, LONB, LATB, HALO_2D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 HALO_2D_4BYTE = REAL(HALO_2D,4)

 ERROR = NF90_PUT_VAR(NCID2, ID_PS_RIGHT, HALO_2D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING PS_RIGHT')

 DEALLOCATE(HALO_2D, HALO_2D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANY(0, LEVSO, Q(:,:,:,1), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_SPHUM_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING SPHUM_RIGHT')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,2), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_O3MR_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING O3MR_RIGHT')

 CALL GL2ANY(0, LEVSO, Q(:,:,:,3), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_CLWMR_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING CLWMR_RIGHT')

 IF (NTRACM > 3) THEN

   CALL GL2ANY(0, LEVSO, Q(:,:,:,4), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_RWMR_RIGHT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING RWMR_RIGHT')
   
   CALL GL2ANY(0, LEVSO, Q(:,:,:,5), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_ICMR_RIGHT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING ICMR_RIGHT')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,6), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_SNMR_RIGHT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING SNMR_RIGHT')

   CALL GL2ANY(0, LEVSO, Q(:,:,:,7), LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
   DO K = 1, LEVSO
     HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_GRLE_RIGHT, HALO_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING GRLE_RIGHT')

 ENDIF

 CALL GL2ANY(0, LEVSO, T, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_T_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING T_RIGHT')

 CALL GL2ANY(0, LEVSO, W, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_W_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING W_RIGHT')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO_P1))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO_P1))

 CALL GL2ANY(0, LEVSO_P1, ZH, LONB, LATB, HALO_3D, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)
 DO K = 1, LEVSO_P1
   HALO_3D_4BYTE(:,:,LEVSO_P1-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_ZH_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING ZH_RIGHT')

 DEALLOCATE(HALO_3D, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = IM - HALO + 1
 IEND   = IM + 1
 JSTART = HALO + 1
 JEND   = JM - HALO

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_W_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_W_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_W_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_W_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = (2*I)-1
     JJ = 2*J
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_W_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_W_RIGHT')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_W_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_W_RIGHT')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

 ISTART = IM - HALO + 1
 IEND   = IM
 JSTART = HALO_P1 + 1
 JEND   = JM + 1 - HALO_P1

 IHALO = IEND - ISTART + 1
 JHALO = JEND - JSTART + 1

 ALLOCATE(IDUM(ISTART:IEND))
 DO I = ISTART, IEND
   IDUM(I) = I
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_I_S_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING I_S_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(IDUM(JSTART:JEND))
 DO J = JSTART, JEND
   IDUM(J) = J
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_J_S_RIGHT, IDUM)
 CALL NETCDF_ERROR(ERROR, 'WRITING J_S_RIGHT')
 DEALLOCATE(IDUM)

 ALLOCATE(GEOLAT_HALO(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_HALO(ISTART:IEND,JSTART:JEND))

 DO J = JSTART, JEND
   DO I = ISTART, IEND
     II = 2*I
     JJ = (2*J) - 1
     GEOLON_HALO(I,J) = GEOLON(II,JJ)
     GEOLAT_HALO(I,J) = GEOLAT(II,JJ)
   ENDDO
 ENDDO

 ALLOCATE(HALO_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(HALO_3D2(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0, LEVSO, U, V, LONB, LATB, HALO_3D, HALO_3D2, IHALO, JHALO, GEOLON_HALO, GEOLAT_HALO)

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_S_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_S_RIGHT')

 DO K = 1, LEVSO
   HALO_3D_4BYTE(:,:,LEVSO-K+1) = REAL(HALO_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_S_RIGHT, HALO_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_S_RIGHT')

 DEALLOCATE(HALO_3D, HALO_3D2, HALO_3D_4BYTE)
 DEALLOCATE(GEOLAT_HALO, GEOLON_HALO)

!----------------------------------------------------------------------------------
! Cleanup and close file.
!----------------------------------------------------------------------------------

 DEALLOCATE(GEOLAT)
 DEALLOCATE(GEOLON)

 ERROR = NF90_CLOSE(NCID2)

 END SUBROUTINE WRITE_FV3_ATMS_BNDY_NETCDF

 SUBROUTINE WRITE_FV3_ATMS_NETCDF(ZS,PS,T,W,U,V,Q,VCOORD,LONB,LATB,&
                                  LEVSO,NTRACM,NVCOORD,NTILES,HALO,&
                                  INPTYP,MODELNAME)

 use netcdf

 IMPLICIT NONE

 CHARACTER(LEN=8), INTENT(IN) :: MODELNAME

 INTEGER,  INTENT(IN)  :: NTILES, LONB, LATB, LEVSO, NTRACM
 INTEGER,  INTENT(IN)  :: NVCOORD, HALO, INPTYP

 REAL, INTENT(IN)      :: PS(LONB,LATB), ZS(LONB,LATB)
 REAL, INTENT(IN)      :: T(LONB,LATB,LEVSO), W(LONB,LATB,LEVSO)
 REAL, INTENT(IN)      :: U(LONB,LATB,LEVSO), V(LONB,LATB,LEVSO)
 REAL, INTENT(IN)      :: Q(LONB,LATB,LEVSO,NTRACM)
 REAL, INTENT(IN)      :: VCOORD(LEVSO+1,NVCOORD)

 CHARACTER(LEN=256)    :: TILEFILE, OUTFILE

 INTEGER               :: ID_DIM, ID_VAR, IM, JM
 INTEGER               :: ERROR, N, NCID, NCID2, NX, NY
 INTEGER               :: INITAL=0, FSIZE=65536
 INTEGER               :: HEADER_BUFFER_VAL = 16384
 INTEGER               :: DIM_LON, DIM_LAT, DIM_LONP, DIM_LATP
 INTEGER               :: DIM_LEV, DIM_LEVP, DIM_TRACER
 INTEGER               :: ID_LON, ID_LAT, ID_PS, ID_T
 INTEGER               :: ID_W, ID_ZH, ID_SPHUM, ID_O3MR
 INTEGER               :: ID_CLWMR, ID_U_W, ID_V_W
 INTEGER               :: ID_RWMR, ID_ICMR, ID_SNMR, ID_GRLE
 INTEGER               :: ID_U_S, ID_V_S, K, LEVSO_P1
 INTEGER               :: I, J, II, JJ
 INTEGER               :: ISTART, IEND, JSTART, JEND, IM_OUT, JM_OUT
 INTEGER               :: START_TILE, END_TILE

 REAL, ALLOCATABLE     :: CUBE_2D(:,:), CUBE_3D(:,:,:), CUBE_3D2(:,:,:)
 REAL, ALLOCATABLE     :: AK(:), BK(:), ZH(:,:,:)
 REAL, ALLOCATABLE     :: GEOLAT(:,:), GEOLAT_W(:,:), GEOLAT_S(:,:)
 REAL, ALLOCATABLE     :: GEOLON(:,:), GEOLON_W(:,:)
 REAL, ALLOCATABLE     :: GEOLON_S(:,:), TMPVAR(:,:)

 REAL(KIND=4), ALLOCATABLE  :: CUBE_2D_4BYTE(:,:)
 REAL(KIND=4), ALLOCATABLE  :: CUBE_3D_4BYTE(:,:,:)

 LEVSO_P1 = LEVSO + 1

 CALL WRITE_FV3_ATMS_HEADER_NETCDF(LEVSO_P1, NTRACM, NVCOORD, VCOORD)

 ALLOCATE(AK(LEVSO_P1))
 ALLOCATE(BK(LEVSO_P1))
 ALLOCATE(ZH(LONB,LATB,(LEVSO_P1)))

 AK = VCOORD(:,1)
 BK = VCOORD(:,2)

 CALL COMPUTE_ZH(LONB,LATB,LEVSO,AK,BK,PS,ZS,T,Q,ZH)
    
 DEALLOCATE(AK, BK)

 PRINT*,''

 IF (HALO == 0) THEN  ! NOT A REGIONAL GRID
   START_TILE = 1
   END_TILE   = NTILES
 ELSE                 ! A REGIONAL GRID.  ASSUME IT IS TILE 7.
   START_TILE = 7
   END_TILE   = 7
 ENDIF

 TILE_LOOP : DO N = START_TILE, END_TILE

 PRINT*,'WRITE FV3 ATMOSPHERIC DATA FOR TILE ',N

 IF (N < 10) THEN
   WRITE(TILEFILE, "(A,I1)") "chgres.fv3.grd.t", N
 ELSE
   WRITE(TILEFILE, "(A,I2)") "chgres.fv3.grd.t", N
 ENDIF

 ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
 CALL NETCDF_ERROR(ERROR, 'OPENING FILE: '//TRIM(TILEFILE) )

 ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NX' )

 ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY ID' )

 ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING NY' )

 IF (MOD(NX,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NX IS NOT EVEN'
   CALL ERREXIT(103)
 ENDIF

 IF (MOD(NY,2) /= 0) THEN
   PRINT*,'FATAL ERROR: NY IS NOT EVEN'
   CALL ERREXIT(104)
 ENDIF

 IM = NX/2
 JM = NY/2

 IF (HALO > 0) THEN
   ISTART     = 1 + HALO
   IEND       = IM - HALO
   JSTART     = 1+ HALO
   JEND       = JM - HALO
   PRINT*,''
   PRINT*,"WILL NOT PROCESS HALO REGION."
   PRINT*,"HALO IS ", HALO, " ROWS/COLUMNS"
   PRINT*,"WILL PROCESS I= ", ISTART, " TO ", IEND
   PRINT*,"WILL PROCESS J= ", JSTART, " TO ", JEND
   PRINT*,''
 ELSE
   ISTART = 1
   IEND   = IM
   JSTART = 1
   JEND   = JM
 ENDIF

 IM_OUT = IEND - ISTART + 1
 JM_OUT = JEND - JSTART +1

 PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(TILEFILE)

 ALLOCATE(TMPVAR(NX+1,NY+1))
 ALLOCATE(GEOLON(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLON_W(ISTART:IEND+1,JSTART:JEND))
 ALLOCATE(GEOLON_S(ISTART:IEND,JSTART:JEND+1))

 ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING X RECORD' )

 DO J = JSTART, JEND
 DO I = ISTART, IEND
   II = 2*I
   JJ = 2*J
   GEOLON(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO

 DO J = JSTART, JEND
 DO I = ISTART, IEND+1
   II = (2*I) - 1
   JJ = 2*J
   GEOLON_W(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO

 DO J = JSTART, JEND+1
 DO I = ISTART, IEND
   II = 2*I
   JJ = (2*J) - 1
   GEOLON_S(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO

 ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR) 
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y ID' )
 ERROR=NF90_GET_VAR(NCID, ID_VAR, TMPVAR)
 CALL NETCDF_ERROR(ERROR, 'ERROR READING Y RECORD' )

 ERROR = NF90_CLOSE(NCID)

 ALLOCATE(GEOLAT(ISTART:IEND,JSTART:JEND))
 ALLOCATE(GEOLAT_W(ISTART:IEND+1,JSTART:JEND))
 ALLOCATE(GEOLAT_S(ISTART:IEND,JSTART:JEND+1))

 DO J = JSTART, JEND
 DO I = ISTART, IEND
   II = 2*I
   JJ = 2*J
   GEOLAT(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO

 DO J = JSTART, JEND
 DO I = ISTART, IEND+1
   II = (2*I) - 1
   JJ = 2*J
   GEOLAT_W(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO

 DO J = JSTART, JEND+1
 DO I = ISTART, IEND
   II = 2*I
   JJ = (2*J) - 1
   GEOLAT_S(I,J) = TMPVAR(II,JJ)
 ENDDO
 ENDDO
 
 DEALLOCATE(TMPVAR)

 IF (N < 10) THEN
   WRITE(OUTFILE, "(A,I1,A)") 'gfs_data.tile', N, '.nc' 
 ELSE
   WRITE(OUTFILE, "(A,I2,A)") 'gfs_data.tile', N, '.nc'
 ENDIF

 ERROR = NF90_CREATE(OUTFILE, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                     NCID2, INITIALSIZE=INITAL, CHUNKSIZE=FSIZE)
 CALL NETCDF_ERROR(ERROR, 'CREATING FILE: '//TRIM(OUTFILE) )

 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'FV3GFS GAUSSIAN NEMSIO FILE')
 ELSEIF (INPTYP == 1) THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'GFS NEMSIO FILE')
 ELSEIF (INPTYP == 2) THEN
   ERROR = NF90_PUT_ATT(NCID2, NF90_GLOBAL, 'source', 'GFS SIGIO FILE')
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING GLOBAL SOURCE ATTRIBUTE')

 ERROR = NF90_DEF_DIM(NCID2, 'lon', IM_OUT, DIM_LON)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LON DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'lat', JM_OUT, DIM_LAT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LAT DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'lonp', (IM_OUT+1), DIM_LONP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LONP DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'latp', (JM_OUT+1), DIM_LATP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LATP DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'lev', LEVSO, DIM_LEV)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LEV DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'levp', LEVSO_P1, DIM_LEVP)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LEVP DIMENSION')

 ERROR = NF90_DEF_DIM(NCID2, 'ntracer', NTRACM, DIM_TRACER)
 CALL NETCDF_ERROR(ERROR, 'DEFINING NTRACER DIMENSION')

 ERROR = NF90_DEF_VAR(NCID2, 'lon', NF90_FLOAT, DIM_LON, ID_LON)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LON VARIABLE')

 ERROR = NF90_PUT_ATT(NCID2, ID_LON, "cartesian_axis", "X")
 CALL NETCDF_ERROR(ERROR, 'DEFINING X-AXIS')

 ERROR = NF90_DEF_VAR(NCID2, 'lat', NF90_FLOAT, DIM_LAT, ID_LAT)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LAT VARIABLE')

 ERROR = NF90_PUT_ATT(NCID2, ID_LAT, "cartesian_axis", "Y")
 CALL NETCDF_ERROR(ERROR, 'DEFINING Y-AXIS')

 ERROR = NF90_DEF_VAR(NCID2, 'ps', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT/), ID_PS)
 CALL NETCDF_ERROR(ERROR, 'DEFINING PS')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS, "long_name", "surface pressure")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PRESSURE ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_PS, "units", "Pa")
 CALL NETCDF_ERROR(ERROR, 'DEFINING PRESSURE UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'w', NF90_FLOAT,  &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_W)
 CALL NETCDF_ERROR(ERROR, 'DEFINING W')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W, "long_name", "vertical velocity")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W, "long_name", "omega")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W ATTRIBUTE')
 IF (TRIM(MODELNAME) == "FV3GFS") THEN
   ERROR = NF90_PUT_ATT(NCID2, ID_W, "units", "m/s")
 ELSE
   ERROR = NF90_PUT_ATT(NCID2, ID_W, "units", "Pa/s")
 ENDIF
 CALL NETCDF_ERROR(ERROR, 'DEFINING W UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'zh', NF90_FLOAT,  &
                             (/DIM_LON, DIM_LAT, DIM_LEVP/), ID_ZH)
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH, "long_name", "height")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_ZH, "units", "m")
 CALL NETCDF_ERROR(ERROR, 'DEFINING ZH UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 't', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_T)
 CALL NETCDF_ERROR(ERROR, 'DEFINING t')
 ERROR = NF90_PUT_ATT(NCID2, ID_T, "long_name", "temperature")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_T, "units", "kelvin")
 CALL NETCDF_ERROR(ERROR, 'DEFINING T UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'sphum', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_SPHUM)
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM, "long_name", "specific humidity")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_SPHUM, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING SPHUM UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'o3mr', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_O3MR)
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR, "long_name", "ozone")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_O3MR, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING O3MR UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'liq_wat', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_CLWMR)
 CALL NETCDF_ERROR(ERROR, 'DEFINING LIQ_WAT')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR, "long_name", "cloud liquid water mixing ratio")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_CLWMR, "units", "kg/kg")
 CALL NETCDF_ERROR(ERROR, 'DEFINING CLWMR UNITS')

 IF (NTRACM > 3) THEN

   ERROR = NF90_DEF_VAR(NCID2, 'rainwat', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_RWMR)
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR, "long_name", "rain water mixing ratio")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_RWMR, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING RWMR UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'ice_wat', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_ICMR)
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR, "long_name", "ice water mixing ratio")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_ICMR, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING ICMR UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'snowwat', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_SNMR)
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR, "long_name", "snow water mixing ratio")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_SNMR, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING SNMR UNITS')

   ERROR = NF90_DEF_VAR(NCID2, 'graupel', NF90_FLOAT, &
                             (/DIM_LON, DIM_LAT, DIM_LEV/), ID_GRLE)
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE, "long_name", "graupel mixing ratio")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE ATTRIBUTE')
   ERROR = NF90_PUT_ATT(NCID2, ID_GRLE, "units", "kg/kg")
   CALL NETCDF_ERROR(ERROR, 'DEFINING GRLE UNITS')

 ENDIF

 ERROR = NF90_DEF_VAR(NCID2, 'u_w', NF90_FLOAT, &
                             (/DIM_LONP, DIM_LAT, DIM_LEV/), ID_U_W)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W, "long_name", "u-component wind on d grid")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_W, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_W UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_w', NF90_FLOAT, &
                             (/DIM_LONP, DIM_LAT, DIM_LEV/), ID_V_W)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W, "long_name", "v-component wind on c grid")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_W, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_W UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'u_s', NF90_FLOAT,  &
                             (/DIM_LON, DIM_LATP, DIM_LEV/), ID_U_S)
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S, "long_name", "u-component wind on c grid")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_U_S, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING U_S UNITS')

 ERROR = NF90_DEF_VAR(NCID2, 'v_s', NF90_FLOAT,  &
                             (/DIM_LON, DIM_LATP, DIM_LEV/), ID_V_S)
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S, "long_name", "v-component wind on d grid")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S ATTRIBUTE')
 ERROR = NF90_PUT_ATT(NCID2, ID_V_S, "units", "m/s")
 CALL NETCDF_ERROR(ERROR, 'DEFINING V_S UNITS')

 ERROR = NF90_ENDDEF(NCID2, HEADER_BUFFER_VAL, 4, 0, 4)
 CALL NETCDF_ERROR(ERROR, 'DEFINING END OF HEADER')

!------------------------------------------------------------------
! Write out data.  fv3 convention: lowest model level is levso.
! top model layer is 1.  this is opposite the gfs convention.
!------------------------------------------------------------------

 ALLOCATE(CUBE_2D(ISTART:IEND,JSTART:JEND), CUBE_2D_4BYTE(ISTART:IEND,JSTART:JEND))

 CUBE_2D_4BYTE = REAL(GEOLON,4)
 ERROR = NF90_PUT_VAR(NCID2, ID_LON, CUBE_2D_4BYTE(:,JSTART))
 CALL NETCDF_ERROR(ERROR, 'WRITING LON')

 CUBE_2D_4BYTE = REAL(GEOLAT,4)
 ERROR = NF90_PUT_VAR(NCID2, ID_LAT, CUBE_2D_4BYTE(ISTART,:))
 CALL NETCDF_ERROR(ERROR, 'WRITING LAT')

 CALL GL2ANY(0,1,PS,LONB,LATB,CUBE_2D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 CUBE_2D_4BYTE = REAL(CUBE_2D,4)

 ERROR = NF90_PUT_VAR(NCID2, ID_PS, CUBE_2D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING PS')
    
 DEALLOCATE(CUBE_2D_4BYTE, CUBE_2D)

 ALLOCATE(CUBE_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO_P1))
 ALLOCATE(CUBE_3D(ISTART:IEND,JSTART:JEND,LEVSO_P1))

 CALL GL2ANY(0,LEVSO_P1,ZH,LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 DO K = 1, LEVSO_P1
   CUBE_3D_4BYTE(:,:,LEVSO_P1-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_ZH, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING ZH')

 DEALLOCATE(CUBE_3D, CUBE_3D_4BYTE)

 ALLOCATE(CUBE_3D_4BYTE(ISTART:IEND,JSTART:JEND,LEVSO))
 ALLOCATE(CUBE_3D(ISTART:IEND,JSTART:JEND,LEVSO))

 CALL GL2ANY(0,LEVSO,W,LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_W, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING W')

 CALL GL2ANY(0,LEVSO,T,LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON,GEOLAT)
 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_T, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING T')

 CALL GL2ANY(0,LEVSO,Q(:,:,:,1),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_SPHUM, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING SPHUM')

 CALL GL2ANY(0,LEVSO,Q(:,:,:,2),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_O3MR, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING O3MR')

 CALL GL2ANY(0,LEVSO,Q(:,:,:,3),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_CLWMR, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING CLWMR')

 IF (NTRACM > 3) THEN

   CALL GL2ANY(0,LEVSO,Q(:,:,:,4),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
   DO K = 1, LEVSO
     CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_RWMR, CUBE_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING RWMR')
 
   CALL GL2ANY(0,LEVSO,Q(:,:,:,5),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
   DO K = 1, LEVSO
     CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_ICMR, CUBE_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING ICMR')

   CALL GL2ANY(0,LEVSO,Q(:,:,:,6),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
   DO K = 1, LEVSO
     CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_SNMR, CUBE_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING SNMR')

   CALL GL2ANY(0,LEVSO,Q(:,:,:,7),LONB,LATB,CUBE_3D,IM_OUT,JM_OUT,GEOLON, GEOLAT)
   DO K = 1, LEVSO
     CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
   ENDDO

   ERROR = NF90_PUT_VAR(NCID2, ID_GRLE, CUBE_3D_4BYTE)
   CALL NETCDF_ERROR(ERROR, 'WRITING GRLE')

 ENDIF

 DEALLOCATE (CUBE_3D, CUBE_3D_4BYTE)

 ALLOCATE(CUBE_3D_4BYTE(ISTART:IEND+1,JSTART:JEND,LEVSO))
 ALLOCATE(CUBE_3D(ISTART:IEND+1,JSTART:JEND,LEVSO))
 ALLOCATE(CUBE_3D2(ISTART:IEND+1,JSTART:JEND,LEVSO))

 CALL GL2ANYV(0,LEVSO,U,V,LONB,LATB,CUBE_3D,CUBE_3D2,(IM_OUT+1),JM_OUT,GEOLON_W, GEOLAT_W)

 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_W, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_W')

 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_W, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_W')

 DEALLOCATE (CUBE_3D, CUBE_3D2, CUBE_3D_4BYTE)

 ALLOCATE(CUBE_3D_4BYTE(ISTART:IEND,JSTART:JEND+1,LEVSO))
 ALLOCATE(CUBE_3D(ISTART:IEND,JSTART:JEND+1,LEVSO))
 ALLOCATE(CUBE_3D2(ISTART:IEND,JSTART:JEND+1,LEVSO))

 CALL GL2ANYV(0,LEVSO,U,V,LONB,LATB,CUBE_3D,CUBE_3D2,IM_OUT,(JM_OUT+1),GEOLON_S, GEOLAT_S)

 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_U_S, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING U_S')

 DO K = 1, LEVSO
   CUBE_3D_4BYTE(:,:,LEVSO-K+1) = REAL(CUBE_3D2(:,:,K),4)
 ENDDO

 ERROR = NF90_PUT_VAR(NCID2, ID_V_S, CUBE_3D_4BYTE)
 CALL NETCDF_ERROR(ERROR, 'WRITING V_S')

 ERROR = NF90_CLOSE(NCID2)

 DEALLOCATE(CUBE_3D, CUBE_3D2, CUBE_3D_4BYTE)
 DEALLOCATE(GEOLON, GEOLON_W, GEOLON_S)
 DEALLOCATE(GEOLAT, GEOLAT_W, GEOLAT_S)

 ENDDO TILE_LOOP

 DEALLOCATE(ZH)

 END SUBROUTINE WRITE_FV3_ATMS_NETCDF

 SUBROUTINE READ_GFS_NSST_DATA_NSTIO(IMI,JMI,NUM_NSST_FIELDS,     &
                                     NSST_INPUT, MASK_INPUT, &
                                     NSST_YEAR, NSST_MON,    &
                                     NSST_DAY, NSST_HOUR,    &
                                     NSST_FHOUR)

 USE NSTIO_MODULE

 IMPLICIT NONE

 INTEGER, INTENT(IN)  :: IMI, JMI, NUM_NSST_FIELDS
 INTEGER, INTENT(OUT) :: NSST_YEAR, NSST_MON
 INTEGER, INTENT(OUT) :: NSST_DAY, NSST_HOUR

 REAL, INTENT(OUT)   :: NSST_FHOUR
 REAL, INTENT(OUT)   :: MASK_INPUT(IMI,JMI)
 REAL, INTENT(OUT)   :: NSST_INPUT(IMI,JMI,NUM_NSST_FIELDS)

 INTEGER(NSTIO_INTKIND) :: NSSTI, IRET

 TYPE(NSTIO_HEAD)        :: NSST_IN_HEAD
 TYPE(NSTIO_DATA)        :: NSST_IN_DATA

 PRINT*,'- READ NSST FILE chgres.inp.nst.'
!  OPEN NSST FILES
 NSSTI=31
 CALL NSTIO_SROPEN(NSSTI,'chgres.inp.nst',IRET)
 IF(IRET/=0)THEN
   PRINT*,'FATAL ERROR OPENING chgres.inp.nst ', IRET
   CALL ERREXIT(105)
 ENDIF
 CALL NSTIO_SRHEAD(NSSTI,NSST_IN_HEAD,IRET)
 IF(IRET/=0)THEN
   PRINT*,'FATAL ERROR READING chgres.inp.nst ', IRET
   CALL ERREXIT(106)
 ENDIF
 CALL NSTIO_ALDATA(NSST_IN_HEAD,NSST_IN_DATA,IRET)
 CALL NSTIO_SRDATA(NSSTI,NSST_IN_HEAD,NSST_IN_DATA,IRET)
 IF(IRET/=0)THEN
   PRINT*,'FATAL ERROR READING chgres.inp.nst ', IRET
   CALL ERREXIT(107)
 ENDIF
 NSST_YEAR=NSST_IN_HEAD%IDATE(4)
 NSST_MON=NSST_IN_HEAD%IDATE(2)
 NSST_DAY=NSST_IN_HEAD%IDATE(3)
 NSST_HOUR=NSST_IN_HEAD%IDATE(1)
 NSST_FHOUR=NSST_IN_HEAD%FHOUR
 NSST_INPUT(:,:,1)=NSST_IN_DATA%XT
 NSST_INPUT(:,:,2)=NSST_IN_DATA%XS
 NSST_INPUT(:,:,3)=NSST_IN_DATA%XU
 NSST_INPUT(:,:,4)=NSST_IN_DATA%XV
 NSST_INPUT(:,:,5)=NSST_IN_DATA%XZ
 NSST_INPUT(:,:,6)=NSST_IN_DATA%ZM
 NSST_INPUT(:,:,7)=NSST_IN_DATA%XTTS
 NSST_INPUT(:,:,8)=NSST_IN_DATA%XZTS
 NSST_INPUT(:,:,9)=NSST_IN_DATA%DT_COOL
 NSST_INPUT(:,:,10)=NSST_IN_DATA%Z_C
 NSST_INPUT(:,:,11)=NSST_IN_DATA%C_0
 NSST_INPUT(:,:,12)=NSST_IN_DATA%C_D
 NSST_INPUT(:,:,13)=NSST_IN_DATA%W_0
 NSST_INPUT(:,:,14)=NSST_IN_DATA%W_D
 NSST_INPUT(:,:,15)=NSST_IN_DATA%D_CONV
 NSST_INPUT(:,:,16)=NSST_IN_DATA%IFD
 NSST_INPUT(:,:,17)=NSST_IN_DATA%TREF
 NSST_INPUT(:,:,18)=NSST_IN_DATA%QRAIN
 MASK_INPUT=NSST_IN_DATA%SLMSK
 CALL NSTIO_AXDATA(NSST_IN_DATA,IRET)
 CALL NSTIO_SRCLOSE(NSSTI,IRET)

 END SUBROUTINE READ_GFS_NSST_DATA_NSTIO

 SUBROUTINE READ_FV3GFS_NSST_DATA_NEMSIO (MASK_INPUT,NSST_INPUT,IMI,JMI, &
                 NUM_NSST_FIELDS,NSST_YEAR,NSST_MON,NSST_DAY,    & 
                 NSST_HOUR,NSST_FHOUR)

!-----------------------------------------------------------------------
! Subroutine: read nsst data from a fv3gfs nemsio file
!
! Author: George Gayno/EMC
!
! Abstract: Reads an fv3gfs nsst file in nemsio format.  Places data
!           in the "nsst_input" array in the order expected by routine
!           nsst_chgres.
!
! Input files: 
!    "chgres.inp.sfc" - input nsst nemsio file.  note: fv3gfs sfc
!                       and nsst fields are in the same file.
!
! Output files:  none
!
! History:
!   2018-05-31   Gayno - Initial version
!
! Condition codes:  all non-zero codes are fatal
!   109 - bad open of nst file "chgres.inp.sfc"
!   110 - bad read of "chgres.inp.sfc" header
!   112 - wrong number of nsst records
!   113 - bad read of landmask record.
!   114 - bad read of an nst file record.
!-----------------------------------------------------------------------

 use nemsio_module

 implicit none

 integer, parameter      :: nrec=18

 character(len=3)        :: levtyp
 character(len=8)        :: recname(nrec)

 integer, intent(in)     :: imi, jmi, num_nsst_fields
 integer, intent(out)    :: nsst_year, nsst_mon
 integer, intent(out)    :: nsst_day, nsst_hour

 real,    intent(out)    :: mask_input(imi,jmi)
 real,    intent(out)    :: nsst_input(imi,jmi,num_nsst_fields)
 real,    intent(out)    :: nsst_fhour

 integer(nemsio_intkind) :: iret, lev, nframe
 integer(nemsio_intkind) :: idate(7), nfhour

 integer                 :: j

 real(nemsio_realkind),allocatable :: dummy(:)

 type(nemsio_gfile)      :: gfile

 data recname   /"xt      ", "xs      ", "xu      ", &
                 "xv      ", "xz      ", "zm      ", &
                 "xtts    ", "xzts    ", "dtcool  ", &
                 "zc      ", "c0      ", "cd      ", &
                 "w0      ", "wd      ", "dconv   ", &
                 "ifd     ", "tref    ", "qrain   " /

 print*,"- READ INPUT NSST DATA IN NEMSIO FORMAT"

 if (nrec /= num_nsst_fields) then
   print*,"- FATAL ERROR: bad number of nsst records."
   call errexit(112)
 endif

! note: fv3gfs surface and nsst fields are in a single file.
  
 call nemsio_open(gfile, "chgres.inp.sfc", "read", iret=iret)
 if (iret /= 0) then
   print*,"- FATAL ERROR: bad open of chgres.inp.sfc."
   print*,"- IRET IS ", iret
   call errexit(109)
 endif

 print*,"- READ FILE HEADER"
 call nemsio_getfilehead(gfile,iret=iret, &
           idate=idate,nfhour=nfhour)
 if (iret /= 0) then
   print*,"- FATAL ERROR: bad read of chgres.inp.sfc header."
   print*,"- IRET IS ", iret
   call errexit(110)
 endif

 nsst_year=idate(1)
 nsst_mon=idate(2)
 nsst_day=idate(3)
 nsst_hour=idate(4)
 nsst_fhour=float(nfhour)

 levtyp='sfc'
 lev=1
 nframe=0

 allocate(dummy(imi*jmi))

!-----------------------------------------------------------------------
! Read land mask into its own variable
!-----------------------------------------------------------------------

 call nemsio_readrecv(gfile,"land",levtyp,lev, &
           dummy,nframe,iret)

 if (iret /= 0) then
   print*,"- FATAL ERROR: bad read of landmask record."
   print*,"- IRET IS ", iret
   call errexit(113)
 endif

 mask_input = reshape (dummy, (/imi,jmi/))

!-----------------------------------------------------------------------
! Read remaining records into nsst_input data structure.
! Note: fv3gfs files do not contain 'ifd' or 'zm' records.  Set
! to default values per recommendation of nsst developer.
!-----------------------------------------------------------------------

 print*,"- READ DATA RECORDS"

 do j = 1, nrec
   if (trim(recname(j)) == 'zm') then  
     nsst_input(:,:,j) = 0.0
     cycle
   endif
   if (trim(recname(j)) == 'ifd') then
     nsst_input(:,:,j) = 1.0
     cycle
   endif
   call nemsio_readrecv(gfile,recname(j),levtyp,lev, &
             dummy,nframe,iret)
   if (iret /= 0) then
     print*,"- FATAL ERROR: bad read of chgres.inp.sfc."
     print*,"- IRET IS ", iret
     call errexit(114)
   endif
   nsst_input(:,:,j) = reshape (dummy, (/imi,jmi/))
 enddo

 deallocate(dummy)

 call nemsio_close(gfile,iret=iret)

 END SUBROUTINE READ_FV3GFS_NSST_DATA_NEMSIO

 SUBROUTINE READ_GFS_NSST_DATA_NEMSIO (MASK_INPUT,NSST_INPUT,IMI,JMI, &
                 NUM_NSST_FIELDS,NSST_YEAR,NSST_MON,NSST_DAY,    & 
                 NSST_HOUR,NSST_FHOUR)

!-----------------------------------------------------------------------
! Subroutine: read nsst data from a gfs nemsio file
!
! Author: George Gayno/EMC
!
! Abstract: Reads an nsst file in nemsio format.  Places data
!           in the "nsst_input" array as expected by routine
!           nsst_chgres.
!
! Input files: 
!    "chgres.inp.nst" - input nsst nemsio file
!
! Output files:  none
!
! History:
!   2016-04-05   Gayno - Initial version
!
! Condition codes:  all non-zero codes are fatal
!   109 - bad open of nst file "chgres.inp.nst"
!   110 - bad read of "chgres.inp.nst" header
!   111 - the program assumes that the resolution of the
!         nst grid matches the input surface grid.  if
!         they are not the same, stop procoessing.
!   112 - the nst file does not have the 19 required records.
!   113 - bad read of landmask record.
!   114 - bad read of an nst file record.
!-----------------------------------------------------------------------

 use nemsio_module

 implicit none

 character(len=3)        :: levtyp
 character(len=8)        :: recname(19)

 integer, intent(in)     :: imi, jmi, num_nsst_fields
 integer, intent(out)    :: nsst_year, nsst_mon
 integer, intent(out)    :: nsst_day, nsst_hour

 real,    intent(out)    :: mask_input(imi,jmi)
 real,    intent(out)    :: nsst_input(imi,jmi,num_nsst_fields)
 real,    intent(out)    :: nsst_fhour

 integer(nemsio_intkind) :: iret, nrec, dimx, dimy, lev, nframe
 integer(nemsio_intkind) :: idate(7), nfhour

 integer                 :: j

 real(nemsio_realkind),allocatable :: dummy(:)

 type(nemsio_gfile)      :: gfile

 data recname   /"land    ", "xt      ", "xs      ", &
                 "xu      ", "xv      ", "xz      ", &
                 "zm      ", "xtts    ", "xzts    ", &
                 "dtcool  ", "zc      ", "c0      ", &
                 "cd      ", "w0      ", "wd      ", &
                 "dconv   ", "ifd     ", "tref    ", &
                 "qrain   " /

 print*,"- READ INPUT NSST DATA IN NEMSIO FORMAT"

 call nemsio_open(gfile, "chgres.inp.nst", "read", iret=iret)
 if (iret /= 0) then
   print*,"- FATAL ERROR: bad open of chgres.inp.nst."
   print*,"- IRET IS ", iret
   call errexit(109)
 endif

 print*,"- READ FILE HEADER"
 call nemsio_getfilehead(gfile,iret=iret,nrec=nrec,dimx=dimx, &
           dimy=dimy,idate=idate,nfhour=nfhour)
 if (iret /= 0) then
   print*,"- FATAL ERROR: bad read of chgres.inp.nst header."
   print*,"- IRET IS ", iret
   call errexit(110)
 endif

 if (dimx /= imi .or. dimy /= jmi) then
   print*,"- FATAL ERROR: nst and sfc file resolution"
   print*,"- must be the same."
   call errexit(111)
 endif

 if (nrec /= 19) then
   print*,"- FATAL ERROR: nst file has wrong number of records."
   call errexit(112)
 endif

 nsst_year=idate(1)
 nsst_mon=idate(2)
 nsst_day=idate(3)
 nsst_hour=idate(4)
 nsst_fhour=float(nfhour)

 levtyp='sfc'
 lev=1
 nframe=0

 allocate(dummy(imi*jmi))

!-----------------------------------------------------------------------
! Read land mask.  Note: older file versions use 'slmsk'
! as the header id.
!-----------------------------------------------------------------------

 call nemsio_readrecv(gfile,recname(1),levtyp,lev, &
           dummy,nframe,iret)
 if (iret /= 0) then
   call nemsio_readrecv(gfile,"slmsk",levtyp,lev, &
             dummy,nframe,iret)
   if (iret /= 0) then
     print*,"- FATAL ERROR: bad read of landmask record."
     print*,"- IRET IS ", iret
     call errexit(113)
   endif
 endif
 mask_input = reshape (dummy, (/imi,jmi/))

 print*,"- READ DATA RECORDS"
 do j = 2, nrec
   call nemsio_readrecv(gfile,recname(j),levtyp,lev, &
             dummy,nframe,iret)
   if (iret /= 0) then
     print*,"- FATAL ERROR: bad read of chgres.inp.nst."
     print*,"- IRET IS ", iret
     call errexit(114)
   endif
   nsst_input(:,:,j-1) = reshape (dummy, (/imi,jmi/))
 enddo

 deallocate(dummy)

 call nemsio_close(gfile,iret=iret)

 END SUBROUTINE READ_GFS_NSST_DATA_NEMSIO

 SUBROUTINE READ_GFS_SFC_HEADER_NEMSIO (IMI,JMI,IVSI,LSOILI, &
                 FCSTHOUR,IDATE4O,KGDS_INPUT)

 USE NEMSIO_MODULE

 IMPLICIT NONE

 INTEGER, INTENT(OUT)     :: IMI,JMI,IVSI,LSOILI,IDATE4O(4)
 INTEGER, INTENT(OUT)     :: KGDS_INPUT(200)

 REAL, INTENT(OUT)        :: FCSTHOUR

 CHARACTER(LEN=8)         :: FILETYPE

 INTEGER(NEMSIO_INTKIND)  :: DIMX, DIMY, IRET, VERSION
 INTEGER(NEMSIO_INTKIND)  :: NSOIL, IDATE(7), NFHOUR

 TYPE(NEMSIO_GFILE)       :: GFILEISFC

 CALL NEMSIO_OPEN(GFILEISFC,'chgres.inp.sfc','read',IRET=IRET)
 IF (IRET /= 0) THEN
   PRINT*,"FATAL ERROR OPENING chgres.inp.sfc"
   PRINT*,"IRET IS: ",IRET
   CALL ERREXIT(119)
 ENDIF

 CALL NEMSIO_GETFILEHEAD(GFILEISFC,GTYPE=FILETYPE,IRET=IRET, &
           VERSION=VERSION, DIMX=DIMX, DIMY=DIMY, NSOIL=NSOIL, &
           IDATE=IDATE, NFHOUR=NFHOUR)
 IF (IRET /= 0) THEN
   PRINT*,"FATAL ERROR READING chgres.inp.sfc FILE HEADER."
   PRINT*,"IRET IS: ",IRET
   CALL ERREXIT(120)
 ENDIF

! check bad status

 CALL NEMSIO_CLOSE(GFILEISFC,IRET=IRET)

 IMI        = DIMX
 JMI        = DIMY
 LSOILI     = NSOIL
 IVSI       = VERSION
 FCSTHOUR   = FLOAT(NFHOUR)
 IDATE4O(1) = IDATE(4)  ! HOUR
 IDATE4O(2) = IDATE(2)  ! MONTH
 IDATE4O(3) = IDATE(3)  ! DAY
 IDATE4O(4) = IDATE(1)  ! YEAR
     
 KGDS_INPUT = 0
 KGDS_INPUT(1) = 4          ! OCT 6 - TYPE OF GRID (GAUSSIAN)
 KGDS_INPUT(2) = IMI        ! OCT 7-8 - # PTS ON LATITUDE CIRCLE
 KGDS_INPUT(3) = JMI        ! OCT 9-10 - # PTS ON LONGITUDE CIRCLE
 KGDS_INPUT(4) = 90000      ! OCT 11-13 - LAT OF ORIGIN
 KGDS_INPUT(5) = 0          ! OCT 14-16 - LON OF ORIGIN
 KGDS_INPUT(6) = 128        ! OCT 17 - RESOLUTION FLAG
 KGDS_INPUT(7) = -90000     ! OCT 18-20 - LAT OF EXTREME POINT
 KGDS_INPUT(8) = NINT(-360000./IMI)  ! OCT 21-23 - LON OF EXTREME POINT
 KGDS_INPUT(9)  = NINT((360.0 / FLOAT(IMI))*1000.0)
                            ! OCT 24-25 - LONGITUDE DIRECTION INCR.
 KGDS_INPUT(10) = JMI /2    ! OCT 26-27 - NUMBER OF CIRCLES POLE TO EQUATOR
 KGDS_INPUT(12) = 255       ! OCT 29 - RESERVED
 KGDS_INPUT(20) = 255       ! OCT 5  - NOT USED, SET TO 255

 END SUBROUTINE READ_GFS_SFC_HEADER_NEMSIO

 SUBROUTINE READ_GFS_SFC_HEADER_SFCIO (NSFCI,IMI,JMI,IVSI,LSOILI, &
                 FCSTHOUR,IDATE4O,KGDS_INPUT)

 USE SFCIO_MODULE

 INTEGER, INTENT(IN)  :: NSFCI
 INTEGER, INTENT(OUT) :: IMI,JMI,IVSI,LSOILI,IDATE4O(4)
 INTEGER, INTENT(OUT) :: KGDS_INPUT(200)
 INTEGER              :: IRET

 REAL, INTENT(OUT)    :: FCSTHOUR

 TYPE(SFCIO_HEAD)     :: SFCHEADI

 CALL SFCIO_SROPEN(NSFCI,'chgres.inp.sfc',IRET)
 IF (IRET /= 0) THEN
   PRINT*,"FATAL ERROR OPENING chgres.inp.sfc"
   PRINT*,"IRET IS: ", IRET
   CALL ERREXIT(121)
 ENDIF

 CALL SFCIO_SRHEAD(NSFCI,SFCHEADI,IRET)
 IF (IRET /= 0) THEN
   PRINT*,"FATAL ERROR READING chgres.inp.sfc HEADER"
   PRINT*,"IRET IS: ", IRET
   CALL ERREXIT(122)
 ENDIF

 CALL SFCIO_SCLOSE(NSFCI,IRET)

 IMI = SFCHEADI%LONB
 JMI = SFCHEADI%LATB
 IVSI = SFCHEADI%IVS
 LSOILI = SFCHEADI%LSOIL
 FCSTHOUR = SFCHEADI%FHOUR
 IDATE4O = SFCHEADI%IDATE

 KGDS_INPUT = 0
 KGDS_INPUT(1) = 4          ! OCT 6 - TYPE OF GRID (GAUSSIAN)
 KGDS_INPUT(2) = IMI        ! OCT 7-8 - # PTS ON LATITUDE CIRCLE
 KGDS_INPUT(3) = JMI        ! OCT 9-10 - # PTS ON LONGITUDE CIRCLE
 KGDS_INPUT(4) = 90000      ! OCT 11-13 - LAT OF ORIGIN
 KGDS_INPUT(5) = 0          ! OCT 14-16 - LON OF ORIGIN
 KGDS_INPUT(6) = 128        ! OCT 17 - RESOLUTION FLAG
 KGDS_INPUT(7) = -90000     ! OCT 18-20 - LAT OF EXTREME POINT
 KGDS_INPUT(8) = NINT(-360000./IMI)  ! OCT 21-23 - LON OF EXTREME POINT
 KGDS_INPUT(9)  = NINT((360.0 / FLOAT(IMI))*1000.0)
                            ! OCT 24-25 - LONGITUDE DIRECTION INCR.
 KGDS_INPUT(10) = JMI /2    ! OCT 26-27 - NUMBER OF CIRCLES POLE TO EQUATOR
 KGDS_INPUT(12) = 255       ! OCT 29 - RESERVED
 KGDS_INPUT(20) = 255       ! OCT 5  - NOT USED, SET TO 255

 END SUBROUTINE READ_GFS_SFC_HEADER_SFCIO

 SUBROUTINE READ_FV3GFS_SFC_DATA_NEMSIO (IMI, JMI, LSOILI, SFCINPUT, &
                                 F10MI, T2MI, Q2MI,  &
                                 UUSTARI, FFMMI, FFHHI, SRFLAGI, &
                                 TPRCPI)

 USE NEMSIO_MODULE
 USE SURFACE_CHGRES

 INTEGER, INTENT(IN)  :: IMI, JMI, LSOILI

 REAL, INTENT(OUT)    :: F10MI(IMI,JMI), T2MI(IMI,JMI)
 REAL, INTENT(OUT)    :: Q2MI(IMI,JMI), UUSTARI(IMI,JMI)
 REAL, INTENT(OUT)    :: FFMMI(IMI,JMI), FFHHI(IMI,JMI)
 REAL, INTENT(OUT)    :: SRFLAGI(IMI,JMI), TPRCPI(IMI,JMI)

 TYPE(SFC2D)              :: SFCINPUT
 TYPE(NEMSIO_GFILE)       :: GFILEISFC

 INTEGER(NEMSIO_INTKIND)  :: IRET

 REAL(NEMSIO_REALKIND)    :: TMP(IMI*JMI)

 CALL NEMSIO_OPEN(GFILEISFC,'chgres.inp.sfc','read',IRET=IRET)
 IF(IRET /= 0)THEN
   PRINT*,"FATAL ERROR OPENING chgres.inp.sfc"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(244)
 ENDIF

 SRFLAGI = 0.0  ! NOT IN FILE.  SET TO ZERO.

 CALL NEMSIO_READRECV(GFILEISFC, 'ffhh', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 FFHHI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'ffmm', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 FFMMI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'f10m', '10 m above gnd', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 F10MI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', '2 m above gnd', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 T2MI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'spfh', '2 m above gnd', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 Q2MI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'fricv', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 UUSTARI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'tprcp', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 TPRCPI = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'alnsf', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%ALNSF = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'alnwf', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%ALNWF = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'alvsf', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%ALVSF = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'cnwat', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%CANOPY_MC = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'veg', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%GREENFRC = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'facsf', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%FACSF = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'facwf', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%FACWF = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SKIN_TEMP = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'land', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%LSMASK = RESHAPE(TMP, (/IMI,JMI/) )

 DO J = 1, JMI
 DO I = 1, IMI
   SFCINPUT%SEA_ICE_FLAG(I,J) = 0
   IF(NINT(SFCINPUT%LSMASK(I,J))==2) THEN
     SFCINPUT%LSMASK(I,J)=0.
     SFCINPUT%SEA_ICE_FLAG(I,J) = 1
   ENDIF
 ENDDO
 ENDDO

 CALL NEMSIO_READRECV(GFILEISFC, 'sfcr', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%Z0 = RESHAPE(TMP, (/IMI,JMI/) )
 SFCINPUT%Z0 = SFCINPUT%Z0 * 100.0  ! convert to cm
 
 CALL NEMSIO_READRECV(GFILEISFC, 'orog', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%OROG = RESHAPE(TMP, (/IMI,JMI/) )
 
 CALL NEMSIO_READRECV(GFILEISFC, 'vtype', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%VEG_TYPE = NINT(RESHAPE(TMP, (/IMI,JMI/) ))
 
 CALL NEMSIO_READRECV(GFILEISFC, 'sotyp', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOIL_TYPE = NINT(RESHAPE(TMP, (/IMI,JMI/) ))
 
 CALL NEMSIO_READRECV(GFILEISFC, 'weasd', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SNOW_LIQ_EQUIV = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'icec', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SEA_ICE_FRACT = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'icetk', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SEA_ICE_DEPTH = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'snoalb', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%MXSNOW_ALB = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'snod', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SNOW_DEPTH = RESHAPE(TMP, (/IMI,JMI/) )
 SFCINPUT%SNOW_DEPTH = SFCINPUT%SNOW_DEPTH * 1000.0 ! convert to mm

 CALL NEMSIO_READRECV(GFILEISFC, 'sltyp', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SLOPE_TYPE = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'shdmin', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%GREENFRC_MIN = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'shdmax', 'sfc', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%GREENFRC_MAX = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soilw', '0-10 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_TOT(:,:,1) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soilw', '10-40 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_TOT(:,:,2) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soilw', '40-100 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_TOT(:,:,3) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soilw', '100-200 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_TOT(:,:,4) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soill', '0-10 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_LIQ(:,:,1) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soill', '10-40 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_LIQ(:,:,2) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soill', '40-100 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_LIQ(:,:,3) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'soill', '100-200 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOILM_LIQ(:,:,4) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', '0-10 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOIL_TEMP(:,:,1) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', '10-40 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOIL_TEMP(:,:,2) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', '40-100 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOIL_TEMP(:,:,3) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_READRECV(GFILEISFC, 'tmp', '100-200 cm down', 1, TMP, IRET=IRET)
 IF (IRET /= 0) GOTO 99
 SFCINPUT%SOIL_TEMP(:,:,4) = RESHAPE(TMP, (/IMI,JMI/) )

 CALL NEMSIO_CLOSE(GFILEISFC, IRET=IRET)

 RETURN

 99 CONTINUE
 PRINT*,"FATAL ERROR READING DATA FROM chgres.inp.sfc"
 PRINT*,"IRET IS ", IRET
 CALL ERREXIT(245)

 END SUBROUTINE READ_FV3GFS_SFC_DATA_NEMSIO

 SUBROUTINE READ_GFS_SFC_DATA_NEMSIO (IMI, JMI, LSOILI, IVSI, SFCINPUT, &
                                 F10MI, T2MI, Q2MI,  &
                                 UUSTARI, FFMMI, FFHHI, SRFLAGI, &
                                 TPRCPI)

 USE NEMSIO_MODULE
 USE NEMSIO_GFS
 USE SURFACE_CHGRES

 IMPLICIT NONE

 INTEGER, INTENT(IN)  :: IMI, JMI, LSOILI, IVSI

 REAL, INTENT(OUT)    :: F10MI(IMI,JMI), T2MI(IMI,JMI)
 REAL, INTENT(OUT)    :: Q2MI(IMI,JMI), UUSTARI(IMI,JMI)
 REAL, INTENT(OUT)    :: FFMMI(IMI,JMI), FFHHI(IMI,JMI)
 REAL, INTENT(OUT)    :: SRFLAGI(IMI,JMI), TPRCPI(IMI,JMI)

 INTEGER(NEMSIO_INTKIND)  :: IRET
 INTEGER                  :: I, J, L

 TYPE(SFC2D)          :: SFCINPUT
 TYPE(NEMSIO_GFILE)   :: GFILEISFC
 TYPE(NEMSIO_DBTA)    :: GFSDATAI

 CALL NEMSIO_OPEN(GFILEISFC,'chgres.inp.sfc','read',IRET=IRET)
 IF(IRET /= 0)THEN
   PRINT*,"FATAL ERROR OPENING chgres.inp.sfc"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(144)
 ENDIF

 CALL NEMSIO_GFS_ALSFC(IMI, JMI, LSOILI, GFSDATAI)

 CALL NEMSIO_GFS_RDSFC(GFILEISFC,GFSDATAI,IRET)
 IF(IRET /= 0)THEN
   PRINT*,"FATAL ERROR READING DATA FROM chgres.inp.sfc"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(145)
 ENDIF

 CALL NEMSIO_CLOSE(GFILEISFC, IRET=IRET)

!$OMP PARALLEL DO PRIVATE(I,J)
 DO J = 1, JMI
   DO I = 1, IMI

     SFCINPUT%ALNSF(I,J) = GFSDATAI%ALNSF(I,J)
     SFCINPUT%ALNWF(I,J) = GFSDATAI%ALNWF(I,J)
     SFCINPUT%ALVSF(I,J) = GFSDATAI%ALVSF(I,J)
     SFCINPUT%ALVWF(I,J) = GFSDATAI%ALVWF(I,J)
     SFCINPUT%CANOPY_MC(I,J) = GFSDATAI%CANOPY(I,J)
     SFCINPUT%GREENFRC(I,J) = GFSDATAI%VFRAC(I,J)
     SFCINPUT%FACSF(I,J) = GFSDATAI%FACSF(I,J)
     SFCINPUT%FACWF(I,J) = GFSDATAI%FACWF(I,J)
     SFCINPUT%SKIN_TEMP(I,J) = GFSDATAI%TSEA(I,J)
     SFCINPUT%LSMASK(I,J) = GFSDATAI%SLMSK(I,J)
     SFCINPUT%SEA_ICE_FLAG(I,J) = 0
     IF(NINT(SFCINPUT%LSMASK(I,J))==2) THEN
       SFCINPUT%LSMASK(I,J)=0.
       SFCINPUT%SEA_ICE_FLAG(I,J) = 1
      ENDIF
     SFCINPUT%Z0(I,J) = GFSDATAI%ZORL(I,J)
     SFCINPUT%OROG(I,J)         = GFSDATAI%OROG(I,J)
     SFCINPUT%VEG_TYPE(I,J)     = NINT(GFSDATAI%VTYPE(I,J))
     SFCINPUT%SOIL_TYPE(I,J)    = NINT(GFSDATAI%STYPE(I,J))
     SFCINPUT%SNOW_LIQ_EQUIV(I,J) = GFSDATAI%SHELEG(I,J)

   ENDDO
 ENDDO
!$OMP END PARALLEL DO
 
 DO L = 1, LSOILI
!$OMP PARALLEL DO PRIVATE(I,J)
   DO J = 1, JMI
     DO I = 1, IMI
       SFCINPUT%SOILM_TOT(I,J,L) = GFSDATAI%SMC(I,J,L)
       SFCINPUT%SOIL_TEMP(I,J,L) = GFSDATAI%STC(I,J,L)
     ENDDO
   ENDDO
!$OMP END PARALLEL DO
 ENDDO

 SRFLAGI = 0.0
 TPRCPI  = 0.0

 IF (IVSI >= 200501) THEN
!$OMP PARALLEL DO PRIVATE(I,J)
   DO J = 1, JMI
     DO I = 1, IMI
       SFCINPUT%SEA_ICE_FRACT(I,J) = GFSDATAI%FICE(I,J)
       SFCINPUT%SEA_ICE_DEPTH(I,J) = GFSDATAI%HICE(I,J)
       SFCINPUT%MXSNOW_ALB(I,J)    = GFSDATAI%SNOALB(I,J)
       SFCINPUT%SNOW_DEPTH(I,J)    = GFSDATAI%SNWDPH(I,J)
       SFCINPUT%SLOPE_TYPE(I,J)    = NINT(GFSDATAI%SLOPE(I,J))
       SFCINPUT%GREENFRC_MAX(I,J)  = GFSDATAI%SHDMAX(I,J)
       SFCINPUT%GREENFRC_MIN(I,J)  = GFSDATAI%SHDMIN(I,J)
       SRFLAGI(I,J)                = GFSDATAI%SRFLAG(I,J)
       TPRCPI(I,J)                 = GFSDATAI%TPRCP(I,J)
     ENDDO
   ENDDO
!$OMP END PARALLEL DO

   DO L=1,LSOILI
!$OMP PARALLEL DO PRIVATE(I,J)
     DO J = 1, JMI
       DO I = 1, IMI
         SFCINPUT%SOILM_LIQ(I,J,L) = GFSDATAI%SLC(I,J,L)
       ENDDO
     ENDDO
   ENDDO

 END IF  ! IVS >= 200501

!$OMP PARALLEL DO PRIVATE(I,J)
 DO J = 1, JMI
   DO I = 1, IMI
     F10MI(I,J) = GFSDATAI%F10M(I,J)
     T2MI(I,J) = GFSDATAI%T2M(I,J)
     Q2MI(I,J) = GFSDATAI%Q2M(I,J)
     UUSTARI(I,J) = GFSDATAI%UUSTAR(I,J)
     FFMMI(I,J) = GFSDATAI%FFMM(I,J)
     FFHHI(I,J) = GFSDATAI%FFHH(I,J)
   ENDDO
 ENDDO
!$OMP END PARALLEL DO

 END SUBROUTINE READ_GFS_SFC_DATA_NEMSIO

 SUBROUTINE READ_GFS_SFC_DATA_SFCIO (NSFCI, IMI, JMI, SFCINPUT,  &
                              F10MI, T2MI, Q2MI, &
                              UUSTARI, FFMMI, FFHHI, SRFLAGI, &
                              TPRCPI)

 USE SFCIO_MODULE
 USE SURFACE_CHGRES

 IMPLICIT NONE

 INTEGER, INTENT(IN)  :: NSFCI, IMI, JMI
 INTEGER              :: I,J,L, IRET

 REAL, INTENT(OUT)    :: F10MI(IMI,JMI), T2MI(IMI,JMI)
 REAL, INTENT(OUT)    :: Q2MI(IMI,JMI), UUSTARI(IMI,JMI)
 REAL, INTENT(OUT)    :: FFMMI(IMI,JMI), FFHHI(IMI,JMI)
 REAL, INTENT(OUT)    :: SRFLAGI(IMI,JMI), TPRCPI(IMI,JMI)

 TYPE(SFC2D)          :: SFCINPUT
 TYPE(SFCIO_HEAD)     :: SFCHEADI
 TYPE(SFCIO_DBTA)     :: SFCDATAI

 CALL SFCIO_SROPEN(NSFCI,'chgres.inp.sfc',IRET)
 IF(IRET /=0) THEN
   PRINT*,"FATAL ERROR OPENING chgres.inp.sfc"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(155)
 ENDIF

 CALL SFCIO_SRHEAD(NSFCI,SFCHEADI,IRET)
 IF(IRET /=0) THEN
   PRINT*,"FATAL ERROR READING chgres.inp.sfc HEADER"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(156)
 ENDIF

 CALL SFCIO_ALDBTA(SFCHEADI,SFCDATAI,IRET)
 IF(IRET.NE.0) THEN
   PRINT*,"FATAL ERROR ALLOCATING SFC DATA STRUCTURE"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(158)
 ENDIF

 CALL SFCIO_SRDBTA(NSFCI,SFCHEADI,SFCDATAI,IRET)
 IF(IRET /=0) THEN
   PRINT*,"FATAL ERROR READING chgres.inp.sfc DATA"
   PRINT*,"IRET IS ", IRET
   CALL ERREXIT(157)
 ENDIF

 CALL SFCIO_SCLOSE(NSFCI,IRET)

!$OMP PARALLEL DO PRIVATE(I,J)

 DO J = 1, SFCHEADI%LATB
 DO I = 1, SFCHEADI%LONB

   SFCINPUT%ALNSF(I,J) = SFCDATAI%ALNSF(I,J)
   SFCINPUT%ALNWF(I,J) = SFCDATAI%ALNWF(I,J)
   SFCINPUT%ALVSF(I,J) = SFCDATAI%ALVSF(I,J)
   SFCINPUT%ALVWF(I,J) = SFCDATAI%ALVWF(I,J)
   SFCINPUT%CANOPY_MC(I,J) = SFCDATAI%CANOPY(I,J)
   SFCINPUT%GREENFRC(I,J) = SFCDATAI%VFRAC(I,J)
   SFCINPUT%FACSF(I,J) = SFCDATAI%FACSF(I,J)
   SFCINPUT%FACWF(I,J) = SFCDATAI%FACWF(I,J)
   SFCINPUT%SKIN_TEMP(I,J) = SFCDATAI%TSEA(I,J)
   SFCINPUT%LSMASK(I,J) = SFCDATAI%SLMSK(I,J)
   SFCINPUT%SEA_ICE_FLAG(I,J) = 0
   IF(NINT(SFCINPUT%LSMASK(I,J))==2) THEN
     SFCINPUT%LSMASK(I,J)=0.
     SFCINPUT%SEA_ICE_FLAG(I,J) = 1
   ENDIF
   SFCINPUT%Z0(I,J) = SFCDATAI%ZORL(I,J)
   SFCINPUT%OROG(I,J)         = SFCDATAI%OROG(I,J)
   SFCINPUT%VEG_TYPE(I,J)     = NINT(SFCDATAI%VTYPE(I,J))
   SFCINPUT%SOIL_TYPE(I,J)    = NINT(SFCDATAI%STYPE(I,J))
   SFCINPUT%SNOW_LIQ_EQUIV(I,J) = SFCDATAI%SHELEG(I,J)

 ENDDO
 ENDDO

!$OMP END PARALLEL DO

 DO L = 1, SFCHEADI%LSOIL
!$OMP PARALLEL DO PRIVATE(I,J)
   DO J = 1, SFCHEADI%LATB
     DO I = 1, SFCHEADI%LONB
       SFCINPUT%SOILM_TOT(I,J,L) = SFCDATAI%SMC(I,J,L)
       SFCINPUT%SOIL_TEMP(I,J,L) = SFCDATAI%STC(I,J,L)
     ENDDO
   ENDDO
!$OMP END PARALLEL DO
 ENDDO

 SRFLAGI = 0.0
 TPRCPI  = 0.0

 IF (SFCHEADI%IVS >= 200501) THEN
!$OMP PARALLEL DO PRIVATE(I,J)
   DO J = 1, SFCHEADI%LATB
     DO I = 1, SFCHEADI%LONB
       SFCINPUT%SEA_ICE_FRACT(I,J) = SFCDATAI%FICE(I,J)
       SFCINPUT%SEA_ICE_DEPTH(I,J) = SFCDATAI%HICE(I,J)
       SFCINPUT%MXSNOW_ALB(I,J)    = SFCDATAI%SNOALB(I,J)
       SFCINPUT%SNOW_DEPTH(I,J)    = SFCDATAI%SNWDPH(I,J)
       SFCINPUT%SLOPE_TYPE(I,J)    = NINT(SFCDATAI%SLOPE(I,J))
       SFCINPUT%GREENFRC_MAX(I,J)  = SFCDATAI%SHDMAX(I,J)
       SFCINPUT%GREENFRC_MIN(I,J)  = SFCDATAI%SHDMIN(I,J)
       SRFLAGI(I,J)                = SFCDATAI%SRFLAG(I,J)
       TPRCPI(I,J)                 = SFCDATAI%TPRCP(I,J)
     ENDDO
   ENDDO
!$OMP END PARALLEL DO

   DO L=1,SFCHEADI%LSOIL
!$OMP PARALLEL DO PRIVATE(I,J)
     DO J = 1, SFCHEADI%LATB
       DO I = 1, SFCHEADI%LONB
         SFCINPUT%SOILM_LIQ(I,J,L) = SFCDATAI%SLC(I,J,L)
       ENDDO
     ENDDO
   ENDDO

 END IF  ! IVS >= 200501

!$OMP PARALLEL DO PRIVATE(I,J)
 DO J = 1, SFCHEADI%LATB
   DO I = 1, SFCHEADI%LONB
     F10MI(I,J) = SFCDATAI%F10M(I,J)
     T2MI(I,J) = SFCDATAI%T2M(I,J)
     Q2MI(I,J) = SFCDATAI%Q2M(I,J)
     UUSTARI(I,J) = SFCDATAI%UUSTAR(I,J)
     FFMMI(I,J) = SFCDATAI%FFMM(I,J)
     FFHHI(I,J) = SFCDATAI%FFHH(I,J)
   ENDDO
 ENDDO
!$OMP END PARALLEL DO

 CALL SFCIO_AXDBTA(SFCDATAI,IRET)

 END SUBROUTINE READ_GFS_SFC_DATA_SFCIO
