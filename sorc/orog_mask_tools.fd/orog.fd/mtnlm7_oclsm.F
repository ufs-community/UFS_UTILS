C> @file
C> Terrain maker for global spectral model.
C> @author Mark Iredell @date 92-04-16
  
C> This program creates 7 terrain-related files computed from the
C> GMTED2010 terrain dataset. The model physics grid parameters and
C> spectral truncation and filter parameters are read by this program as
C> input.
C>      
C> The 7 files produced are:
C> 1. sea-land mask on model physics grid
C> 2. gridded orography on model physics grid
C> 3. mountain std dev on model physics grid
C> 4. spectral orography in spectral domain
C> 5. unfiltered gridded orography on model physics grid
C> 6. grib sea-land mask on model physics grid
C> 7. grib gridded orography on model physics grid
C>
C> PROGRAM HISTORY LOG:
C> -  92-04-16  IREDELL
C> -  98-02-02  IREDELL  FILTER
C> -  98-05-31  HONG Modified for subgrid orography used in Kim's scheme
C> -  98-12-31  HONG Modified for high-resolution GTOPO orography
C> -  99-05-31  HONG Modified for getting OL4 (mountain fraction)
C>  - 00-02-10  Moorthi's modifications
C> -  00-04-11  HONG Modified for reduced grids
C> -  00-04-12  Iredell Modified for reduced grids
C> -  02-01-07  (*j*) modified for principal axes of orography
C>             There are now 14 files, 4 additional for lm mb
C>  - 04-04-04  (*j*) re-Test on IST/ilen calc for sea-land mask(*j*)
C>  - 04-09-04   minus sign here in MAKEOA IST and IEN as in MAKEMT!
C>  - 05-09-05   if test on HK and HLPRIM for GAMMA SQRT
C>  - 07-08-07   replace 8' with 30" incl GICE, conintue w/ S-Y. lake slm
C>  - 08-08-07  All input 30", UMD option, and filter as described below
C>       
C>   INPUT FILES:
C>  -   UNIT5      - PHYSICS LONGITUDES (IM), PHYSICS LATITUDES (JM),
C>                  SPECTRAL TRUNCATION (NM), RHOMBOIDAL FLAG (NR),
C>                  AND FIRST AND SECOND FILTER PARAMETERS (NF0,NF1).
C>                  RESPECTIVELY READ IN FREE FORMAT.
C>  -   NCID       - GMTED2010 USGS orography (NetCDF)
C>  -   NCID       - 30" UMD land cover mask. (NetCDF)
C>  -   NCID       - GICE Grumbine 30" RAMP Antarctica orog IMNx3601. (NetCDF)
C>
C>   OUTPUT FILES:
C>  -   UNIT51     - SEA-LAND MASK (IM,JM)
C>  -   UNIT52     - GRIDDED OROGRAPHY (IM,JM)
C>  -   UNIT54     - SPECTRAL OROGRAPHY ((NM+1)*((NR+1)*NM+2))
C>  -   UNIT55     - UNFILTERED GRIDDED OROGRAPHY (IM,JM)
C>  -   UNIT57     - GRIB GRIDDED OROGRAPHY (IM,JM)
C>
C>   SUBPROGRAMS CALLED:
C>  -   UNIQUE:
C>  -   TERSUB     - MAIN SUBPROGRAM
C>  -   SPLAT      - COMPUTE GAUSSIAN LATITUDES OR EQUALLY-SPACED LATITUDES
C>  -   LIBRARY:
C>  -   SPTEZ      - SPHERICAL TRANSFORM
C>
C> @return 0 for success, error code otherwise.
      implicit none
      integer imn, jmn
      character(len=256) :: MDL_GRID_FILE = "none"
      character(len=256) :: merge_file = "none"
      logical :: mask_only = .false.
      integer :: IM,JM,EFAC
      READ(5,*) MDL_GRID_FILE
      READ(5,*) mask_only
      READ(5,*) merge_file
      EFAC=0
      print*, "MASK_ONLY", mask_only
      print*, "MERGE_FILE ", trim(merge_file)
      print*, EFAC
      IMN = 360*120
      JMN = 180*120
      print *, ' Starting terr12 mtnlm7_slm30.f  IMN,JMN:',IMN,JMN

      call read_mdl_dims(mdl_grid_file, im, jm)
         
      CALL TERSUB(IMN,JMN,IM,JM,EFAC,
     &            MDL_GRID_FILE,MASK_ONLY,MERGE_FILE)
      STOP
      END

!> Driver routine to compute terrain.
!!
!! @param[in] IMN "i" dimension of the input terrain dataset.
!! @param[in] JMN "j" dimension of the input terrain dataset.
!! @param[in] IM "i" dimension of the model grid tile.
!! @param[in] JM "j" dimension of the model grid tile.
!! @param[in] EFAC Factor to adjust orography by its variance.
!! @param[in] OUTGRID The 'grid' file for the model tile.
!! grid. When specified, will be interpolated to model tile.
!! When not specified, program will create fields from
!! raw high-resolution topography data.
!! @param[in] MASK_ONLY Flag to generate the Land Mask only
!! @param[in] MERGE_FILE Ocean merge file 
!! @author Jordan Alpert NOAA/EMC
      SUBROUTINE TERSUB(IMN,JMN,IM,JM,EFAC,
     &     OUTGRID,MASK_ONLY,MERGE_FILE)
      implicit none
      include 'netcdf.inc'
C
      integer                      :: IMN,JMN,IM,JM
      character(len=*), intent(in) :: OUTGRID
      character(len=*), intent(in) :: MERGE_FILE

      logical, intent(in) :: mask_only

      real, parameter :: MISSING_VALUE=-9999.
      real, PARAMETER :: PI=3.1415926535897931
      integer, PARAMETER :: NMT=14

      integer :: efac,zsave1,zsave2
      integer :: i,j,nx,ny,ncid,js,jn,iw,ie,k,it,jt,error,id_dim
      integer :: id_var,fsize,wgta,IN,INW,INE,IS,ISW,ISE
      integer :: itest,jtest
      integer :: i_south_pole,j_south_pole,i_north_pole,j_north_pole
      integer(1) :: i3save
      integer(2) :: i2save

      integer, allocatable :: JST(:),JEN(:)

      integer, allocatable :: IST(:,:),IEN(:,:),ZSLMX(:,:)
      integer, allocatable :: ZAVG(:,:),ZSLM(:,:)
      integer(1), allocatable :: UMD(:,:)
      integer(2), allocatable :: glob(:,:)

      integer, allocatable :: IWORK(:,:,:)

      real :: DEGRAD,timef,tbeg,tend,tbeg1
      real :: DELXN,slma,oroa,vara,var4a,xn,XS

      real, allocatable :: XLAT(:),XLON(:),oaa(:),ola(:),GLAT(:)

      real, allocatable :: GEOLON(:,:),GEOLON_C(:,:),DX(:,:)
      real, allocatable :: GEOLAT(:,:),GEOLAT_C(:,:),DY(:,:)
      real, allocatable :: SLM(:,:),ORO(:,:),VAR(:,:)
      real, allocatable :: land_frac(:,:),lake_frac(:,:)
      real, allocatable :: THETA(:,:),GAMMA(:,:),SIGMA(:,:),ELVMAX(:,:)
      real, allocatable :: VAR4(:,:)
      real, allocatable :: WORK1(:,:),WORK2(:,:),WORK3(:,:),WORK4(:,:)
      real, allocatable :: WORK5(:,:),WORK6(:,:)
      real, allocatable :: tmpvar(:,:)
      real(4), allocatable::  GICE(:,:)

      real, allocatable :: OA(:,:,:),OL(:,:,:),HPRIME(:,:,:)

      logical :: fexist,opened
      logical :: is_south_pole(IM,JM), is_north_pole(IM,JM)

      tbeg1=timef()
      tbeg=timef()
      fsize = 65536
! integers
      allocate (JST(JM),JEN(JM))
      allocate (IST(IM,jm),IEN(IM,jm),ZSLMX(2700,1350))
      allocate (glob(IMN,JMN))

! reals
      allocate (oaa(4),ola(4),GLAT(JMN))

      allocate (ZAVG(IMN,JMN))
      allocate (ZSLM(IMN,JMN))
      allocate (UMD(IMN,JMN))

!
!  SET CONSTANTS AND ZERO FIELDS
!
      DEGRAD = 180./PI

      print *,' In TERSUB'

      print *,' Attempt to open/read UMD 30sec slmsk.'

      error=NF__OPEN("./landcover.umd.30s.nc",NF_NOWRITE,fsize,ncid)
      call netcdf_err(error, 'Open file landcover.umd.30s.nc' )
      error=nf_inq_varid(ncid, 'land_mask', id_var)
      call netcdf_err(error, 'Inquire varid of land_mask')
      error=nf_get_var_int1(ncid, id_var, UMD)
      call netcdf_err(error, 'Inquire data of land_mask')
      error = nf_close(ncid)

      print *,' UMD lake, UMD(50,50)=',UMD(50,50)
C
C- READ_G for global 30" terrain 
C
       print *,' Call read_g to read global topography'
       call read_g(glob)
! --- transpose even though glob 30" is from S to N and NCEP std is N to S
       do j=1,jmn/2
       do I=1,imn
        jt=jmn - j + 1
        i2save = glob(I,j)
        glob(I,j)=glob(I,jt)
        glob(I,jt) = i2save
       enddo
       enddo 
! --- transpose glob as USGS 30" is from dateline and NCEP std is 0
       do j=1,jmn
       do I=1,imn/2
        it=imn/2 + i 
        i2save = glob(i,J)
        glob(i,J)=glob(it,J)
        glob(it,J) = i2save
       enddo
       enddo 
       print *,' After read_g, glob(500,500)=',glob(500,500)
!

!  --- IMN,JMN
      print*, ' IM, JM, EFAC'
      print*, IM,JM,EFAC
       print *,'  imn,jmn,glob(imn,jmn)=',imn,jmn,glob(imn,jmn)
       print *,' UBOUND ZAVG=',UBOUND(ZAVG)
       print *,' UBOUND glob=',UBOUND(glob)
       print *,' UBOUND ZSLM=',UBOUND(ZSLM)
       print *,' UBOUND GICE=',IMN+1,3601
!
! ---  0 is ocean and 1 is land for slm
!
C
! --- ZSLM initialize with all land 1, ocean 0
      ZSLM=1
! --- ZAVG initialize from glob
      ZAVG=glob

! --- transpose mask even though glob 30" is from N to S and NCEP std is S to N
      do j=1,jmn/2
      do I=1,imn
        jt=jmn - j + 1
        i3save = UMD(I,j)
        UMD(I,j)=UMD(I,jt)
        UMD(I,jt) = i3save
      enddo
      enddo 
! --- transpose UMD as USGS 30" is from dateline and NCEP std is 0
      do j=1,jmn
      do i=1,imn/2
        it=imn/2 + i 
        i3save = UMD(i,J)
        UMD(i,J)=UMD(it,J)
        UMD(it,J) = i3save
      enddo
      enddo
! ---  Non-land is 0.
      do j=1,jmn
      do i=1,imn
         if ( UMD(i,j) .eq. 0 ) ZSLM(i,j) = 0
      enddo
      enddo

      deallocate (ZSLMX,UMD,glob)
! ---
! ---  Fixing an error in the topo 30" data set at pole (-9999).  
            do i=1,imn
            ZSLM(i,1)=0
            ZSLM(i,JMN)=1
            enddo
!
      allocate (GICE(IMN+1,3601))
!
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
      ENDDO
      print *,
     & ' Before GICE ZAVG(1,2)=',ZAVG(1,2),ZSLM(1,2)
      print *,
     & ' Before GICE ZAVG(1,12)=',ZAVG(1,12),ZSLM(1,12)
      print *,
     & ' Before GICE ZAVG(1,52)=',ZAVG(1,52),ZSLM(1,52)
      print *,
     & ' Before GICE ZAVG(1,112)=',ZAVG(1,JMN-112),ZSLM(1,112)

! Read 30-sec Antarctica RAMP data. Points scan from South
! to North, and from Greenwich to Greenwich.

      error=NF__OPEN("./topography.antarctica.ramp.30s.nc",
     &               NF_NOWRITE,fsize,ncid)
      call netcdf_err(error, 'Opening RAMP topo file' )
      error=nf_inq_varid(ncid, 'topo', id_var)
      call netcdf_err(error, 'Inquire varid of RAMP topo')
      error=nf_get_var_real(ncid, id_var, GICE)
      call netcdf_err(error, 'Inquire data of RAMP topo')
      error = nf_close(ncid)
      
      print *,' GICE 30" Antarctica RAMP orog 43201x3601 read OK'
      print *,' Processing! '
      print *,' Processing! '
      print *,' Processing! '
      do j = 1, 3601 
      do i = 1, IMN
         zsave1 = ZAVG(i,j)
         zsave2 = ZSLM(i,j)
       if( GICE(i,j) .ne. -99. .and.  GICE(i,j) .ne. -1.0 ) then
         if ( GICE(i,j) .gt. 0.) then 
              ZAVG(i,j) = int( GICE(i,j) + 0.5 )
!! --- for GICE values less than or equal to 0 (0, -1, or -99) then
!! --- radar-sat (RAMP) values are not valid and revert back to old orog 
                ZSLM(i,j) = 1
         endif
       endif
  152 format(1x,' ZAVG(i=',i4,' j=',i4,')=',i5,i3,
     &' orig:',i5,i4,' Lat=',f7.3,f8.2,'E',' GICE=',f8.1)
      enddo
      enddo

      deallocate (GICE)
C
C     COMPUTE MOUNTAIN DATA : ORO SLM VAR (Std Dev) OC
C
      allocate (GEOLON(IM,JM),GEOLON_C(IM+1,JM+1),DX(IM,JM))
      allocate (GEOLAT(IM,JM),GEOLAT_C(IM+1,JM+1),DY(IM,JM))
      allocate (SLM(IM,JM),ORO(IM,JM),VAR(IM,JM),VAR4(IM,JM))
      allocate (land_frac(IM,JM),lake_frac(IM,JM))

!--- reading grid file.

      inquire(file=trim(OUTGRID), exist=fexist)
      if(.not. fexist) then
         print*, "FATAL ERROR: file "//trim(OUTGRID)
         print*, "does not exist."
         CALL ERREXIT(4)
      endif
      do ncid = 103, 512
        inquire( ncid,OPENED=opened )
        if( .NOT.opened )exit
      end do

      print*, "outgrid=", trim(outgrid)
      error=NF__OPEN(trim(OUTGRID),NF_NOWRITE,fsize,ncid)
      call netcdf_err(error, 'Open file '//trim(OUTGRID) )
      error=nf_inq_dimid(ncid, 'nx', id_dim)
      call netcdf_err(error, 'inquire dimension nx from file '//
     &                   trim(OUTGRID) )
      nx = 2*IM
      ny = 2*JM
      print*, "Read the grid from file "//trim(OUTGRID)

      allocate(tmpvar(nx+1,ny+1))

      error=nf_inq_varid(ncid, 'x', id_var)
      call netcdf_err(error, 'inquire varid of x from file '
     &                   //trim(OUTGRID) )
      error=nf_get_var_double(ncid, id_var, tmpvar)
      call netcdf_err(error, 'inquire data of x from file '
     &                   //trim(OUTGRID) )
      !--- adjust lontitude to be between 0 and 360.
      do j = 1,ny+1; do i = 1,nx+1
         if(tmpvar(i,j) .NE. MISSING_VALUE) then
            if(tmpvar(i,j) .GT. 360) tmpvar(i,j) = tmpvar(i,j) - 360
            if(tmpvar(i,j) .LT. 0) tmpvar(i,j) = tmpvar(i,j) + 360
         endif
      enddo; enddo

      geolon(1:IM,1:JM) = tmpvar(2:nx:2,2:ny:2)
      geolon_c(1:IM+1,1:JM+1) = tmpvar(1:nx+1:2,1:ny+1:2)
         
      error=nf_inq_varid(ncid, 'y', id_var)
      call netcdf_err(error, 'inquire varid of y from file '
     &                   //trim(OUTGRID) )
      error=nf_get_var_double(ncid, id_var, tmpvar)
      call netcdf_err(error, 'inquire data of y from file '
     &                   //trim(OUTGRID) )
      geolat(1:IM,1:JM) = tmpvar(2:nx:2,2:ny:2)
      geolat_c(1:IM+1,1:JM+1) = tmpvar(1:nx+1:2,1:ny+1:2)

      call find_poles(tmpvar, nx, ny, i_north_pole, j_north_pole,
     &                     i_south_pole, j_south_pole)

      deallocate(tmpvar)
         
      call find_nearest_pole_points(i_north_pole, j_north_pole,
     &  i_south_pole, j_south_pole, im, jm, is_north_pole,
     &  is_south_pole)

      allocate(tmpvar(nx,ny))
      error=nf_inq_varid(ncid, 'area', id_var)
      call netcdf_err(error, 'inquire varid of area from file '
     &                   //trim(OUTGRID) )
      error=nf_get_var_double(ncid, id_var, tmpvar)
      call netcdf_err(error, 'inquire data of area from file '
     &                   //trim(OUTGRID) )

      do j = 1, jm
         do i = 1, im
            dx(i,j) = sqrt(tmpvar(2*i-1,2*j-1)+tmpvar(2*i,2*j-1)
     &                     +tmpvar(2*i-1,2*j  )+tmpvar(2*i,2*j  ))
            dy(i,j) = dx(i,j)
         enddo
      enddo
      deallocate(tmpvar)         

      tend=timef()
      write(6,*)' Timer 1 time= ',tend-tbeg
                                !
       tbeg=timef()

         IF (MERGE_FILE == 'none') then
           CALL MAKE_MASK(ZSLM,SLM,land_frac,GLAT,
     &      IM,JM,IMN,JMN,geolon_c,geolat_c)
            lake_frac=9999.9
         ELSE
           print*,'Read in external mask ',merge_file
           CALL READ_MASK(MERGE_FILE,SLM,land_frac,lake_frac,im,jm)
         ENDIF

         IF (MASK_ONLY) THEN
           print*,'Computing mask only.'
           CALL WRITE_MASK_NETCDF(IM,JM,SLM,land_frac,
     1                  1,1,GEOLON,GEOLAT)

           print*,' DONE.'
           STOP
         END IF

         CALL MAKEMT2(ZAVG,ZSLM,ORO,SLM,VAR,VAR4,GLAT,
     & IM,JM,IMN,JMN,geolon_c,geolat_c,lake_frac,land_frac)

      tend=timef()
      write(6,*)' MAKEMT2 time= ',tend-tbeg

       call minmxj(IM,JM,ORO,'     ORO')
       call minmxj(IM,JM,SLM,'     SLM')
       call minmxj(IM,JM,VAR,'     VAR')
       call minmxj(IM,JM,VAR4,'    VAR4')
C
C ===  Compute mtn principal coord HTENSR: THETA,GAMMA,SIGMA
C
      allocate (THETA(IM,JM),GAMMA(IM,JM),SIGMA(IM,JM),ELVMAX(IM,JM))
      tbeg=timef()
      CALL MAKEPC2(ZAVG,ZSLM,THETA,GAMMA,SIGMA,GLAT,
     1            IM,JM,IMN,JMN,geolon_c,geolat_c,SLM)
      tend=timef()
      write(6,*)' MAKEPC2 time= ',tend-tbeg

       call minmxj(IM,JM,THETA,'   THETA')
       call minmxj(IM,JM,GAMMA,'   GAMMA')
       call minmxj(IM,JM,SIGMA,'   SIGMA')
C
C     COMPUTE MOUNTAIN DATA : OA OL
C
       allocate (IWORK(IM,JM,4))
       allocate (OA(IM,JM,4),OL(IM,JM,4),HPRIME(IM,JM,14))
       allocate (WORK1(IM,JM),WORK2(IM,JM),WORK3(IM,JM),WORK4(IM,JM))
       allocate (WORK5(IM,JM),WORK6(IM,JM))

       call minmxj(IM,JM,ORO,'     ORO')

       print*, "calling MAKEOA2 to compute OA, OL"
       tbeg=timef()
       CALL MAKEOA2(ZAVG,zslm,VAR,GLAT,OA,OL,IWORK,ELVMAX,ORO,
     1            WORK1,WORK2,WORK3,WORK4,WORK5,WORK6,
     2            IM,JM,IMN,JMN,geolon_c,geolat_c,
     3            geolon,geolat,dx,dy,is_south_pole,is_north_pole)
       tend=timef()
       write(6,*)' MAKEOA2 time= ',tend-tbeg

!      Deallocate 2d vars
       deallocate(IST,IEN)
       deallocate (ZSLM,ZAVG)
       deallocate (dx,dy)
       deallocate (WORK2,WORK3,WORK4,WORK5,WORK6)

!      Deallocate 3d vars
       deallocate(IWORK)

       tbeg=timef()
       call minmxj(IM,JM,OA,'      OA')
       call minmxj(IM,JM,OL,'      OL')
       call minmxj(IM,JM,ELVMAX,'  ELVMAX')
       call minmxj(IM,JM,ORO,'     ORO')

       print *,' ===> Replacing ELVMAX with ELVMAX-ORO <=== ' 
       print *,' ===> if ELVMAX<=ORO replace with proxy <=== ' 
       print *,' ===> the sum of mean orog (ORO) and std dev <=== ' 
      DO J = 1,JM
      DO I = 1,IM
        if (ELVMAX(I,J) .lt. ORO(I,J) ) then
C---  subtracting off ORO leaves std dev (this should never happen)
       ELVMAX(I,J) = MAX(  3. * VAR(I,J),0.)
        else
       ELVMAX(I,J) = MAX( ELVMAX(I,J) - ORO(I,J),0.)
        endif
      ENDDO
      ENDDO
c
       call mnmxja(IM,JM,ELVMAX,itest,jtest,'  ELVMAX')
C
C     ZERO OVER OCEAN
C
      print *,' Testing at point (itest,jtest)=',itest,jtest
      print *,' SLM(itest,jtest)=',slm(itest,jtest),itest,jtest
      print *,' ORO(itest,jtest)=',oro(itest,jtest),itest,jtest
      DO J = 1,JM
        DO I = 1,IM
          IF(SLM(I,J).EQ.0.) THEN
C           VAR(I,J) = 0.
            VAR4(I,J) = 0.
            OA(I,J,1) = 0.
            OA(I,J,2) = 0.
            OA(I,J,3) = 0.
            OA(I,J,4) = 0.
            OL(I,J,1) = 0.
            OL(I,J,2) = 0.
            OL(I,J,3) = 0.
            OL(I,J,4) = 0.
C           THETA(I,J) =0.
C           GAMMA(I,J) =0.
C           SIGMA(I,J) =0.
C           ELVMAX(I,J)=0.
! --- the sub-grid scale parameters for mtn blocking and gwd retain 
! --- properties even if over ocean but there is elevation within the 
! --- gaussian grid box.
          ENDIF
       ENDDO
      ENDDO
C
      IF (MERGE_FILE == 'none') then

C  REMOVE ISOLATED POINTS
      iso_loop : DO J=2,JM-1
        JN=J-1
        JS=J+1
        DO I=1,IM
          IW=MOD(I+IM-2,IM)+1
          IE=MOD(I,IM)+1
          SLMA=SLM(IW,J)+SLM(IE,J)
          OROA=ORO(IW,J)+ORO(IE,J)
          VARA=VAR(IW,J)+VAR(IE,J)
          VAR4A=VAR4(IW,J)+VAR4(IE,J)
          DO K=1,4
            OAA(K)=OA(IW,J,K)+OA(IE,J,K)
! --- (*j*) fix typo:
            OLA(K)=OL(IW,J,K)+OL(IE,J,K)
          ENDDO
          WGTA=2
          XN=(I-1)+1
          IF(ABS(XN-NINT(XN)).LT.1.E-2) THEN
            IN=MOD(NINT(XN)-1,IM)+1
            INW=MOD(IN+IM-2,IM)+1
            INE=MOD(IN,IM)+1
            SLMA=SLMA+SLM(INW,JN)+SLM(IN,JN)+SLM(INE,JN)
            OROA=OROA+ORO(INW,JN)+ORO(IN,JN)+ORO(INE,JN)
            VARA=VARA+VAR(INW,JN)+VAR(IN,JN)+VAR(INE,JN)
            VAR4A=VAR4A+VAR4(INW,JN)+VAR4(IN,JN)+VAR4(INE,JN)
            DO K=1,4
              OAA(K)=OAA(K)+OA(INW,JN,K)+OA(IN,JN,K)+OA(INE,JN,K)
              OLA(K)=OLA(K)+OL(INW,JN,K)+OL(IN,JN,K)+OL(INE,JN,K)
            ENDDO
            WGTA=WGTA+3
          ELSE
            INW=INT(XN)
            INE=MOD(INW,IM)+1
            SLMA=SLMA+SLM(INW,JN)+SLM(INE,JN)
            OROA=OROA+ORO(INW,JN)+ORO(INE,JN)
            VARA=VARA+VAR(INW,JN)+VAR(INE,JN)
            VAR4A=VAR4A+VAR4(INW,JN)+VAR4(INE,JN)
            DO K=1,4
              OAA(K)=OAA(K)+OA(INW,JN,K)+OA(INE,JN,K)
              OLA(K)=OLA(K)+OL(INW,JN,K)+OL(INE,JN,K)
            ENDDO
            WGTA=WGTA+2
          ENDIF
          XS=(I-1)+1
          IF(ABS(XS-NINT(XS)).LT.1.E-2) THEN
            IS=MOD(NINT(XS)-1,IM)+1
            ISW=MOD(IS+IM-2,IM)+1
            ISE=MOD(IS,IM)+1
            SLMA=SLMA+SLM(ISW,JS)+SLM(IS,JS)+SLM(ISE,JS)
            OROA=OROA+ORO(ISW,JS)+ORO(IS,JS)+ORO(ISE,JS)
            VARA=VARA+VAR(ISW,JS)+VAR(IS,JS)+VAR(ISE,JS)
            VAR4A=VAR4A+VAR4(ISW,JS)+VAR4(IS,JS)+VAR4(ISE,JS)
            DO K=1,4
              OAA(K)=OAA(K)+OA(ISW,JS,K)+OA(IS,JS,K)+OA(ISE,JS,K)
              OLA(K)=OLA(K)+OL(ISW,JS,K)+OL(IS,JS,K)+OL(ISE,JS,K)
            ENDDO
            WGTA=WGTA+3
          ELSE
            ISW=INT(XS)
            ISE=MOD(ISW,IM)+1
            SLMA=SLMA+SLM(ISW,JS)+SLM(ISE,JS)
            OROA=OROA+ORO(ISW,JS)+ORO(ISE,JS)
            VARA=VARA+VAR(ISW,JS)+VAR(ISE,JS)
            VAR4A=VAR4A+VAR4(ISW,JS)+VAR4(ISE,JS)
            DO K=1,4
              OAA(K)=OAA(K)+OA(ISW,JS,K)+OA(ISE,JS,K)
              OLA(K)=OLA(K)+OL(ISW,JS,K)+OL(ISE,JS,K)
            ENDDO
            WGTA=WGTA+2
          ENDIF
          OROA=OROA/WGTA
          VARA=VARA/WGTA
          VAR4A=VAR4A/WGTA
          DO K=1,4
            OAA(K)=OAA(K)/WGTA
            OLA(K)=OLA(K)/WGTA
          ENDDO
          IF(SLM(I,J).EQ.0..AND.SLMA.EQ.WGTA) THEN
            PRINT '("SEA ",2F8.0," MODIFIED TO LAND",2F8.0," AT ",2I8)',
     &       ORO(I,J),VAR(I,J),OROA,VARA,I,J
            SLM(I,J)=1.
            ORO(I,J)=OROA
            VAR(I,J)=VARA
            VAR4(I,J)=VAR4A
            DO K=1,4
              OA(I,J,K)=OAA(K)
              OL(I,J,K)=OLA(K)
            ENDDO
          ELSEIF(SLM(I,J).EQ.1..AND.SLMA.EQ.0.) THEN
            PRINT '("LAND",2F8.0," MODIFIED TO SEA ",2F8.0," AT ",2I8)',
     &       ORO(I,J),VAR(I,J),OROA,VARA,I,J
            SLM(I,J)=0.
            ORO(I,J)=OROA
            VAR(I,J)=VARA
            VAR4(I,J)=VAR4A
            DO K=1,4
              OA(I,J,K)=OAA(K)
              OL(I,J,K)=OLA(K)
            ENDDO
          ENDIF
        ENDDO
      ENDDO iso_loop
C--- print for testing after isolated points removed
      print *,' after isolated points removed'
       call minmxj(IM,JM,ORO,'     ORO')
      print *,' ORO(itest,jtest)=',oro(itest,jtest)
      print *,' VAR(itest,jtest)=',var(itest,jtest)
      print *,' VAR4(itest,jtest)=',var4(itest,jtest)
      print *,' OA(itest,jtest,1)=',oa(itest,jtest,1)
      print *,' OA(itest,jtest,2)=',oa(itest,jtest,2)
      print *,' OA(itest,jtest,3)=',oa(itest,jtest,3)
      print *,' OA(itest,jtest,4)=',oa(itest,jtest,4)
      print *,' OL(itest,jtest,1)=',ol(itest,jtest,1)
      print *,' OL(itest,jtest,2)=',ol(itest,jtest,2)
      print *,' OL(itest,jtest,3)=',ol(itest,jtest,3)
      print *,' OL(itest,jtest,4)=',ol(itest,jtest,4)
      print *,' Testing at point (itest,jtest)=',itest,jtest
      print *,' THETA(itest,jtest)=',theta(itest,jtest)
      print *,' GAMMA(itest,jtest)=',GAMMA(itest,jtest)
      print *,' SIGMA(itest,jtest)=',SIGMA(itest,jtest)
      print *,' ELVMAX(itest,jtest)=',ELVMAX(itest,jtest)
      print *,' EFAC=',EFAC

      endif

C
      DO J=1,JM
        DO I=1,IM
          ORO(I,J) = ORO(I,J) + EFAC*VAR(I,J)
          HPRIME(I,J,1) = VAR(I,J)
          HPRIME(I,J,2) = VAR4(I,J)
          HPRIME(I,J,3) = oa(I,J,1)
          HPRIME(I,J,4) = oa(I,J,2)
          HPRIME(I,J,5) = oa(I,J,3)
          HPRIME(I,J,6) = oa(I,J,4)
          HPRIME(I,J,7) = ol(I,J,1)
          HPRIME(I,J,8) = ol(I,J,2)
          HPRIME(I,J,9) = ol(I,J,3)
          HPRIME(I,J,10)= ol(I,J,4)
          HPRIME(I,J,11)= THETA(I,J)
          HPRIME(I,J,12)= GAMMA(I,J)
          HPRIME(I,J,13)= SIGMA(I,J)
          HPRIME(I,J,14)= ELVMAX(I,J)
        ENDDO
      ENDDO
!
      deallocate(VAR4)
      deallocate (WORK1)

      call mnmxja(IM,JM,ELVMAX,itest,jtest,'  ELVMAX')
      call minmxj(IM,JM,ORO,'     ORO')

      print *,' ORO(itest,jtest),itest,jtest:',
     &          ORO(itest,jtest),itest,jtest
      print *,' ELVMAX(',itest,jtest,')=',ELVMAX(itest,jtest)

      tend=timef()
      write(6,*)' Timer 5 time= ',tend-tbeg
C
      DELXN = 360./IM

      allocate(xlat(jm), xlon(im))
      do j = 1, jm
         xlat(j) = geolat(1,j)
      enddo
      do i = 1, im
         xlon(i) = geolon(i,1)
      enddo

      tbeg=timef()
      CALL WRITE_NETCDF(IM,JM,SLM,land_frac,ORO,HPRIME,1,1,
     1                  GEOLON(1:IM,1:JM),GEOLAT(1:IM,1:JM), XLON,XLAT)
      tend=timef()
      write(6,*)' WRITE_NETCDF time= ',tend-tbeg
      print *,' wrote netcdf file out.oro.tile?.nc'

      print *,' ===== Deallocate Arrays and ENDING MTN VAR OROG program'

!     Deallocate 1d vars
      deallocate(JST,JEN)
      deallocate(XLAT,XLON,oaa,ola,GLAT)

!     Deallocate 2d vars
      deallocate (GEOLON,GEOLON_C,GEOLAT,GEOLAT_C)
      deallocate (SLM,ORO,VAR,land_frac)
      deallocate (THETA,GAMMA,SIGMA,ELVMAX)

      tend=timef()
      write(6,*)' Total runtime time= ',tend-tbeg1
      RETURN
      END SUBROUTINE TERSUB

!> Determine the location of a cubed-sphere point within
!! the high-resolution orography data.  The location is
!! described by the range of i/j indices on the high-res grid.
!!
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data set.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data set.
!! @param[in] npts Number of vertices to describe the cubed-sphere point.
!! @param[in] lonO The longitudes of the cubed-sphere vertices.
!! @param[in] latO The latitudes of the cubed-sphere vertices.
!! @param[in] delxn Resolution of the high-resolution orography
!! data set.
!! @param[out] jst Starting 'j' index on the high-resolution grid.
!! @param[out] jen Ending 'j' index on the high-resolution grid.
!! @param[out] ilist List of 'i' indices on the high-resolution grid.
!! @param[out] numx The number of 'i' indices on the high-resolution
!! grid.
!! @author GFDL programmer
      SUBROUTINE get_index(IMN,JMN,npts,lonO,latO,DELXN,
     &           jst,jen,ilist,numx)
        implicit none
        integer, intent(in)  :: IMN,JMN
        integer              :: npts
        real,    intent(in)  :: LONO(npts), LATO(npts)
        real,    intent(in)  :: DELXN
        integer, intent(out) :: jst,jen
        integer, intent(out) :: ilist(IMN)
        integer, intent(out) :: numx
        real    minlat,maxlat,minlon,maxlon
        integer i2, ii, ist, ien
        
         minlat = minval(LATO)
         maxlat = maxval(LATO)
         minlon = minval(LONO)
         maxlon = maxval(LONO)
         ist = minlon/DELXN+1
         ien = maxlon/DELXN+1
         jst = (minlat+90)/DELXN+1 
         jen = (maxlat+90)/DELXN 
         !--- add a few points to both ends of j-direction
         jst = jst - 5
         if(jst<1) jst = 1
         jen = jen + 5
         if(jen>JMN) jen = JMN

         !--- when around the pole, just search through all the points.
         if((jst == 1 .OR. jen == JMN) .and. 
     &            (ien-ist+1 > IMN/2) )then
            numx = IMN
            do i2 = 1, IMN
               ilist(i2) = i2
            enddo
         else if( ien-ist+1 > IMN/2 ) then  ! cross longitude = 0
            !--- find the minimum that greater than IMN/2 
            !--- and maximum that less than IMN/2
            ist = 0
            ien = IMN+1
            do i2 = 1, npts
               ii = LONO(i2)/DELXN+1
               if(ii <0 .or. ii>IMN) print*,"ii=",ii,IMN,LONO(i2),DELXN
               if( ii < IMN/2 ) then
                  ist = max(ist,ii)
               else if( ii > IMN/2 ) then
                  ien = min(ien,ii)
               endif
            enddo 
            if(ist<1 .OR. ist>IMN) then
               print*, "FATAL ERROR: ist<1 .or. ist>IMN"
               call ABORT()
            endif
           if(ien<1 .OR. ien>IMN) then
               print*, "FATAL ERROR: iend<1 .or. iend>IMN"
               call ABORT()
            endif

            numx = IMN - ien + 1
            do i2 = 1, numx
               ilist(i2) = ien + (i2-1)           
            enddo
            do i2 = 1, ist
               ilist(numx+i2) = i2
            enddo
            numx = numx+ist
         else
            numx = ien-ist+1
            do i2 = 1, numx
               ilist(i2) = ist + (i2-1)
            enddo
         endif

      END SUBROUTINE get_index
      
!> Create the land-mask, land fraction.
!! This routine is used for the FV3GFS model.
!!
!! @param[in] zslm The high-resolution input land-mask dataset.
!! @param[out] slm Land-mask on the model tile.
!! @param[out] land_frac Land fraction on the model tile.
!! @param[out] glat Latitude of each row of the high-resolution 
!! orography and land-mask datasets.
!! @param[in] im "i" dimension of the model grid.
!! @param[in] jm "j" dimension of the model grid.
!! @param[in] imn "i" dimension of the hi-res input orog/mask datasets.
!! @param[in] jmn "j" dimension of the hi-res input orog/mask datasets.
!! @param[in] lon_c Longitude of the model grid corner points.
!! @param[in] lat_c Latitude on the model grid corner points.
!! @author GFDL Programmer
      SUBROUTINE MAKE_MASK(zslm,SLM,land_frac,
     1 GLAT,IM,JM,IMN,JMN,lon_c,lat_c)
      implicit none
      real, parameter :: D2R = 3.14159265358979/180.
      integer, parameter :: MAXSUM=20000000
      integer IM, JM, IMN, JMN, jst, jen
      real GLAT(JMN), GLON(IMN)
      INTEGER ZSLM(IMN,JMN)
      real land_frac(IM,JM)
      real SLM(IM,JM)
      real lon_c(IM+1,JM+1), lat_c(IM+1,JM+1)
      real    LONO(4),LATO(4),LONI,LATI
      integer JM1,i,j,nsum,nsum_all,ii,jj,numx,i2
      integer ilist(IMN)
      real    DELXN,XNSUM,XLAND,XWATR,XL1,XS1,XW1
      real    XNSUM_ALL,XLAND_ALL,XWATR_ALL
      logical inside_a_polygon
C
      print *,' _____ SUBROUTINE MAKE_MASK '
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      JM1 = JM - 1
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
      ENDDO
      DO I=1,IMN
         GLON(I) = 0. + (I-1) * DELXN + DELXN * 0.5
      ENDDO
 
      land_frac(:,:) = 0.0     
C
C---- FIND THE AVERAGE OF THE MODES IN A GRID BOX
C
C  (*j*)  for hard wired zero offset (lambda s =0) for terr05 
!$omp parallel do
!$omp* private (j,i,xnsum,xland,xwatr,nsum,xl1,xs1,xw1,lono,
!$omp*          lato,jst,jen,ilist,numx,jj,i2,ii,loni,lati,
!$omp*          xnsum_all,xland_all,xwatr_all,nsum_all)
!$omp*          
      DO J=1,JM
!       print*, "J=", J
       DO I=1,IM
         XNSUM = 0.0
         XLAND = 0.0
         XWATR = 0.0
         nsum = 0
         XNSUM_ALL = 0.0
         XLAND_ALL = 0.0
         XWATR_ALL = 0.0
         nsum_all = 0
         
         LONO(1) = lon_c(i,j) 
         LONO(2) = lon_c(i+1,j) 
         LONO(3) = lon_c(i+1,j+1) 
         LONO(4) = lon_c(i,j+1) 
         LATO(1) = lat_c(i,j) 
         LATO(2) = lat_c(i+1,j) 
         LATO(3) = lat_c(i+1,j+1) 
         LATO(4) = lat_c(i,j+1) 
         call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)
         do jj = jst, jen; do i2 = 1, numx
            ii = ilist(i2)
            LONI = ii*DELXN
            LATI = -90 + jj*DELXN

            XLAND_ALL = XLAND_ALL + FLOAT(ZSLM(ii,jj))
            XWATR_ALL = XWATR_ALL + FLOAT(1-ZSLM(ii,jj))
            XNSUM_ALL = XNSUM_ALL + 1.
            nsum_all = nsum_all+1
            if(nsum_all > MAXSUM) then
              print*, "FATAL ERROR: nsum_all is greater than MAXSUM,"  
              print*, "increase MAXSUM."
              call ABORT()
            endif

            if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &          LONO*D2R,LATO*D2R))then

               XLAND = XLAND + FLOAT(ZSLM(ii,jj))
               XWATR = XWATR + FLOAT(1-ZSLM(ii,jj))
               XNSUM = XNSUM + 1.
               nsum = nsum+1
               if(nsum > MAXSUM) then
                 print*, "FATAL ERROR: nsum is greater than MAXSUM,"
                 print*, "increase MAXSUM."
                 call ABORT()
               endif
            endif
         enddo ; enddo

         
         IF(XNSUM.GT.1.) THEN
               land_frac(i,j) = XLAND/XNSUM  
               SLM(I,J) = FLOAT(NINT(XLAND/XNSUM))
         ELSEIF(XNSUM_ALL.GT.1.) THEN
               land_frac(i,j) = XLAND_ALL/XNSUM _ALL 
               SLM(I,J) = FLOAT(NINT(XLAND_ALL/XNSUM_ALL))
         ELSE
               print*, "FATAL ERROR: no source points in MAKE_MASK."
               call ABORT()
         ENDIF
       ENDDO
      ENDDO
!$omp end parallel do
      WRITE(6,*) "! MAKE_MASK DONE"
C
      RETURN
      END SUBROUTINE MAKE_MASK
!> Create the orography, land-mask, land fraction, standard 
!! deviation of orography and the convexity on a model 
!! cubed-sphere tile. This routine is used for the FV3GFS model.
!!
!! @param[in] zavg The high-resolution input orography dataset.
!! @param[in] zslm The high-resolution input land-mask dataset.
!! @param[out] oro Orography on the model tile.
!! @param[in] slm Land-mask on the model tile.
!! @param[out] var Standard deviation of orography on the model tile.
!! @param[out] var4 Convexity on the model tile.
!! @param[out] glat Latitude of each row of the high-resolution 
!! orography and land-mask datasets.
!! @param[in] im "i" dimension of the model grid.
!! @param[in] jm "j" dimension of the model grid.
!! @param[in] imn "i" dimension of the hi-res input orog/mask datasets.
!! @param[in] jmn "j" dimension of the hi-res input orog/mask datasets.
!! @param[in] lon_c Longitude of the model grid corner points.
!! @param[in] lat_c Latitude on the model grid corner points.
!! @param[in] lake_frac Fractional lake within the grid
!! @param[in] land_frac Fractional land within the grid
!! @author GFDL Programmer
      SUBROUTINE MAKEMT2(ZAVG,ZSLM,ORO,SLM,VAR,VAR4,
     1 GLAT,IM,JM,IMN,JMN,lon_c,lat_c,lake_frac,land_frac)
      implicit none
      real, parameter :: D2R = 3.14159265358979/180.
      integer, parameter :: MAXSUM=20000000
      real, dimension(:), allocatable ::  hgt_1d, hgt_1d_all
      integer IM, JM, IMN, JMN
      real GLAT(JMN), GLON(IMN)
      INTEGER ZAVG(IMN,JMN),ZSLM(IMN,JMN)
      real ORO(IM,JM),VAR(IM,JM),VAR4(IM,JM)
      real, intent(in) :: SLM(IM,JM), lake_frac(im,jm),land_frac(im,jm)
      integer JST, JEN
      real lon_c(IM+1,JM+1), lat_c(IM+1,JM+1)
      real    LONO(4),LATO(4),LONI,LATI
      real    HEIGHT
      integer JM1,i,j,nsum,nsum_all,ii,jj,i1,numx,i2
      integer ilist(IMN)
      real    DELXN,XNSUM,XLAND,XWATR,XL1,XS1,XW1,XW2,XW4
      real    XNSUM_ALL,XLAND_ALL,XWATR_ALL,HEIGHT_ALL
      real    XL1_ALL,XS1_ALL,XW1_ALL,XW2_ALL,XW4_ALL
      logical inside_a_polygon
C
      print *,' _____ SUBROUTINE MAKEMT2 '
      allocate(hgt_1d(MAXSUM))
      allocate(hgt_1d_all(MAXSUM))
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      JM1 = JM - 1
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
      ENDDO
      DO I=1,IMN
         GLON(I) = 0. + (I-1) * DELXN + DELXN * 0.5
      ENDDO
 
!     land_frac(:,:) = 0.0     
C
C---- FIND THE AVERAGE OF THE MODES IN A GRID BOX
C
C  (*j*)  for hard wired zero offset (lambda s =0) for terr05 
!$omp parallel do
!$omp* private (j,i,xnsum,xland,xwatr,nsum,xl1,xs1,xw1,xw2,xw4,lono,
!$omp*          lato,jst,jen,ilist,numx,jj,i2,ii,loni,lati,height,
!$omp*          hgt_1d,
!$omp*          xnsum_all,xland_all,xwatr_all,nsum_all,
!$omp*          xl1_all,xs1_all,xw1_all,xw2_all,xw4_all,
!$omp*          height_all,hgt_1d_all)
      DO J=1,JM
!       print*, "J=", J
       DO I=1,IM
         ORO(I,J)  = 0.0
         VAR(I,J)  = 0.0
         VAR4(I,J) = 0.0
         XNSUM = 0.0
         XLAND = 0.0
         XWATR = 0.0
         nsum = 0
         XL1 = 0.0
         XS1 = 0.0
         XW1 = 0.0
         XW2 = 0.0
         XW4 = 0.0
         XNSUM_ALL = 0.0
         XLAND_ALL = 0.0
         XWATR_ALL = 0.0
         nsum_all = 0
         XL1_ALL = 0.0
         XS1_ALL = 0.0
         XW1_ALL = 0.0
         XW2_ALL = 0.0
         XW4_ALL = 0.0
         
         LONO(1) = lon_c(i,j) 
         LONO(2) = lon_c(i+1,j) 
         LONO(3) = lon_c(i+1,j+1) 
         LONO(4) = lon_c(i,j+1) 
         LATO(1) = lat_c(i,j) 
         LATO(2) = lat_c(i+1,j) 
         LATO(3) = lat_c(i+1,j+1) 
         LATO(4) = lat_c(i,j+1) 
         call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)
         do jj = jst, jen; do i2 = 1, numx
            ii = ilist(i2)
            LONI = ii*DELXN
            LATI = -90 + jj*DELXN

            XLAND_ALL = XLAND_ALL + FLOAT(ZSLM(ii,jj))
            XWATR_ALL = XWATR_ALL + FLOAT(1-ZSLM(ii,jj))
            XNSUM_ALL = XNSUM_ALL + 1.
            HEIGHT_ALL = FLOAT(ZAVG(ii,jj)) 
            nsum_all = nsum_all+1
            if(nsum_all > MAXSUM) then
              print*, "FATAL ERROR: nsum_all is greater than MAXSUM,"
              print*, "increase MAXSUM."
              call ABORT()
            endif
            hgt_1d_all(nsum_all) = HEIGHT_ALL
            IF(HEIGHT_ALL.LT.-990.) HEIGHT_ALL = 0.0
            XL1_ALL = XL1_ALL + HEIGHT_ALL * FLOAT(ZSLM(ii,jj))
            XS1_ALL = XS1_ALL + HEIGHT_ALL * FLOAT(1-ZSLM(ii,jj))
            XW1_ALL = XW1_ALL + HEIGHT_ALL
            XW2_ALL = XW2_ALL + HEIGHT_ALL ** 2

            if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &          LONO*D2R,LATO*D2R))then

               XLAND = XLAND + FLOAT(ZSLM(ii,jj))
               XWATR = XWATR + FLOAT(1-ZSLM(ii,jj))
               XNSUM = XNSUM + 1.
               HEIGHT = FLOAT(ZAVG(ii,jj)) 
               nsum = nsum+1
               if(nsum > MAXSUM) then
                 print*, "FATAL ERROR: nsum is greater than MAXSUM,"
                 print*, "increase MAXSUM."
                 call ABORT()
               endif
               hgt_1d(nsum) = HEIGHT
               IF(HEIGHT.LT.-990.) HEIGHT = 0.0
               XL1 = XL1 + HEIGHT * FLOAT(ZSLM(ii,jj))
               XS1 = XS1 + HEIGHT * FLOAT(1-ZSLM(ii,jj))
               XW1 = XW1 + HEIGHT
               XW2 = XW2 + HEIGHT ** 2
            endif
         enddo ; enddo
         
         IF(XNSUM.GT.1.) THEN
               IF(SLM(I,J) .NE. 0. .OR. LAND_FRAC(I,J) > 0.) THEN
                  IF (XLAND > 0) THEN
                    ORO(I,J)= XL1 / XLAND
                  ELSE
                    ORO(I,J)= XS1 / XWATR
                  ENDIF
               ELSE
                  IF (XWATR > 0) THEN
                    ORO(I,J)= XS1 / XWATR
                  ELSE
                    ORO(I,J)= XL1 / XLAND
                  ENDIF
               ENDIF

               VAR(I,J)=SQRT(MAX(XW2/XNSUM-(XW1/XNSUM)**2,0.))
               do I1 = 1, NSUM
                  XW4 = XW4 + (hgt_1d(I1) - ORO(i,j)) ** 4
               enddo   

               IF(VAR(I,J).GT.1.) THEN
                  VAR4(I,J) = MIN(XW4/XNSUM/VAR(I,J) **4,10.)
               ENDIF

         ELSEIF(XNSUM_ALL.GT.1.) THEN

               !IF(SLM(I,J).NE.0.) THEN
               IF(SLM(I,J) .NE. 0. .OR. LAND_FRAC(I,J) > 0.) THEN
                  IF (XLAND_ALL > 0) THEN
                    ORO(I,J)= XL1_ALL / XLAND_ALL
                  ELSE
                    ORO(I,J)= XS1_ALL / XWATR_ALL
                  ENDIF
               ELSE
                  IF (XWATR_ALL > 0) THEN
                    ORO(I,J)= XS1_ALL / XWATR_ALL
                  ELSE
                    ORO(I,J)= XL1_ALL / XLAND_ALL
                  ENDIF
               ENDIF

               VAR(I,J)=SQRT(MAX(XW2_ALL/XNSUM_ALL-
     &                         (XW1_ALL/XNSUM_ALL)**2,0.))
               do I1 = 1, NSUM_ALL
                  XW4_ALL = XW4_ALL + 
     &                     (hgt_1d_all(I1) - ORO(i,j)) ** 4
               enddo   

               IF(VAR(I,J).GT.1.) THEN
                  VAR4(I,J) = MIN(XW4_ALL/XNSUM_ALL/VAR(I,J) **4,10.)
               ENDIF
         ELSE
               print*, "FATAL ERROR: no source points in MAKEMT2."
               call ABORT()
         ENDIF
   
! set orog to 0 meters at ocean.
!        IF (LAKE_FRAC(I,J) .EQ. 0. .AND. SLM(I,J) .EQ. 0.)THEN
         IF (LAKE_FRAC(I,J) .EQ. 0. .AND. LAND_FRAC(I,J) .EQ. 0.)THEN 
         ORO(I,J) = 0.0
         ENDIF

       ENDDO
      ENDDO
!$omp end parallel do
      WRITE(6,*) "! MAKEMT2 ORO SLM VAR VAR4 DONE"
C
      deallocate(hgt_1d)
      deallocate(hgt_1d_all)
      RETURN
      END SUBROUTINE MAKEMT2

!> Make the principle coordinates - slope of orography, 
!! anisotropy, angle of mountain range with respect to east.
!! This routine is used for the FV3GFS cubed-sphere grid.
!!
!! @param[in] zavg The high-resolution input orography dataset.
!! @param[in] zslm The high-resolution input land-mask dataset.
!! @param[out] theta Angle of mountain range with respect to
!! east for each model point.
!! @param[out] gamma Anisotropy for each model point.
!! @param[out] sigma Slope of orography for each model point.
!! @param[out] glat Latitude of each row of the high-resolution 
!! orography and land-mask datasets.
!! @param[in] im "i" dimension of the model grid tile.
!! @param[in] jm "j" dimension of the model grid tile.
!! @param[in] imn "i" dimension of the hi-res input orog/mask datasets.
!! @param[in] jmn "j" dimension of the hi-res input orog/mask datasets.
!! @param[in] lon_c Longitude of model grid corner points.
!! @param[in] lat_c Latitude of the model grid corner points.
!! @param[in] SLM mask
!! @author GFDL Programmer
      SUBROUTINE MAKEPC2(ZAVG,ZSLM,THETA,GAMMA,SIGMA,
     1           GLAT,IM,JM,IMN,JMN,lon_c,lat_c,SLM)
C
C===  PC: principal coordinates of each Z avg orog box for L&M
C
      implicit none
      real, parameter :: REARTH=6.3712E+6
      real, parameter :: D2R = 3.14159265358979/180. 
      integer :: IM,JM,IMN,JMN
      real  :: GLAT(JMN),DELTAX(JMN)
      INTEGER ZAVG(IMN,JMN),ZSLM(IMN,JMN)
      real lon_c(IM+1,JM+1), lat_c(IM+1,JM+1)
      real, intent(in) :: SLM(IM,JM)
      real HL(IM,JM),HK(IM,JM)
      real HX2(IM,JM),HY2(IM,JM),HXY(IM,JM),HLPRIM(IM,JM)
      real THETA(IM,JM),GAMMA(IM,JM),SIGMA2(IM,JM),SIGMA(IM,JM)
      real PI,CERTH,DELXN,DELTAY,XNSUM,XLAND
      real xfp,yfp,xfpyfp,xfp2,yfp2
      real hi0,hip1,hj0,hjp1,hijax,hi1j1
      real LONO(4),LATO(4),LONI,LATI
      integer i,j,i1,j1,i2,jst,jen,numx,i0,ip1,ijax
      integer ilist(IMN)
      logical inside_a_polygon
      LOGICAL DEBUG
C===  DATA DEBUG/.TRUE./
      DATA DEBUG/.FALSE./
C
      PI = 4.0 * ATAN(1.0)
      CERTH = PI * REARTH
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
      DELTAY =  CERTH / FLOAT(JMN)
      print *, 'MAKEPC2: DELTAY=',DELTAY
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
         DELTAX(J) = DELTAY * COS(GLAT(J)*D2R)
      ENDDO
C
C---- FIND THE AVERAGE OF THE MODES IN A GRID BOX
C

C... DERIVITIVE TENSOR OF HEIGHT
C
!$omp parallel do
!$omp* private (j,i,xnsum,xland,xfp,yfp,xfpyfp,
!$omp*          xfp2,yfp2,lono,lato,jst,jen,ilist,numx,j1,i2,i1,
!$omp*          loni,lati,i0,ip1,hi0,hip1,hj0,hjp1,ijax,
!$omp*          hijax,hi1j1)
      JLOOP : DO J=1,JM
!        print*, "J=", J
        ILOOP : DO I=1,IM
          HX2(I,J) = 0.0
          HY2(I,J) = 0.0
          HXY(I,J) = 0.0
          XNSUM = 0.0
          XLAND = 0.0
            xfp = 0.0
            yfp = 0.0
            xfpyfp = 0.0
            xfp2 = 0.0
            yfp2 = 0.0
            HL(I,J) = 0.0
            HK(I,J) = 0.0
            HLPRIM(I,J) = 0.0
            THETA(I,J) = 0.0 
            GAMMA(I,J) = 0.
            SIGMA2(I,J) = 0.
            SIGMA(I,J) = 0.

            LONO(1) = lon_c(i,j) 
            LONO(2) = lon_c(i+1,j) 
            LONO(3) = lon_c(i+1,j+1) 
            LONO(4) = lon_c(i,j+1) 
            LATO(1) = lat_c(i,j) 
            LATO(2) = lat_c(i+1,j) 
            LATO(3) = lat_c(i+1,j+1) 
            LATO(4) = lat_c(i,j+1) 
            call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)

            do j1 = jst, jen; do i2 = 1, numx
              i1 = ilist(i2)         
              LONI = i1*DELXN
              LATI = -90 + j1*DELXN
              INSIDE : if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &           LONO*D2R,LATO*D2R))then

C===  set the rest of the indexs for ave: 2pt staggered derivitive
C
                i0 = i1 - 1
                if (i1 - 1 .le. 0 )   i0 = i0 + imn
                if (i1 - 1 .gt. imn)  i0 = i0 - imn
C
                ip1 = i1 + 1
                if (i1 + 1 .le. 0 )   ip1 = ip1 + imn
                if (i1 + 1 .gt. imn)  ip1 = ip1 - imn

                  XLAND = XLAND + FLOAT(ZSLM(I1,J1))
                  XNSUM = XNSUM + 1.
C
                  hi0 =  float(zavg(i0,j1))
                  hip1 =  float(zavg(ip1,j1))
C
                  if(hi0 .lt. -990.)  hi0 = 0.0
                  if(hip1 .lt. -990.)  hip1 = 0.0
C........           xfp = xfp + 0.5 * ( hip1 - hi0 ) / DELTAX(J1)
                  xfp = 0.5 * ( hip1 - hi0 ) / DELTAX(J1)
                  xfp2 = xfp2 + 0.25 * ( ( hip1 - hi0 )/DELTAX(J1) )** 2 
C
! --- not at boundaries
!RAB                 if ( J1 .ne. JST(1)  .and. J1 .ne. JEN(JM) ) then
                 if ( J1 .ne. 1  .and. J1 .ne. JMN ) then
                  hj0 =  float(zavg(i1,j1-1))
                  hjp1 =  float(zavg(i1,j1+1))
                  if(hj0 .lt. -990.)  hj0 = 0.0
                  if(hjp1 .lt. -990.)  hjp1 = 0.0
C.......          yfp = yfp + 0.5 * ( hjp1 - hj0 ) / DELTAY
                  yfp = 0.5 * ( hjp1 - hj0 ) / DELTAY
                  yfp2 = yfp2 + 0.25 * ( ( hjp1 - hj0 )/DELTAY )**2   
C
C..............elseif ( J1 .eq. JST(J) .or. J1 .eq. JEN(JM) ) then
C ===     the NH pole: NB J1 goes from High at NP to Low toward SP
C
!RAB                 elseif ( J1 .eq. JST(1) ) then
                 elseif ( J1 .eq. 1 ) then
		   ijax = i1 + imn/2 
                   if (ijax .le. 0 )   ijax = ijax + imn
                   if (ijax .gt. imn)  ijax = ijax - imn
C..... at N pole we stay at the same latitude j1 but cross to opp side
                   hijax = float(zavg(ijax,j1))
                   hi1j1 = float(zavg(i1,j1))
                   if(hijax .lt. -990.)  hijax = 0.0
                   if(hi1j1 .lt. -990.)  hi1j1 = 0.0
C.......        yfp = yfp + 0.5 * ( ( 0.5 * ( hijax + hi1j1) ) - hi1j1 )/DELTAY
                   yfp = 0.5 * ( ( 0.5 * ( hijax - hi1j1 ) ) )/DELTAY
                   yfp2 = yfp2 + 0.25 * ( ( 0.5 *  ( hijax - hi1j1) ) 
     1                                              / DELTAY )**2
C
C ===     the SH pole: NB J1 goes from High at NP to Low toward SP
C
!RAB                 elseif ( J1 .eq. JEN(JM) ) then
                 elseif ( J1 .eq. JMN ) then
		  ijax = i1 + imn/2 
                  if (ijax .le. 0 )   ijax = ijax + imn
                  if (ijax .gt. imn)  ijax = ijax - imn
                  hijax = float(zavg(ijax,j1))
                  hi1j1 = float(zavg(i1,j1))
                  if(hijax  .lt. -990.)  hijax = 0.0
                  if(hi1j1  .lt. -990.)  hi1j1 = 0.0
                  if ( i1 .lt. 5 )print *,' S.Pole i1,j1 :',i1,j1,
     &                      hijax,hi1j1
C.....        yfp = yfp + 0.5 *  (0.5 * ( hijax - hi1j1) )/DELTAY  
        yfp = 0.5 *  (0.5 * ( hijax - hi1j1) )/DELTAY  
        yfp2 = yfp2 + 0.25 * (  (0.5 * (hijax - hi1j1) )
     1                                                 / DELTAY )**2  
                 endif
C
C ===    The above does an average across the pole for the bndry in j.
C23456789012345678901234567890123456789012345678901234567890123456789012......
C
                  xfpyfp = xfpyfp + xfp * yfp
               ENDIF INSIDE
C
C === average the HX2, HY2 and HXY
C === This will be done over all land
C
               ENDDO
            ENDDO
C
C ===  HTENSR 
C
         XNSUM_GT_1 : IF(XNSUM.GT.1.) THEN
               IF(SLM(I,J).NE.0.) THEN
                  IF (XLAND > 0) THEN
                    HX2(I,J) =  xfp2  / XLAND
                    HY2(I,J) =  yfp2  / XLAND
                    HXY(I,J) =  xfpyfp / XLAND
                  ELSE
                    HX2(I,J) =  xfp2  / XNSUM
                    HY2(I,J) =  yfp2  / XNSUM
                    HXY(I,J) =  xfpyfp / XNSUM
                  ENDIF
               ENDIF
C=== degub testing
      if (debug) then
          print *," I,J,i1,j1:", I,J,i1,j1,
     1         XLAND,SLM(i,j)
          print *," xfpyfp,xfp2,yfp2:",xfpyfp,xfp2,yfp2
          print *," HX2,HY2,HXY:",HX2(I,J),HY2(I,J),HXY(I,J)
      ENDIF
C
C === make the principal axes, theta, and the degree of anisotropy, 
C === and sigma2, the slope parameter
C
               HK(I,J) = 0.5 * ( HX2(I,J) + HY2(I,J) )
               HL(I,J) = 0.5 * ( HX2(I,J) - HY2(I,J) )
               HLPRIM(I,J) = SQRT(HL(I,J)*HL(I,J) + HXY(I,J)*HXY(I,J))
           IF( HL(I,J).NE. 0. .AND. SLM(I,J) .NE. 0. ) THEN
C
             THETA(I,J) = 0.5 * ATAN2(HXY(I,J),HL(I,J)) / D2R
C ===   for testing print out in degrees
C            THETA(I,J) = 0.5 * ATAN2(HXY(I,J),HL(I,J))
            ENDIF
             SIGMA2(I,J) =  ( HK(I,J) + HLPRIM(I,J) )
        if ( SIGMA2(I,J) .GE. 0. ) then 
             SIGMA(I,J) =  SQRT(SIGMA2(I,J) )
             if (sigma2(i,j) .ne. 0. .and. 
     &        HK(I,J) .GE. HLPRIM(I,J) ) 
     1       GAMMA(I,J) = sqrt( (HK(I,J) - HLPRIM(I,J)) / SIGMA2(I,J) )
        else
             SIGMA(I,J)=0.
        endif
           ENDIF XNSUM_GT_1
                  if (debug) then
       print *," I,J,THETA,SIGMA,GAMMA,",I,J,THETA(I,J),
     1                                       SIGMA(I,J),GAMMA(I,J)
       print *," HK,HL,HLPRIM:",HK(I,J),HL(I,J),HLPRIM(I,J)
                  endif
        ENDDO ILOOP
      ENDDO JLOOP
!$omp end parallel do
      WRITE(6,*) "! MAKE Principal Coord  DONE"
C
      RETURN
      END SUBROUTINE MAKEPC2
      
!> Convert the 'x' direction distance of a cubed-sphere grid
!! point to the corresponding distance in longitude.
!!
!! @param[in] dx Distance along the 'x' direction of a 
!! cubed-sphere grid point.
!! @param[in] lat Latitude of the cubed-sphere point.
!! @param[in] degrad  Conversion from radians to degrees.
!! @return get_lon_angle Corresponding distance in longitude.
!! @author GFDL programmer
      function get_lon_angle(dx,lat, DEGRAD)
      implicit none
      real dx, lat, DEGRAD
      
      real get_lon_angle
         real, parameter :: RADIUS = 6371200

         get_lon_angle = 2*asin( sin(dx/RADIUS*0.5)/cos(lat) )*DEGRAD
         
      end function get_lon_angle

!> Convert the 'y' direction distance of a cubed-sphere grid
!! point to the corresponding distance in latitude.
!!
!! @param[in] dy Distance along the 'y' direction of a cubed-sphere
!! point.
!! @param[in] degrad Conversion from radians to degrees.
!! @return get_lat_angle Corresponding distance in latitude.
!! @author GFDL programmer
      function get_lat_angle(dy, DEGRAD)
      implicit none
      real dy, DEGRAD
      
      real get_lat_angle
         real, parameter :: RADIUS = 6371200

         get_lat_angle = dy/RADIUS*DEGRAD
         
      end function get_lat_angle
      
!> Create orographic asymmetry and orographic length scale on
!! the model grid.  This routine is used for the cubed-sphere
!! grid.
!!
!! @param[in] zavg High-resolution orography data.
!! @param[in] zslm High-resolution land-mask data.
!! @param[in] var Standard deviation of orography on the model grid.
!! @param[out] glat Latitude of each row of input terrain dataset.
!! @param[out] oa4 Orographic asymmetry on the model grid. Four
!! directional components - W/S/SW/NW
!! @param[out] ol Orographic length scale on the model grid. Four
!! directional components - W/S/SW/NW
!! @param[out] ioa4 Count of oa4 values between certain thresholds.
!! @param[out] elvmax Maximum elevation within a model grid box.
!! @param[in] oro Orography on the model grid.
!! @param[out] oro1 Save array for model grid orography.
!! @param[out] xnsum Not used.
!! @param[out] xnsum1 Not used.
!! @param[out] xnsum2 Not used.
!! @param[out] xnsum3 Not used.
!! @param[out] xnsum4 Not used.
!! @param[in] im "i" dimension of the model grid tile.
!! @param[in] jm "j" dimension of the model grid tile.
!! @param[in] imn "i" dimension of the high-resolution orography and
!! mask data.
!! @param[in] jmn "j" dimension of the high-resolution orography and
!! mask data.
!! @param[in] lon_c Corner point longitudes of the model grid points.
!! @param[in] lat_c Corner point latitudes of the model grid points.
!! @param[in] lon_t Center point longitudes of the model grid points.
!! @param[in] lat_t Center point latitudes of the model grid points.
!! @param[in] dx Length of model grid points in the 'x' direction.
!! @param[in] dy Length of model grid points in the 'y' direction.
!! @param[in] is_south_pole Is the model point at the south pole?
!! @param[in] is_north_pole is the model point at the north pole?
!! @author GFDL Programmer
      SUBROUTINE MAKEOA2(ZAVG,zslm,VAR,GLAT,OA4,OL,IOA4,ELVMAX,
     1           ORO,oro1,XNSUM,XNSUM1,XNSUM2,XNSUM3,XNSUM4,
     2           IM,JM,IMN,JMN,lon_c,lat_c,lon_t,lat_t,dx,dy,
     3           is_south_pole,is_north_pole )
      implicit none
      real, parameter :: MISSING_VALUE = -9999.
      real, parameter :: D2R = 3.14159265358979/180.
      real, PARAMETER :: R2D=180./3.14159265358979
      integer IM,JM,IMN,JMN
      real    GLAT(JMN)
      INTEGER ZAVG(IMN,JMN),ZSLM(IMN,JMN)
      real    ORO(IM,JM),ORO1(IM,JM),ELVMAX(IM,JM),ZMAX(IM,JM)
      real    OA4(IM,JM,4)
      integer IOA4(IM,JM,4)
      real    lon_c(IM+1,JM+1), lat_c(IM+1,JM+1)
      real    lon_t(IM,JM), lat_t(IM,JM)
      real    dx(IM,JM), dy(IM,JM)
      logical is_south_pole(IM,JM), is_north_pole(IM,JM)
      real    XNSUM(IM,JM),XNSUM1(IM,JM),XNSUM2(IM,JM)
      real    XNSUM3(IM,JM),XNSUM4(IM,JM)
      real    VAR(IM,JM),OL(IM,JM,4)
      integer i,j,ilist(IMN),numx,i1,j1,ii1
      integer KWD
      real    LONO(4),LATO(4),LONI,LATI
      real    DELXN,HC,HEIGHT,XNPU,XNPD,T
      integer NS0,NS1,NS2,NS3,NS4,NS5,NS6
      logical inside_a_polygon
      real    lon,lat,dlon,dlat,dlat_old
      real    lon1,lat1,lon2,lat2
      real    xnsum11,xnsum12,xnsum21,xnsum22
      real    HC_11, HC_12, HC_21, HC_22
      real    xnsum1_11,xnsum1_12,xnsum1_21,xnsum1_22
      real    xnsum2_11,xnsum2_12,xnsum2_21,xnsum2_22
      real    get_lon_angle, get_lat_angle, get_xnsum
      integer jst, jen
C   
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
      ENDDO
      print *,' IM=',IM,' JM=',JM,' IMN=',IMN,' JMN=',JMN
C
C---- FIND THE AVERAGE OF THE MODES IN A GRID BOX
C
C
      DO J=1,JM
        DO I=1,IM
          XNSUM(I,J) = 0.0
          ELVMAX(I,J) = ORO(I,J)
          ZMAX(I,J)   = 0.0
C---- COUNT NUMBER OF MODE. HIGHER THAN THE HC, CRITICAL HEIGHT
C     IN A GRID BOX
          XNSUM1(I,J) = 0.0
          XNSUM2(I,J) = 0.0
          XNSUM3(I,J) = 0.0
          XNSUM4(I,J) = 0.0
          ORO1(I,J) = ORO(I,J)
          ELVMAX(I,J) = ZMAX(I,J) 
        ENDDO
      ENDDO

! --- # of peaks > ZAVG value and ZMAX(IM,JM) -- ORO is already avg.
! ---  to JM or to JM1
!$omp parallel do
!$omp* private (j,i,hc,lono,lato,jst,jen,ilist,numx,j1,ii1,i1,loni,
!$omp*          lati,height)
      DO J=1,JM
!        print*, "J=", J 
        DO I=1,IM
          HC = 1116.2 - 0.878 * VAR(I,J) 
          LONO(1) = lon_c(i,j) 
          LONO(2) = lon_c(i+1,j) 
          LONO(3) = lon_c(i+1,j+1) 
          LONO(4) = lon_c(i,j+1) 
          LATO(1) = lat_c(i,j) 
          LATO(2) = lat_c(i+1,j) 
          LATO(3) = lat_c(i+1,j+1) 
          LATO(4) = lat_c(i,j+1) 
          call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)
          do j1 = jst, jen; do ii1 = 1, numx          
            i1 = ilist(ii1)
            LONI = i1*DELXN
            LATI = -90 + j1*DELXN
            if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &          LONO*D2R,LATO*D2R))then

              HEIGHT = FLOAT(ZAVG(I1,J1))
              IF(HEIGHT.LT.-990.) HEIGHT = 0.0
              IF ( HEIGHT .gt. ORO(I,J) ) then
                 if ( HEIGHT .gt. ZMAX(I,J) )ZMAX(I,J) = HEIGHT
              ENDIF   
            endif
          ENDDO ; ENDDO
        ENDDO
      ENDDO
!$omp end parallel do      
C
! ---  this will make work1 array take on oro's values on return
! ---  this will make work1 array take on oro's values on return
      DO J=1,JM
        DO I=1,IM

          ORO1(I,J) = ORO(I,J)
          ELVMAX(I,J) = ZMAX(I,J) 
        ENDDO
      ENDDO
      
      DO KWD = 1, 4
        DO J=1,JM
          DO I=1,IM
            OA4(I,J,KWD) = 0.0
            OL(I,J,KWD) = 0.0
          ENDDO
        ENDDO
      ENDDO
                                !
! --- # of peaks > ZAVG value and ZMAX(IM,JM) -- ORO is already avg.
C
C---- CALCULATE THE 3D OROGRAPHIC ASYMMETRY FOR 4 WIND DIRECTIONS
C---- AND THE 3D OROGRAPHIC SUBGRID OROGRAPHY FRACTION
C     (KWD = 1  2  3  4)
C     ( WD = W  S SW NW)
C
C
!$omp parallel do
!$omp* private (j,i,lon,lat,kwd,dlon,dlat,lon1,lon2,lat1,lat2,
!$omp*          xnsum11,xnsum12,xnsum21,xnsum22,xnpu,xnpd,
!$omp*          xnsum1_11,xnsum2_11,hc_11, xnsum1_12,xnsum2_12,
!$omp*          hc_12,xnsum1_21,xnsum2_21,hc_21, xnsum1_22,
!$omp*          xnsum2_22,hc_22)
      DO J=1,JM
!       print*, "j = ", j
        DO I=1,IM
          lon = lon_t(i,j)
          lat = lat_t(i,j)
          !--- for around north pole, oa and ol are all 0
          
          if(is_north_pole(i,j)) then
             print*, "set oa1 = 0 and ol=0 at i,j=", i,j
             do kwd = 1, 4
                  oa4(i,j,kwd) = 0.
                  ol(i,j,kwd) = 0.
             enddo
          else if(is_south_pole(i,j)) then
             print*, "set oa1 = 0 and ol=1 at i,j=", i,j
             do kwd = 1, 4
                oa4(i,j,kwd) = 0.
                ol(i,j,kwd) = 1.
             enddo    
          else
             
          !--- for each point, find a lat-lon grid box with same dx and dy as the cubic grid box
          dlon = get_lon_angle(dx(i,j), lat*D2R, R2D  )
          dlat = get_lat_angle(dy(i,j), R2D)
          !--- adjust dlat if the points are close to pole.
          if( lat-dlat*0.5<-90.) then
             print*, "at i,j =", i,j, lat, dlat, lat-dlat*0.5
             print*, "FATAL ERROR: lat-dlat*0.5<-90."
             call ERREXIT(4)
          endif
          if( lat+dlat*2 > 90.) then
             dlat_old = dlat
             dlat = (90-lat)*0.5
             print*, "at i,j=",i,j," adjust dlat from ",
     &              dlat_old, " to ", dlat
          endif   
          !--- lower left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5

          if(lat1<-90 .or. lat2>90) then
             print*, "at upper left i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum11 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- upper left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower left i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum12 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- lower right
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at upper right i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum21 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- upper right 
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower right i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          
          xnsum22 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          
          
           XNPU = xnsum11 + xnsum12
           XNPD = xnsum21 + xnsum22
           IF (XNPD .NE. XNPU) OA4(I,J,1) = 1. - XNPD / MAX(XNPU , 1.)

           XNPU = xnsum11 + xnsum21
           XNPD = xnsum12 + xnsum22
           IF (XNPD .NE. XNPU) OA4(I,J,2) = 1. - XNPD / MAX(XNPU , 1.)

           XNPU = xnsum11 + (xnsum12+xnsum21)*0.5
           XNPD = xnsum22 + (xnsum12+xnsum21)*0.5
           IF (XNPD .NE. XNPU) OA4(I,J,3) = 1. - XNPD / MAX(XNPU , 1.)

           XNPU = xnsum12 + (xnsum11+xnsum22)*0.5
           XNPD = xnsum21 + (xnsum11+xnsum22)*0.5
           IF (XNPD .NE. XNPU) OA4(I,J,4) = 1. - XNPD / MAX(XNPU , 1.)

           
          !--- calculate OL3 and OL4
          !--- lower left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at upper left i=,j=", i, j, lat, dlat,lat1,lat2
          endif          
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_11, xnsum2_11, HC_11)          

          !--- upper left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower left i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_12, xnsum2_12, HC_12)          

          !--- lower right
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at upper right i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_21, xnsum2_21, HC_21)          

          !--- upper right 
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower right i=,j=", i, j, lat, dlat,lat1,lat2
          endif          
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_22, xnsum2_22, HC_22)           
                  
          OL(i,j,3) = (XNSUM1_22+XNSUM1_11)/(XNSUM2_22+XNSUM2_11)
          OL(i,j,4) = (XNSUM1_12+XNSUM1_21)/(XNSUM2_12+XNSUM2_21)

          !--- calculate OL1 and OL2
          !--- lower left 
          lon1 = lon-dlon*2.0
          lon2 = lon-dlon
          lat1 = lat
          lat2 = lat+dlat
          if(lat1<-90 .or. lat2>90) then
             print*, "at upper left i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_11, xnsum2_11, HC_11)          

          !--- upper left 
          lon1 = lon-dlon*2.0
          lon2 = lon-dlon
          lat1 = lat+dlat
          lat2 = lat+dlat*2.0
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower left i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_12, xnsum2_12, HC_12)          

          !--- lower right
          lon1 = lon-dlon
          lon2 = lon
          lat1 = lat
          lat2 = lat+dlat
          if(lat1<-90 .or. lat2>90) then
             print*, "at upper right i=,j=", i, j, lat, dlat,lat1,lat2
          endif          
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_21, xnsum2_21, HC_21)          

          !--- upper right 
          lon1 = lon-dlon
          lon2 = lon
          lat1 = lat+dlat
          lat2 = lat+dlat*2.0
          if(lat1<-90 .or. lat2>90) then
             print*, "at lower right i=,j=", i, j, lat, dlat,lat1,lat2
          endif
          
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_22, xnsum2_22, HC_22)           
                  
          OL(i,j,1) = (XNSUM1_11+XNSUM1_21)/(XNSUM2_11+XNSUM2_21)
          OL(i,j,2) = (XNSUM1_21+XNSUM1_22)/(XNSUM2_21+XNSUM2_22)         
          ENDIF          
        ENDDO
      ENDDO
!$omp end parallel do
      DO KWD=1,4
        DO J=1,JM
          DO I=1,IM
            T = OA4(I,J,KWD)
            OA4(I,J,KWD) = SIGN( MIN( ABS(T), 1. ), T )
          ENDDO
        ENDDO
      ENDDO
C
      NS0 = 0
      NS1 = 0
      NS2 = 0
      NS3 = 0
      NS4 = 0
      NS5 = 0
      NS6 = 0
      DO KWD=1,4
      DO J=1,JM
      DO I=1,IM
         T = ABS( OA4(I,J,KWD) )
         IF(T .EQ. 0.) THEN
            IOA4(I,J,KWD) = 0
            NS0 = NS0 + 1
         ELSE IF(T .GT. 0. .AND. T .LE. 1.) THEN
            IOA4(I,J,KWD) = 1
            NS1 = NS1 + 1
         ELSE IF(T .GT. 1. .AND. T .LE. 10.) THEN
            IOA4(I,J,KWD) = 2
            NS2 = NS2 + 1
         ELSE IF(T .GT. 10. .AND. T .LE. 100.) THEN
            IOA4(I,J,KWD) = 3
            NS3 = NS3 + 1
         ELSE IF(T .GT. 100. .AND. T .LE. 1000.) THEN
            IOA4(I,J,KWD) = 4
            NS4 = NS4 + 1
         ELSE IF(T .GT. 1000. .AND. T .LE. 10000.) THEN
            IOA4(I,J,KWD) = 5
            NS5 = NS5 + 1
         ELSE IF(T .GT. 10000.) THEN
            IOA4(I,J,KWD) = 6
            NS6 = NS6 + 1
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      WRITE(6,*) "! MAKEOA2 EXIT"
C
      RETURN

      END SUBROUTINE MAKEOA2
 
!> Print out the maximum and minimum values of
!! an array.
!!
!! @param[in] im The 'i' dimension of the array.
!! @param[in] jm The 'i' dimension of the array.
!! @param[in] a The array to check.
!! @param[in] title Name of the data to be checked.
!! @author Jordan Alpert NOAA/EMC
      SUBROUTINE minmxj(IM,JM,A,title)
      implicit none

      real A(IM,JM),rmin,rmax
      integer i,j,IM,JM
      character*8 title

      rmin=1.e+10
      rmax=-rmin
csela....................................................
csela if(rmin.eq.1.e+10)return
csela....................................................
      DO j=1,JM
      DO i=1,IM
        if(A(i,j).ge.rmax)rmax=A(i,j)
        if(A(i,j).le.rmin)rmin=A(i,j)
      ENDDO
      ENDDO
      write(6,150)rmin,rmax,title
150   format('rmin=',e13.4,2x,'rmax=',e13.4,2x,a8,' ')
C
      RETURN
      END SUBROUTINE minmxj

!> Print out the maximum and minimum values of
!! an array. Pass back the i/j location of the
!! maximum value.
!!
!! @param[in] im The 'i' dimension of the array.
!! @param[in] jm The 'i' dimension of the array.
!! @param[in] a The array to check.
!! @param[out] imax 'i' location of maximum
!! @param[out] jmax 'j' location of maximum
!! @param[in] title Name of the data to be checked.
!! @author Jordan Alpert NOAA/EMC
      SUBROUTINE mnmxja(IM,JM,A,imax,jmax,title)
      implicit none

      real A(IM,JM),rmin,rmax
      integer i,j,IM,JM,imax,jmax
      character*8 title

      rmin=1.e+10
      rmax=-rmin
csela....................................................
csela if(rmin.eq.1.e+10)return
csela....................................................
      DO j=1,JM
      DO i=1,IM
        if(A(i,j).ge.rmax)then
                         rmax=A(i,j)
                        imax=i
                        jmax=j
        endif
        if(A(i,j).le.rmin)rmin=A(i,j)
      ENDDO
      ENDDO
      write(6,150)rmin,rmax,title
150   format('rmin=',e13.4,2x,'rmax=',e13.4,2x,a8,' ')
C
      RETURN
      END SUBROUTINE mnmxja

!> Read input global 30-arc second orography data.
!!
!! @param[out] glob The orography data.
!! @author Jordan Alpert NOAA/EMC
      subroutine read_g(glob)
      implicit none

      include 'netcdf.inc'

      integer*2, intent(out) :: glob(360*120,180*120)

      integer :: ncid, error, id_var, fsize

      fsize=65536

      error=NF__OPEN("./topography.gmted2010.30s.nc",
     &               NF_NOWRITE,fsize,ncid)
      call netcdf_err(error, 'Open file topography.gmted2010.30s.nc' )
      error=nf_inq_varid(ncid, 'topo', id_var)
      call netcdf_err(error, 'Inquire varid of topo')
      error=nf_get_var_int2(ncid, id_var, glob)
      call netcdf_err(error, 'Read topo')
      error = nf_close(ncid)

      print*,' '
      call maxmin     (glob,360*120*180*120,'global0')

      return
      end subroutine read_g

!> Print the maximum, mininum, mean and
!! standard deviation of an array.
!! 
!! @param [in] ia The array to be checked.
!! @param [in] len The number of points to be checked.
!! @param [in] tile A name associated with the array.
!! @author Jordan Alpert NOAA/EMC
      subroutine maxmin(ia,len,tile)
ccmr
      implicit none
ccmr
      integer*2 ia(len)
      character*7 tile
      integer iaamax, iaamin, len, m, ja, kount
      integer(8) sum2,std,mean,isum
      integer i_count_notset,kount_9
! --- missing is -9999
c
      isum = 0
      sum2 = 0
      kount = 0
      kount_9 = 0
      iaamax = -9999999
ccmr  iaamin = 1
      iaamin =  9999999
      i_count_notset=0
           do 10 m=1,len
      ja=ia(m)
ccmr  if ( ja .lt. 0 ) print *,' ja < 0:',ja
ccmr  if ( ja .eq. -9999 ) goto 10
      if ( ja .eq. -9999 ) then
           kount_9=kount_9+1
           goto 10
      endif
      if ( ja .eq. -12345 ) i_count_notset=i_count_notset+1
ccmr  if ( ja .eq. 0 ) goto 11
      iaamax = max0( iaamax, ja )
      iaamin = min0( iaamin, ja )
!     iaamax = max0( iaamax, ia(m,j) )
!     iaamin = min0( iaamin, ia(m,j) )
  11  continue
      kount = kount + 1
      isum = isum + ja
ccmr  sum2 = sum2 + ifix( float(ja) * float(ja) )
      sum2 = sum2 + ja*ja
  10  continue
!
      mean = isum/kount
      std = ifix(sqrt(float((sum2/(kount))-mean**2)))
      print*,tile,' max=',iaamax,' min=',iaamin,' sum=',isum,
     &       ' i_count_notset=',i_count_notset
      print*,tile,' mean=',mean,' std.dev=',std,
     &       ' ko9s=',kount,kount_9,kount+kount_9 
      return
      end subroutine maxmin

!> Convert from latitude and longitude to x,y,z coordinates.
!!
!! @param[in] siz Number of points to convert.
!! @param[in] lon Longitude of points to convert.
!! @param[in] lat Latitude of points to convert.
!! @param[out] x 'x' coordinate of the converted points.
!! @param[out] y 'y' coordinate of the converted points.
!! @param[out] z 'z' coordinate of the converted points.
!! @author GFDL programmer
      subroutine latlon2xyz(siz,lon, lat, x, y, z)
      implicit none
      integer, intent(in) :: siz
      real, intent(in) :: lon(siz), lat(siz)
      real, intent(out) :: x(siz), y(siz), z(siz)
      
      integer n

      do n = 1, siz
        x(n) = cos(lat(n))*cos(lon(n))
        y(n) = cos(lat(n))*sin(lon(n))
        z(n) = sin(lat(n))
      enddo
      end subroutine latlon2xyz

!> Compute spherical angle.
!!
!! @param[in] v1 Vector 1.
!! @param[in] v2 Vector 2.
!! @param[in] v3 Vector 3.
!! @return spherical_angle Spherical Angle.
!! @author GFDL programmer
      FUNCTION spherical_angle(v1, v2, v3)
        implicit none
        real, parameter :: EPSLN30 = 1.e-30
        real, parameter :: PI=3.1415926535897931
        real v1(3), v2(3), v3(3)
        real  spherical_angle
 
        real px, py, pz, qx, qy, qz, ddd;
  
        ! vector product between v1 and v2 
        px = v1(2)*v2(3) - v1(3)*v2(2)
        py = v1(3)*v2(1) - v1(1)*v2(3)
        pz = v1(1)*v2(2) - v1(2)*v2(1)
        ! vector product between v1 and v3 
        qx = v1(2)*v3(3) - v1(3)*v3(2);
        qy = v1(3)*v3(1) - v1(1)*v3(3);
        qz = v1(1)*v3(2) - v1(2)*v3(1);

        ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
        if ( ddd <= 0.0 ) then
          spherical_angle = 0. 
        else 
          ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
          if( abs(ddd-1) < EPSLN30 ) ddd = 1;
          if( abs(ddd+1) < EPSLN30 ) ddd = -1;
          if ( ddd>1. .or. ddd<-1. ) then
            !FIX to correctly handle co-linear points (angle near pi or 0) */
            if (ddd < 0.) then
              spherical_angle = PI
            else
              spherical_angle = 0.
            endif
          else
            spherical_angle = acos( ddd )
          endif
        endif  

        return
      END FUNCTION spherical_angle
      
!> Check if a point is inside a polygon.
!!
!! @param[in] lon1 Longitude of the point to check.
!! @param[in] lat1 Latitude of the point to check.
!! @param[in] npts Number of polygon vertices.
!! @param[in] lon2 Longitude of the polygon vertices.
!! @param[in] lat2 Latitude of the polygon vertices.
!! @return inside_a_polygon When true, point is within
!! the polygon.
!! @author GFDL programmer
      FUNCTION inside_a_polygon(lon1, lat1, npts, lon2, lat2)
        implicit none
        real, parameter :: EPSLN10 = 1.e-10
        real, parameter :: EPSLN8 = 1.e-8
        real, parameter :: PI=3.1415926535897931
        real, parameter :: RANGE_CHECK_CRITERIA=0.05
        real :: anglesum, angle, spherical_angle
        integer i, ip1
        real lon1, lat1
        integer npts
        real lon2(npts), lat2(npts)
        real x2(npts), y2(npts), z2(npts)
        real lon1_1d(1), lat1_1d(1)
        real x1(1), y1(1), z1(1)
        real pnt0(3),pnt1(3),pnt2(3)
        logical inside_a_polygon
        real max_x2,min_x2,max_y2,min_y2,max_z2,min_z2
        !first convert to cartesian grid */
        call latlon2xyz(npts,lon2, lat2, x2, y2, z2);
        lon1_1d(1) = lon1
        lat1_1d(1) = lat1
        call latlon2xyz(1,lon1_1d, lat1_1d, x1, y1, z1);
        inside_a_polygon = .false.
        max_x2 = maxval(x2)
        if( x1(1) > max_x2+RANGE_CHECK_CRITERIA ) return
        min_x2 = minval(x2)
        if( x1(1)+RANGE_CHECK_CRITERIA < min_x2 ) return
        max_y2 = maxval(y2)
        if( y1(1) > max_y2+RANGE_CHECK_CRITERIA ) return
        min_y2 = minval(y2)
        if( y1(1)+RANGE_CHECK_CRITERIA < min_y2 ) return
        max_z2 = maxval(z2)
        if( z1(1) > max_z2+RANGE_CHECK_CRITERIA ) return
        min_z2 = minval(z2)
        if( z1(1)+RANGE_CHECK_CRITERIA < min_z2 ) return

        pnt0(1) = x1(1)
        pnt0(2) = y1(1)
        pnt0(3) = z1(1)
        
        anglesum = 0;
        do i = 1, npts
           if(abs(x1(1)-x2(i)) < EPSLN10 .and.
     &          abs(y1(1)-y2(i)) < EPSLN10 .and.
     &         abs(z1(1)-z2(i)) < EPSLN10 ) then ! same as the corner point
              inside_a_polygon = .true.
              return
           endif
           ip1 = i+1
           if(ip1>npts) ip1 = 1
           pnt1(1) = x2(i)
           pnt1(2) = y2(i)
           pnt1(3) = z2(i)
           pnt2(1) = x2(ip1)
           pnt2(2) = y2(ip1)
           pnt2(3) = z2(ip1)

           angle = spherical_angle(pnt0, pnt2, pnt1);
!           anglesum = anglesum + spherical_angle(pnt0, pnt2, pnt1);
           anglesum = anglesum + angle
        enddo

        if(abs(anglesum-2*PI) < EPSLN8) then
           inside_a_polygon = .true.
        else
           inside_a_polygon = .false.
        endif

        return
        
      end function inside_a_polygon

!> Count the number of high-resolution orography points that
!! are higher than the model grid box average orography height.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model 
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] zslm  The high-resolution land mask.
!! @param[in] delxn Resolution of the high-res orography data.
!! @return get_xnsum The number of high-res points above the
!! mean orography.
!! @author GFDL Programmer
      function get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,
     &                   glat,zavg,zslm,delxn)
        implicit none

        real get_xnsum
        real, intent(in) :: lon1,lat1,lon2,lat2,delxn
        integer, intent(in) ::  IMN,JMN
        real, intent(in) :: glat(JMN)
        integer, intent(in) ::  zavg(IMN,JMN),zslm(IMN,JMN)
        integer i, j, ist, ien, jst, jen, i1
        real    oro, HEIGHT
        real    xland,xwatr,xl1,xs1,slm,xnsum
        !---figure out ist,ien,jst,jen
        do j = 1, JMN
           if( GLAT(J) .GT. lat1 ) then
              jst = j
              exit
           endif
        enddo
        do j = 1, JMN
           if( GLAT(J) .GT. lat2 ) then
              jen = j
              exit
           endif
        enddo

        
        ist = lon1/delxn + 1
        ien = lon2/delxn
        if(ist .le.0) ist = ist + IMN
        if(ien < ist) ien = ien + IMN

        !--- compute average oro
          oro = 0.0
          xnsum = 0
          xland = 0
          xwatr = 0
          xl1 = 0
          xs1 = 0
          do j = jst,jen
             do i1 = 1, ien - ist + 1
                i = ist + i1 -1
                if( i .LE. 0) i = i + imn
                if( i .GT. IMN) i = i - imn
                XLAND = XLAND + FLOAT(ZSLM(I,J))
                XWATR = XWATR + FLOAT(1-ZSLM(I,J))
                XNSUM = XNSUM + 1.
                HEIGHT = FLOAT(ZAVG(I,J)) 
                IF(HEIGHT.LT.-990.) HEIGHT = 0.0
                XL1 = XL1 + HEIGHT * FLOAT(ZSLM(I,J))
                XS1 = XS1 + HEIGHT * FLOAT(1-ZSLM(I,J))
             enddo
          enddo
          if( XNSUM > 1.) THEN
             SLM = FLOAT(NINT(XLAND/XNSUM))
               IF(SLM.NE.0.) THEN
                  ORO= XL1 / XLAND
               ELSE
                  ORO = XS1 / XWATR
               ENDIF
          ENDIF
          
         get_xnsum = 0
         do j = jst, jen
            do i1= 1, ien-ist+1
               i = ist + i1 -1
               if( i .LE. 0) i = i + imn
               if( i .GT. IMN) i = i - imn
               HEIGHT = FLOAT(ZAVG(I,J))
               IF(HEIGHT.LT.-990.) HEIGHT = 0.0
               IF ( HEIGHT .gt. ORO ) get_xnsum = get_xnsum + 1
            enddo       
         enddo
         
      end function get_xnsum  
      
!> Count the number of high-resolution orography points that
!! are higher than a critical value inside a model grid box
!! (or a portion of a model grid box). The critical value is a 
!! function of the standard deviation of orography.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model 
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] delxn Resolution of the high-res orography data.
!! @param[out] xnsum1 The number of high-resolution orography
!! above the critical value inside a model grid box.
!! @param[out] xnsum2 The number of high-resolution orography
!! points inside a model grid box.
!! @param[out] hc Critical height.
!! @author GFDL Programmer
      subroutine get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,
     &                   glat,zavg,delxn,xnsum1,xnsum2,HC)
        implicit none

        real, intent(out) :: xnsum1,xnsum2,HC
        real lon1,lat1,lon2,lat2,delxn
        integer IMN,JMN
        real    glat(JMN)
        integer zavg(IMN,JMN)
        integer i, j, ist, ien, jst, jen, i1
        real    HEIGHT, var
        real    XW1,XW2,xnsum
        !---figure out ist,ien,jst,jen
        do j = 1, JMN
           if( GLAT(J) .GT. lat1 ) then
              jst = j
              exit
           endif
        enddo
        do j = 1, JMN
           if( GLAT(J) .GT. lat2 ) then
              jen = j
              exit
           endif
        enddo

        
        ist = lon1/delxn + 1
        ien = lon2/delxn
        if(ist .le.0) ist = ist + IMN
        if(ien < ist) ien = ien + IMN

        !--- compute average oro
          xnsum = 0
          XW1 = 0
          XW2 = 0
          do j = jst,jen
             do i1 = 1, ien - ist + 1
                i = ist + i1 -1
                if( i .LE. 0) i = i + imn
                if( i .GT. IMN) i = i - imn
                XNSUM = XNSUM + 1.
                HEIGHT = FLOAT(ZAVG(I,J)) 
                IF(HEIGHT.LT.-990.) HEIGHT = 0.0
                XW1 = XW1 + HEIGHT
                XW2 = XW2 + HEIGHT ** 2
             enddo
          enddo
          var = SQRT(MAX(XW2/XNSUM-(XW1/XNSUM)**2,0.))
          HC = 1116.2 - 0.878 * VAR
         xnsum1 = 0
         xnsum2 = 0
         do j = jst, jen
            do i1= 1, ien-ist+1
               i = ist + i1 -1
               if( i .LE. 0) i = i + imn
               if( i .GT. IMN) i = i - imn
               HEIGHT = FLOAT(ZAVG(I,J))
               IF ( HEIGHT .gt. HC ) xnsum1 = xnsum1 + 1
                xnsum2 = xnsum2 + 1
            enddo       
         enddo
         
      end subroutine get_xnsum2 

!> Count the number of high-resolution orography points that
!! are higher than a critical value inside a model grid box
!! (or a portion of a model grid box). Unlike routine
!! get_xnsum2(), this routine does not compute the critical
!! value. Rather, it is passed in.
!!
!! @param[in] lon1 Longitude of corner point 1 of the model 
!! grid box.
!! @param[in] lat1 Latitude of corner point 1 of the model
!! grid box.
!! @param[in] lon2 Longitude of corner point 2 of the model
!! grid box.
!! @param[in] lat2 Latitude of corner point 2 of the model
!! grid box.
!! @param[in] imn 'i' dimension of the high-resolution orography
!! data.
!! @param[in] jmn 'j' dimension of the high-resolution orography
!! data.
!! @param[in] glat Latitude of each row of the high-resolution
!! orography data.
!! @param[in] zavg The high-resolution orography.
!! @param[in] delxn Resolution of the high-res orography data.
!! @param[out] xnsum1 The number of high-resolution orography
!! above the critical value inside a model grid box.
!! @param[out] xnsum2 The number of high-resolution orography
!! points inside a model grid box.
!! @param[in] hc Critical height.
!! @author GFDL Programmer
      subroutine get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,
     &                   glat,zavg,delxn,xnsum1,xnsum2,HC)
        implicit none

        real, intent(out) :: xnsum1,xnsum2
        real lon1,lat1,lon2,lat2,delxn
        integer IMN,JMN
        real    glat(JMN)
        integer zavg(IMN,JMN)
        integer i, j, ist, ien, jst, jen, i1
        real    HEIGHT, HC
        !---figure out ist,ien,jst,jen
        ! if lat1 or lat 2 is 90 degree. set jst = JMN
        jst = JMN
        jen = JMN
        do j = 1, JMN
           if( GLAT(J) .GT. lat1 ) then
              jst = j
              exit
           endif
        enddo
        do j = 1, JMN
           if( GLAT(J) .GT. lat2 ) then
              jen = j
              exit
           endif
        enddo

        
        ist = lon1/delxn + 1
        ien = lon2/delxn
        if(ist .le.0) ist = ist + IMN
        if(ien < ist) ien = ien + IMN

         xnsum1 = 0
         xnsum2 = 0
         do j = jst, jen
            do i1= 1, ien-ist+1
               i = ist + i1 -1
               if( i .LE. 0) i = i + imn
               if( i .GT. IMN) i = i - imn
               HEIGHT = FLOAT(ZAVG(I,J))
               IF ( HEIGHT .gt. HC ) xnsum1 = xnsum1 + 1
                xnsum2 = xnsum2 + 1
            enddo       
         enddo
         
      end subroutine get_xnsum3
!> Get the date/time for the system clock.
!!
!! @author Mark Iredell
!! @return timef
      real function timef()
      character(8) :: date
      character(10) :: time
      character(5) :: zone
      integer,dimension(8) :: values
      integer :: total
      real :: elapsed
      call date_and_time(date,time,zone,values)
      total=(3600*values(5))+(60*values(6))
     *      +values(7)
      elapsed=float(total) + (1.0e-3*float(values(8)))
      timef=elapsed
      return
      end function timef

!> Find the point on the model grid tile closest to the
!! north and south pole.
!!
!! @param[in] geolat Latitude on the supergrid.
!! @param[in] nx i-dimension of the supergrid.
!! @param[in] ny j-dimension of the supergrid.
!! @param[out] i_north_pole 'i' index of north pole. '0' if
!! pole is outside of grid.
!! @param[out] j_north_pole 'j' index of north pole. '0' if
!! pole is outside of grid.
!! @param[out] i_south_pole 'i' index of south pole. '0' if
!! pole is outside of grid.
!! @param[out] j_south_pole 'j' index of south pole. '0' if
!! pole is outside of grid.
!! @author GFDL Programmer
      subroutine find_poles(geolat, nx, ny, i_north_pole, j_north_pole,
     &                     i_south_pole, j_south_pole)

      implicit none
      
      integer, intent(in)  :: nx, ny

      real, intent(in)     :: geolat(nx+1,ny+1)

      integer, intent(out) :: i_north_pole, j_north_pole
      integer, intent(out) :: i_south_pole, j_south_pole

      integer              :: i, j

      real                 :: maxlat, minlat

      print*,'- CHECK IF THE TILE CONTAINS A POLE.'
    !--- figure out pole location.
      maxlat = -90
      minlat = 90
      i_north_pole = 0
      j_north_pole = 0
      i_south_pole = 0
      j_south_pole = 0
      do j = 1, ny+1; do i = 1, nx+1
         if( geolat(i,j) > maxlat ) then
            i_north_pole=i
            j_north_pole=j
            maxlat = geolat(i,j)
         endif
         if( geolat(i,j) < minlat ) then
            i_south_pole=i
            j_south_pole=j
            minlat = geolat(i,j)
         endif
      enddo ; enddo
    !--- only when maxlat is close to 90. the point is north pole
      if(maxlat < 89.9 ) then
         i_north_pole = 0
         j_north_pole = 0
      endif
      if(minlat > -89.9 ) then
         i_south_pole = 0
         j_south_pole = 0
      endif
      print*, "- MINLAT=", minlat, "MAXLAT=", maxlat
      print*, "- NORTH POLE SUPERGRID INDEX IS ",
     &           i_north_pole, j_north_pole
      print*, "- SOUTH POLE SUPERGRID INDEX IS ",
     &           i_south_pole, j_south_pole

      end subroutine find_poles

!> Find the point on the model grid tile closest to the
!! north and south pole.
!!
!! @param[in] i_north_pole 'i' index of north pole. '0' if
!! pole is outside of grid.
!! @param[in] j_north_pole 'j' index of north pole. '0' if
!! pole is outside of grid.
!! @param[in] i_south_pole 'i' index of south pole. '0' if
!! pole is outside of grid.
!! @param[in] j_south_pole 'j' index of south pole. '0' if
!! pole is outside of grid.
!! @param[in] im i-dimension of model tile
!! @param[in] jm j-dimension of model tile
!! @param[out] is_north_pole 'true' for points surrounding the north pole.
!! @param[out] is_south_pole 'true' for points surrounding the south pole.
!! @author GFDL Programmer
      subroutine find_nearest_pole_points(i_north_pole, j_north_pole,
     &  i_south_pole, j_south_pole, im, jm, is_north_pole,
     &  is_south_pole)

      implicit none

      integer, intent(in)     :: im, jm
      integer, intent(in)     :: i_north_pole, j_north_pole
      integer, intent(in)     :: i_south_pole, j_south_pole

      logical, intent(out)    :: is_north_pole(im,jm)
      logical, intent(out)    :: is_south_pole(im,jm)

      integer                 :: i, j

      print*,'- FIND NEAREST POLE POINTS.'

      is_north_pole=.false.
      is_south_pole=.false.

      if(i_south_pole >0 .and. j_south_pole > 0) then
        if(mod(i_south_pole,2)==0) then ! stretched grid
           do j = 1, JM; do i = 1, IM
             if(i==i_south_pole/2 .and. (j==j_south_pole/2 
     &              .or. j==j_south_pole/2+1) ) then
               is_south_pole(i,j) = .true.
               print*, "- SOUTH POLE AT I,J= ", i, j
             endif
           enddo; enddo
        else
           do j = 1, JM; do i = 1, IM
              if((i==i_south_pole/2 .or. i==i_south_pole/2+1)
     &             .and. (j==j_south_pole/2 .or.
     &              j==j_south_pole/2+1) ) then
                 is_south_pole(i,j) = .true.
                 print*, "- SOUTH POLE AT I,J= ", i, j
              endif
           enddo; enddo
        endif            
      endif

      if(i_north_pole >0 .and. j_north_pole > 0) then
        if(mod(i_north_pole,2)==0) then ! stretched grid
           do j = 1, JM; do i = 1, IM
              if(i==i_north_pole/2 .and. (j==j_north_pole/2 .or.
     &                j==j_north_pole/2+1) ) then
                is_north_pole(i,j) = .true.
                print*, "- NORTH POLE AT I,J= ", i, j
              endif
           enddo; enddo
        else
           do j = 1, JM; do i = 1, IM
              if((i==i_north_pole/2 .or. i==i_north_pole/2+1)
     &             .and. (j==j_north_pole/2 .or.
     &                j==j_north_pole/2+1) ) then
                is_north_pole(i,j) = .true.
                print*, "- NORTH POLE AT I,J= ", i, j
              endif
           enddo; enddo
         endif            
      endif

      end subroutine find_nearest_pole_points
