C> @file
C> Terrain maker for the ufs weather model.
C> @author Mark Iredell @date 92-04-16
  
C> This program creates landmask, land fraction, terrain and
C> and fields required for the model's gravity wave drag 
C> (GWD) scheme.
C>
C> Specifically:
C>
C> - Land mask (yes/no flag)
C> - Land fraction
C> - Terrain (orography)
C> - Maximum elevation
C> - Standard deviation of terrain
C> - Convexity
C> - Orographic Asymetry - W/S/SW/NW directional components.
C> - Orographic Length Scale - W/S/SW/NW directional components.
C> - Anisotropy
C> - Slope of terrain
C> - Angle of mountain range with respect to East.
C>      
C> This program operates on a single cubed-sphere tile.
C>
C> Optionally, the program can compute and output only the
C> land mask and land fraction. Or, it can read in the mask
C> and fraction from an external file, then compute the 
C> terrain and GWD fields using that mask. These options
C> are used to support coupled (atm/oceann) runs of the UFS.
C> The process is:
C>   - Run this program and output the mask/fraction only.
C>   - Adjust or merge the mask/fraction with the ocean
C>     mask (using another program).
C>   - Read in this 'merged' mask/fraction and compute the
C>     terrain and GWD fields.
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
C>  - 24-08-15  Remove old code used by spectral GFS.
C>       
C>   INPUT FILES:
C>  -   UNIT5      - PROGRAM CONTROL NAMELIST.
C>  -   NCID       - MODEL 'GRID' FILE
C>  -   NCID       - GMTED2010 USGS orography (NetCDF)
C>  -   NCID       - 30" UMD land cover mask. (NetCDF)
C>  -   NCID       - GICE Grumbine 30" RAMP Antarctica orog. (NetCDF)
C>  -   NCID       - MERGE FILE. CONTAINS LAND MASK, FRACTION AND
C>                   LAKE FRACTION THAT HAS BEEN MERGED WITH AN
C>                   OCEAN GRID. (NetCDF)
C>
C>   OUTPUT FILES (ALL ON A SINGLE CUBED-SPHERE TILE) :
C>  -   NCID       - OROGRAPHY FILE (NetCDF) IF MASK_ONLY=FALSE
C>  -   NCID       - MASK FILE (NetCDF) IF MASK_ONLY=TRUE
C>  -                CONTAINS ONLY LAND MASK AND FRACTION.
C>
C> @return 0 for success, error code otherwise.
      implicit none

      character(len=256)    :: mdl_grid_file = "none"
      character(len=256)    :: external_mask_file = "none"
      integer               :: imn, jmn, im, jm, efac
      logical               :: mask_only = .false.

      print*,"- BEGIN OROGRAPHY PROGRAM."

      read(5,*) mdl_grid_file
      read(5,*) mask_only
      read(5,*) external_mask_file

      efac = 0
      imn  = 360*120
      jmn  = 180*120

      if (mask_only) then
        print*,"- WILL COMPUTE LANDMASK ONLY."
      endif

      if (trim(external_mask_file) /= "none") then
        print*,"- WILL USE EXTERNAL LANDMASK FROM FILE: ", 
     &          trim(external_mask_file)
      endif

      call read_mdl_dims(mdl_grid_file, im, jm)
         
      call tersub(imn,jmn,im,jm,efac,
     &            mdl_grid_file,mask_only,external_mask_file)

      print*,"- NORMAL TERMINATION."

      stop
      end

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
!! @param[in] EXTERNAL_MASK_FILE File containing an externally
!! generated land mask/fraction.
!! @author Jordan Alpert NOAA/EMC
      SUBROUTINE TERSUB(IMN,JMN,IM,JM,EFAC,
     &     OUTGRID,MASK_ONLY,EXTERNAL_MASK_FILE)

      use orog_utils, only : minmax, timef, remove_isolated_pts

      implicit none
      include 'netcdf.inc'
C
      integer                      :: IMN,JMN,IM,JM
      character(len=*), intent(in) :: OUTGRID
      character(len=*), intent(in) :: EXTERNAL_MASK_FILE

      logical, intent(in) :: mask_only

      integer :: efac
      integer :: i,j
      integer :: itest,jtest

      integer, allocatable :: ZAVG(:,:),ZSLM(:,:)
      integer(1), allocatable :: UMD(:,:)
      integer(2), allocatable :: glob(:,:)

      integer, allocatable :: IWORK(:,:,:)

      real :: tbeg,tend,tbeg1

      real, allocatable :: XLAT(:),XLON(:)
      real, allocatable :: GEOLON(:,:),GEOLON_C(:,:),DX(:,:)
      real, allocatable :: GEOLAT(:,:),GEOLAT_C(:,:),DY(:,:)
      real, allocatable :: SLM(:,:),ORO(:,:),VAR(:,:)
      real, allocatable :: land_frac(:,:),lake_frac(:,:)
      real, allocatable :: THETA(:,:),GAMMA(:,:),SIGMA(:,:),ELVMAX(:,:)
      real, allocatable :: VAR4(:,:)
      real, allocatable :: WORK1(:,:),WORK2(:,:),WORK3(:,:),WORK4(:,:)
      real, allocatable :: WORK5(:,:),WORK6(:,:)
      real, allocatable :: OA(:,:,:),OL(:,:,:),HPRIME(:,:,:)

      logical :: is_south_pole(IM,JM), is_north_pole(IM,JM)

      tbeg1=timef()
      tbeg=timef()

      allocate (glob(IMN,JMN))
      allocate (ZAVG(IMN,JMN))
      allocate (ZSLM(IMN,JMN))
      allocate (UMD(IMN,JMN))

! Read global mask data.

      call read_global_mask(imn,jmn,umd)
 
! Read global orography data.
 
      call read_global_orog(imn,jmn,glob)

! ZSLM initialize with all land (1). Ocean is '0'.

      ZSLM=1

! ZAVG initialize from glob

      ZAVG=glob

      do j=1,jmn
      do i=1,imn
        if ( UMD(i,j) .eq. 0 ) ZSLM(i,j) = 0
      enddo
      enddo

      deallocate (UMD,glob)

! Fixing an error in the topo 30" data set at pole (-9999).  

      do i=1,imn
        ZSLM(i,1)=0
        ZSLM(i,JMN)=1
      enddo
 
! Quality control the global topography data over Antarctica
! using RAMP data.

      call qc_orog_by_ramp(imn, jmn, zavg, zslm)

      allocate (GEOLON(IM,JM),GEOLON_C(IM+1,JM+1),DX(IM,JM))
      allocate (GEOLAT(IM,JM),GEOLAT_C(IM+1,JM+1),DY(IM,JM))
      allocate (SLM(IM,JM))
      allocate (land_frac(IM,JM),lake_frac(IM,JM))

! Reading grid file.

      call read_mdl_grid_file(outgrid,im,jm,geolon,geolon_c,
     &     geolat,geolat_c,dx,dy,is_north_pole,is_south_pole)

      tend=timef()
      print*,"- TIMING: READING INPUT DATA ",tend-tbeg
                                !
      tbeg=timef()

      IF (EXTERNAL_MASK_FILE == 'none') then
        CALL MAKE_MASK(ZSLM,SLM,land_frac,
     &       IM,JM,IMN,JMN,geolon_c,geolat_c)
        lake_frac=9999.9
      ELSE
        CALL READ_MASK(EXTERNAL_MASK_FILE,SLM,land_frac,
     &       lake_frac,im,jm)
      ENDIF

      IF (MASK_ONLY) THEN
        print*,'- WILL COMPUTE LANDMASK ONLY.'
        CALL WRITE_MASK_NETCDF(IM,JM,SLM,land_frac,
     &                         1,1,GEOLON,GEOLAT)

        DEALLOCATE(ZAVG, ZSLM, SLM, LAND_FRAC, LAKE_FRAC)
        DEALLOCATE(GEOLON, GEOLON_C, GEOLAT, GEOLAT_C)
        print*,'- NORMAL TERMINATION.'
        STOP
      END IF

      allocate (VAR(IM,JM),VAR4(IM,JM),ORO(IM,JM))

      CALL MAKEMT2(ZAVG,ZSLM,ORO,SLM,VAR,VAR4,
     &     IM,JM,IMN,JMN,geolon_c,geolat_c,lake_frac,land_frac)

      tend=timef()
      print*,"- TIMING: MASK AND OROG CREATION ", tend-tbeg

      call minmax(IM,JM,ORO,'ORO     ')
      call minmax(IM,JM,SLM,'SLM     ')
      call minmax(IM,JM,VAR,'VAR     ')
      call minmax(IM,JM,VAR4,'VAR4    ')
 
! Compute mtn principal coord HTENSR: THETA,GAMMA,SIGMA
 
      allocate (THETA(IM,JM),GAMMA(IM,JM),SIGMA(IM,JM),ELVMAX(IM,JM))

      tbeg=timef()
      CALL MAKEPC2(ZAVG,ZSLM,THETA,GAMMA,SIGMA,
     1            IM,JM,IMN,JMN,geolon_c,geolat_c,SLM)
      tend=timef()

      print*,"- TIMING: CREATE PRINCIPLE COORDINATE ",tend-tbeg

      call minmax(IM,JM,THETA,'THETA   ')
      call minmax(IM,JM,GAMMA,'GAMMA   ')
      call minmax(IM,JM,SIGMA,'SIGMA   ')
 
! COMPUTE MOUNTAIN DATA : OA OL
 
       allocate (IWORK(IM,JM,4))
       allocate (OA(IM,JM,4),OL(IM,JM,4))
       allocate (WORK1(IM,JM),WORK2(IM,JM),WORK3(IM,JM),WORK4(IM,JM))
       allocate (WORK5(IM,JM),WORK6(IM,JM))

       tbeg=timef()
       CALL MAKEOA2(ZAVG,zslm,VAR,OA,OL,IWORK,ELVMAX,ORO,
     1            WORK1,WORK2,WORK3,WORK4,WORK5,WORK6,
     2            IM,JM,IMN,JMN,geolon_c,geolat_c,
     3            geolon,geolat,dx,dy,is_south_pole,is_north_pole)
       tend=timef()

       print*,"- TIMING: CREATE ASYMETRY AND LENGTH SCALE ",tend-tbeg

!      Deallocate 2d vars
       deallocate (ZSLM,ZAVG)
       deallocate (dx,dy)
       deallocate (WORK2,WORK3,WORK4,WORK5,WORK6)

!      Deallocate 3d vars
       deallocate(IWORK)

       tbeg=timef()
       call minmax(IM,JM,OA,'OA      ')
       call minmax(IM,JM,OL,'OL      ')
       call minmax(IM,JM,ELVMAX,'ELVMAX  ')
       call minmax(IM,JM,ORO,'ORO     ')

!  Replace maximum elevation with max elevation minus orography.
!  If maximum elevation is less than the orography, replace with
!  a proxy.

      print*,"- QC MAXIMUM ELEVATION."
      DO J = 1,JM
      DO I = 1,IM
        if (ELVMAX(I,J) .lt. ORO(I,J) ) then
          ELVMAX(I,J) = MAX(  3. * VAR(I,J),0.)
        else
          ELVMAX(I,J) = MAX( ELVMAX(I,J) - ORO(I,J),0.)
        endif
      ENDDO
      ENDDO
 
      call minmax(IM,JM,ELVMAX,'ELVMAX  ',itest,jtest)
 
      print*,"- ZERO FIELDS OVER OCEAN."
 
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
 
      IF (EXTERNAL_MASK_FILE == 'none') then

        call remove_isolated_pts(im,jm,slm,oro,var,var4,oa,ol)

      endif

      allocate(hprime(im,jm,14))

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
 
      deallocate(VAR4)
      deallocate (WORK1)

      call minmax(IM,JM,ELVMAX,'ELVMAX  ',itest,jtest)
      call minmax(IM,JM,ORO,'ORO     ')

      print *,'- ORO(itest,jtest),itest,jtest:',
     &          ORO(itest,jtest),itest,jtest
      print *,'- ELVMAX(',itest,jtest,')=',ELVMAX(itest,jtest)

      tend=timef()
      print*,"- TIMING: FINAL QUALITY CONTROL ", tend-tbeg

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
      print*,"- TIMING: WRITE OUTPUT FILE ", tend-tbeg

      deallocate(XLAT,XLON)
      deallocate (GEOLON,GEOLON_C,GEOLAT,GEOLAT_C)
      deallocate (SLM,ORO,VAR,land_frac)
      deallocate (THETA,GAMMA,SIGMA,ELVMAX,HPRIME)

      tend=timef()
      print*,"- TIMING: TOTAL RUNTIME ", tend-tbeg1

      return
      END SUBROUTINE TERSUB

!> Create the land-mask, land fraction.
!! This routine is used for the FV3GFS model.
!!
!! @param[in] zslm The high-resolution input land-mask dataset.
!! @param[out] slm Land-mask on the model tile.
!! @param[out] land_frac Land fraction on the model tile.
!! @param[in] im "i" dimension of the model grid.
!! @param[in] jm "j" dimension of the model grid.
!! @param[in] imn "i" dimension of the hi-res input orog/mask datasets.
!! @param[in] jmn "j" dimension of the hi-res input orog/mask datasets.
!! @param[in] lon_c Longitude of the model grid corner points.
!! @param[in] lat_c Latitude on the model grid corner points.
!! @author GFDL Programmer
      SUBROUTINE MAKE_MASK(zslm,SLM,land_frac,
     1 IM,JM,IMN,JMN,lon_c,lat_c)
      use orog_utils, only : inside_a_polygon, get_index
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
      real    LONO_RAD(4), LATO_RAD(4)
      integer JM1,i,j,nsum,nsum_all,ii,jj,numx,i2
      integer ilist(IMN)
      real    DELXN,XNSUM,XLAND,XWATR,XL1,XS1,XW1
      real    XNSUM_ALL,XLAND_ALL,XWATR_ALL
C
      print *,'- CREATE LANDMASK AND LAND FRACTION.'
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
!$omp*          lato,lono_rad,lato_rad,jst,jen,ilist,numx,jj,i2,ii,loni,lati,
!$omp*          xnsum_all,xland_all,xwatr_all,nsum_all)
!$omp*          
      DO J=1,JM
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
         LONO_RAD=LONO*D2R
         LATO_RAD=LATO*D2R
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
     &          LONO_RAD,LATO_RAD))then

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
 
      RETURN
      END SUBROUTINE MAKE_MASK
!> Create the orography, standard deviation of orography
!! and the convexity on a model tile.
!!
!! @param[in] zavg The high-resolution input orography dataset.
!! @param[in] zslm The high-resolution input land-mask dataset.
!! @param[out] oro Orography on the model tile.
!! @param[in] slm Land-mask on the model tile.
!! @param[out] var Standard deviation of orography on the model tile.
!! @param[out] var4 Convexity on the model tile.
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
     1 IM,JM,IMN,JMN,lon_c,lat_c,lake_frac,land_frac)
      use orog_utils, only : inside_a_polygon, get_index
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
      real    LONO_RAD(4), LATO_RAD(4)
      real    HEIGHT
      integer JM1,i,j,nsum,nsum_all,ii,jj,i1,numx,i2
      integer ilist(IMN)
      real    DELXN,XNSUM,XLAND,XWATR,XL1,XS1,XW1,XW2,XW4
      real    XNSUM_ALL,XLAND_ALL,XWATR_ALL,HEIGHT_ALL
      real    XL1_ALL,XS1_ALL,XW1_ALL,XW2_ALL,XW4_ALL
C
      print*,'- CREATE OROGRAPHY AND CONVEXITY.'
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
!$omp*          lato_rad,lono_rad,hgt_1d,
!$omp*          xnsum_all,xland_all,xwatr_all,nsum_all,
!$omp*          xl1_all,xs1_all,xw1_all,xw2_all,xw4_all,
!$omp*          height_all,hgt_1d_all)
      DO J=1,JM
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
         LONO_RAD = LONO*D2R
         LATO_RAD = LATO*D2R
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
     &          LONO_RAD,LATO_RAD))then

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
!! @param[in] im "i" dimension of the model grid tile.
!! @param[in] jm "j" dimension of the model grid tile.
!! @param[in] imn "i" dimension of the hi-res input orog/mask datasets.
!! @param[in] jmn "j" dimension of the hi-res input orog/mask datasets.
!! @param[in] lon_c Longitude of model grid corner points.
!! @param[in] lat_c Latitude of the model grid corner points.
!! @param[in] SLM mask
!! @author GFDL Programmer
      SUBROUTINE MAKEPC2(ZAVG,ZSLM,THETA,GAMMA,SIGMA,
     1           IM,JM,IMN,JMN,lon_c,lat_c,SLM)
C
C===  PC: principal coordinates of each Z avg orog box for L&M
C
      use orog_utils, only : get_index, inside_a_polygon
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
      real LONO_RAD(4), LATO_RAD(4)
      integer i,j,i1,j1,i2,jst,jen,numx,i0,ip1,ijax
      integer ilist(IMN)
      LOGICAL DEBUG
C===  DATA DEBUG/.TRUE./
      DATA DEBUG/.FALSE./
C
      print*,"- CREATE PRINCIPLE COORDINATES."
      PI = 4.0 * ATAN(1.0)
      CERTH = PI * REARTH
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
      DELTAY =  CERTH / FLOAT(JMN)
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
!$omp*          hijax,hi1j1,lono_rad,lato_rad)
      JLOOP : DO J=1,JM
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
            LATO_RAD = LATO *D2R
            LONO_RAD = LONO *D2R
            call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)

            do j1 = jst, jen; do i2 = 1, numx
              i1 = ilist(i2)         
              LONI = i1*DELXN
              LATI = -90 + j1*DELXN
              INSIDE : if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &           LONO_RAD,LATO_RAD))then

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
                 elseif ( J1 .eq. JMN ) then
		  ijax = i1 + imn/2 
                  if (ijax .le. 0 )   ijax = ijax + imn
                  if (ijax .gt. imn)  ijax = ijax - imn
                  hijax = float(zavg(ijax,j1))
                  hi1j1 = float(zavg(i1,j1))
                  if(hijax  .lt. -990.)  hijax = 0.0
                  if(hi1j1  .lt. -990.)  hi1j1 = 0.0
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
 
      RETURN
      END SUBROUTINE MAKEPC2
      
!> Create orographic asymmetry and orographic length scale on
!! the model grid.  This routine is used for the cubed-sphere
!! grid.
!!
!! @param[in] zavg High-resolution orography data.
!! @param[in] zslm High-resolution land-mask data.
!! @param[in] var Standard deviation of orography on the model grid.
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
      SUBROUTINE MAKEOA2(ZAVG,zslm,VAR,OA4,OL,IOA4,ELVMAX,
     1           ORO,oro1,XNSUM,XNSUM1,XNSUM2,XNSUM3,XNSUM4,
     2           IM,JM,IMN,JMN,lon_c,lat_c,lon_t,lat_t,dx,dy,
     3           is_south_pole,is_north_pole )
      use orog_utils, only : get_lat_angle, get_lon_angle,
     &                       get_index, inside_a_polygon,
     &                       get_xnsum
      implicit none
      real, parameter :: MISSING_VALUE = -9999.
      real, parameter :: D2R = 3.14159265358979/180.
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
      real    LONO_RAD(4), LATO_RAD(4)
      real    DELXN,HC,HEIGHT,XNPU,XNPD,T
      integer NS0,NS1,NS2,NS3,NS4,NS5,NS6
      real    lon,lat,dlon,dlat,dlat_old
      real    lon1,lat1,lon2,lat2
      real    xnsum11,xnsum12,xnsum21,xnsum22
      real    HC_11, HC_12, HC_21, HC_22
      real    xnsum1_11,xnsum1_12,xnsum1_21,xnsum1_22
      real    xnsum2_11,xnsum2_12,xnsum2_21,xnsum2_22
      integer jst, jen

      print*,"- CREATE ASYMETRY AND LENGTH SCALE."
C   
C---- GLOBAL XLAT AND XLON ( DEGREE )
C
      DELXN = 360./IMN      ! MOUNTAIN DATA RESOLUTION
C
      DO J=1,JMN
         GLAT(J) = -90. + (J-1) * DELXN + DELXN * 0.5
      ENDDO
      print*,'- IM=',IM,' JM=',JM,' IMN=',IMN,' JMN=',JMN
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
!$omp*          lati,height,lono_rad,lato_rad)
      DO J=1,JM
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
          LONO_RAD = LONO * D2R
          LATO_RAD = LATO * D2R
          call get_index(IMN,JMN,4,LONO,LATO,DELXN,jst,jen,ilist,numx)
          do j1 = jst, jen; do ii1 = 1, numx          
            i1 = ilist(ii1)
            LONI = i1*DELXN
            LATI = -90 + j1*DELXN
            if(inside_a_polygon(LONI*D2R,LATI*D2R,4,
     &          LONO_RAD,LATO_RAD))then

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
        DO I=1,IM
          lon = lon_t(i,j)
          lat = lat_t(i,j)
          !--- for around north pole, oa and ol are all 0
          
          if(is_north_pole(i,j)) then
             print*, "- SET OA1 = 0 AND OL=0 AT I,J=", i,j
             do kwd = 1, 4
                  oa4(i,j,kwd) = 0.
                  ol(i,j,kwd) = 0.
             enddo
          else if(is_south_pole(i,j)) then
             print*, "- SET OA1 = 0 AND OL=1 AT I,J=", i,j
             do kwd = 1, 4
                oa4(i,j,kwd) = 0.
                ol(i,j,kwd) = 1.
             enddo    
          else
             
          !--- for each point, find a lat-lon grid box with same dx and dy as the cubic grid box
          dlon = get_lon_angle(dx(i,j), lat )
          dlat = get_lat_angle(dy(i,j))
          !--- adjust dlat if the points are close to pole.
          if( lat-dlat*0.5<-90.) then
             print*, "- AT I,J =", i,j, lat, dlat, lat-dlat*0.5
             print*, "FATAL ERROR: lat-dlat*0.5<-90."
             call ERREXIT(4)
          endif
          if( lat+dlat*2 > 90.) then
             dlat_old = dlat
             dlat = (90-lat)*0.5
             print*, "- AT I,J=",i,j," ADJUST DLAT FROM ",
     &              dlat_old, " TO ", dlat
          endif   
          !--- lower left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5

          if(lat1<-90 .or. lat2>90) then
             print*, "- AT UPPER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum11 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- upper left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum12 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- lower right
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT UPPER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          xnsum21 = get_xnsum(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,zslm,delxn)          

          !--- upper right 
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
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
             print*, "- AT UPPER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif          
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_11, xnsum2_11, HC_11)          

          !--- upper left 
          lon1 = lon-dlon*1.5
          lon2 = lon-dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_12, xnsum2_12, HC_12)          

          !--- lower right
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat-dlat*0.5
          lat2 = lat+dlat*0.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT UPPER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum2(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_21, xnsum2_21, HC_21)          

          !--- upper right 
          lon1 = lon-dlon*0.5
          lon2 = lon+dlon*0.5
          lat1 = lat+dlat*0.5
          lat2 = lat+dlat*1.5
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
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
             print*, "- AT UPPER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_11, xnsum2_11, HC_11)          

          !--- upper left 
          lon1 = lon-dlon*2.0
          lon2 = lon-dlon
          lat1 = lat+dlat
          lat2 = lat+dlat*2.0
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER LEFT I=,J=", i, j, lat, dlat,lat1,lat2
          endif
          
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_12, xnsum2_12, HC_12)          

          !--- lower right
          lon1 = lon-dlon
          lon2 = lon
          lat1 = lat
          lat2 = lat+dlat
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT UPPER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
          endif          
          call get_xnsum3(lon1,lat1,lon2,lat2,IMN,JMN,GLAt,
     &     zavg,delxn, xnsum1_21, xnsum2_21, HC_21)          

          !--- upper right 
          lon1 = lon-dlon
          lon2 = lon
          lat1 = lat+dlat
          lat2 = lat+dlat*2.0
          if(lat1<-90 .or. lat2>90) then
             print*, "- AT LOWER RIGHT I=,J=", i, j, lat, dlat,lat1,lat2
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
 
      RETURN

      END SUBROUTINE MAKEOA2
 
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
