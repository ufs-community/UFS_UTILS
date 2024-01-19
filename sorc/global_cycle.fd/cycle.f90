!> @file
!! @brief Update surface and NSST fields
!! @author Mark Iredell NCEP/EMC
  
!>  Stand alone surface/NSST cycle driver for the cubed-sphere grid.
!!  Each cubed-sphere tile runs independently on its own mpi task.
!!  The surface update component runs with threads.  The NSST
!!  update component in not threaded.
!!
!!  There are three main options (which can be called in combination):
!!  1. Update the surface fields with sfccylce (do_sfccycle = .true.)
!!  2. Update the land states with increments read in from file (do_lndinc = .true.)
!!     Designed to work with either: 
!!   2a. A land increment file  created by the GSI on the Gaussian
!!     grid. The increments are interpolated here to the model grid, using the
!!     same method as for the NST increments. This is currently implemented for
!!     applying soil temperature increments calculated from the EnKF
!!     assimilation of T2m (but this is not a requirement -  any
!!     GSI-generated soil temperature increment file can be applied here).
!!   2b. A land increment file created by JEDI, on the native model grid 
!!      (cube sphere tiles). This is currently implemented for snow depth 
!!      updates for the Noah model. 
!!  3. Update the NSST field, several options:
!!
!!  3a. Update the NSST TREF field using
!!     GSI increments on the Gaussian grid.  All other NSST
!!     fields are cycled.  Invoke this option by setting
!!     namelist variable DONST=.true. and NST_FILE to
!!     the name of the GSI increment file.
!!
!!  3b. Run with NSST, but postpone the TREF update.
!!     Here all NSST fields are cycled.  But the NSST IFD field is
!!     used to flag points that flipped from ice to open water.
!!     To invoke this option, set DONST=.true. and NST_FILE="NULL".
!!
!!  INPUT FILES:
!!  - fngrid.$NNN      The cubed-sphere grid file (contains
!!                     grid point latitude and longitdue).
!!  - fnorog.$NNN      The cubed-sphere orography file (contains
!!                     land mask and orography).
!!  - fnbgsi.$NNN      The cubed-sphere input sfc/nsst restart
!!                     file.
!!  - $NST_FILE        Gaussian GSI file which contains NSST
!!                     TREF increments
!!  - $LND_SOI_FILE.$NNN    Gaussian GSI file which contains soil state
!!                     increments
!!  - xainc.$NNN       The cubed-sphere increment file (contains 
!!                     increments calculated by JEDI on the native 
!!                     model grid). 
!!  
!!  OUTPUT FILES:
!!  - fnbgso.$NNN        The updated sfc/nsst restart file.
!!
!!  NOTE: $NNN corresponds to (mpi rank + 1)
!!
!!  NAMELIST VARIABLE DEFINITIONS:
!!
!!  - IDIM,JDIM    i/j dimension of a cubed-sphere tile.
!!  - LUGB         Unit number used in the sfccycle subprogram
!!                 to read input datasets.
!!  - LSOIL        Number of soil layers.
!!  - IY,IM,ID,IH  Year, month, day, and hour of initial state.
!!  - FH           Forecast hour
!!  - DELTSFC      Cycling frequency in hours.
!!  - IALB         Use modis albedo when '1'. Use brigleb when '0'.
!!  - USE_UFO      Adjust sst and soil substrate temperature for
!!                 differences between the filtered and unfiltered
!!                 terrain.
!!  -DONST         Process NSST records.
!!  -DO_SFCCYCLE   Call sfccycle routine to update surface fields
!!  -DO_LNDINC     Read in land increment files, and add increments to
!!                 relevant states.
!!  -DO_SOI_INC     Do land increments to soil states.
!!  -DO_SNO_INC     Do land increments to snow states.
!!  - ISOT         Use statsgo soil type when '1'. Use zobler when '0'.
!!  - IVEGSRC      Use igbp veg type when '1'.  Use sib when '2'.
!!  - ZSEA1/2_MM   When running with NSST model, this is the lower/
!!                 upper bound of depth of sea temperature.  In
!!                 whole mm.
!!  - MAX_TASKS    Normally, program should be run with a number of mpi 
!!                 tasks equal to the number of cubed-sphere tiles 
!!                 being processed. However, the current parallel 
!!                 scripts may over-specify the number of tasks.
!!                 Set this variable to not process any ranks >
!!                 (max_tasks-1).
!!  -NST_FILE       path/name of the gaussian GSI file which contains NSST
!!                 TREF increments.
!!  -LND_SOI_FILE  path/name of the gaussian GSI file which contains soil
!!                 state increments.
!!
!!  -2005-02-03:  Iredell   for global_analysis
!!  -2014-11-30:  xuli      add nst_anl
!!  -2015-05-26:  Hang Lei  Added NEMSIO read/write function in the code
!!  -2017-08-08:  Gayno     Modify to work on cubed-sphere grid.
!!                         Added processing of NSST and TREF update.
!!                         Added mpi directives.
!!  -2020-02-17:  Clara Draper Added soil state increments capability.
!!
!! @author Mark Iredell NOAA/EMC
!! @return 0 for success, error code otherwise.
 PROGRAM SFC_DRV

 use mpi

 IMPLICIT NONE
!
 CHARACTER(LEN=3) :: DONST
 INTEGER :: IDIM, JDIM, LSOIL, LUGB, IY, IM, ID, IH, IALB
 INTEGER :: ISOT, IVEGSRC, LENSFC, ZSEA1_MM, ZSEA2_MM, IERR
 INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS, MAX_TASKS
 REAL    :: FH, DELTSFC, ZSEA1, ZSEA2
 LOGICAL :: USE_UFO, DO_NSST, DO_LNDINC, DO_SFCCYCLE
!
 NAMELIST/NAMCYC/ IDIM,JDIM,LSOIL,LUGB,IY,IM,ID,IH,FH,&
                  DELTSFC,IALB,USE_UFO,DONST,             &
                  DO_SFCCYCLE,ISOT,IVEGSRC,ZSEA1_MM,      &
                  ZSEA2_MM, MAX_TASKS, DO_LNDINC
!
 DATA IDIM,JDIM,LSOIL/96,96,4/
 DATA IY,IM,ID,IH,FH/1997,8,2,0,0./
 DATA LUGB/51/, DELTSFC/0.0/, IALB/1/, MAX_TASKS/99999/
 DATA ISOT/1/, IVEGSRC/2/, ZSEA1_MM/0/, ZSEA2_MM/0/
!
 CALL MPI_INIT(IERR)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

 if (myrank==0) call w3tagb('GLOBAL_CYCLE',2018,0179,0055,'NP20')

 NUM_THREADS = NUM_PARTHDS()

 PRINT*
 PRINT*,"STARTING CYCLE PROGRAM ON RANK ", MYRANK
 PRINT*,"RUNNING WITH ", NPROCS, "TASKS"
 PRINT*,"AND WITH ", NUM_THREADS, " THREADS."

 USE_UFO = .FALSE.
 DONST   = "NO"
 DO_LNDINC   = .FALSE.
 DO_SFCCYCLE = .TRUE.

 PRINT*
 PRINT*,"READ NAMCYC NAMELIST."

 CALL BAOPENR(36, "fort.36", IERR)
 READ(36, NML=NAMCYC)
 !IF (MYRANK==0) WRITE(6,NAMCYC)

 IF (MAX_TASKS < 99999 .AND. MYRANK > (MAX_TASKS - 1)) THEN
   PRINT*,"USER SPECIFIED MAX NUMBER OF TASKS: ", MAX_TASKS
   PRINT*,"WILL NOT RUN CYCLE PROGRAM ON RANK: ", MYRANK
   GOTO 333
 ENDIF

 LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE

 ZSEA1 = FLOAT(ZSEA1_MM) / 1000.0  ! CONVERT FROM MM TO METERS
 ZSEA2 = FLOAT(ZSEA2_MM) / 1000.0

 IF (DONST == "YES") THEN
   DO_NSST=.TRUE.
 ELSE
   DO_NSST=.FALSE.
 ENDIF

 PRINT*
 IF (MYRANK==0) PRINT*,"LUGB,IDIM,JDIM,ISOT,IVEGSRC,LSOIL,DELTSFC,IY,IM,ID,IH,FH: ", &
              LUGB,IDIM,JDIM,ISOT,IVEGSRC,LSOIL,DELTSFC,IY,IM,ID,IH,FH

 CALL SFCDRV(LUGB,IDIM,JDIM,LENSFC,LSOIL,DELTSFC,  &
             IY,IM,ID,IH,FH,IALB,                  &
             USE_UFO,DO_NSST,DO_SFCCYCLE,DO_LNDINC, &
             ZSEA1,ZSEA2,ISOT,IVEGSRC,MYRANK)
 
 PRINT*
 PRINT*,'CYCLE PROGRAM COMPLETED NORMALLY ON RANK: ', MYRANK

 333 CONTINUE

 CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

 if (myrank==0) call w3tage('GLOBAL_CYCLE')

 CALL MPI_FINALIZE(IERR)

 STOP

 END PROGRAM SFC_DRV

 !> Driver routine for updating the surface fields.
 !!
 !!  This program runs in two different modes:
 !!
 !!  1.  Analysis mode (FH=0.)
 !!
 !!      This program merges climatology, analysis and forecast guess to create
 !!      new surface fields.  If analysis file is given, the program 
 !!      uses it if date of the analysis matches with IY,IM,ID,IH (see Note
 !!      below).
 !!
 !!  2.  Forecast mode (FH.GT.0.)
 !!   
 !!      This program interpolates climatology to the date corresponding to the 
 !!      forecast hour.  If surface analysis file is given, for the corresponding
 !!      dates, the program will use it.  This is forcing-by-observation experiment.
 !!
 !!  If the date of the analysis does not match given IY,IM,ID,IH, (and FH),
 !!  the program searches an old analysis by going back 6 hours, then 12 hours,
 !!  then one day upto NREPMX days (parameter statement in the SUBROTINE FIXRD. 
 !!  Now defined as 15).  This allows the user to provide non-daily analysis to 
 !!  be used.  If matching field is not found, the forecast guess will be used.
 !!
 !!  Variable naming convention for this program:
 !!
 !!   - OROG .. Orography
 !!   - ALB  .. Snow-free albedo
 !!   - SWE  .. Snow water equivalent
 !!   - ZOR  .. Surface roughness length
 !!   - VET  .. Vegetation type
 !!   - TSF  .. Surface skin temperature.  Sea surface temp. over ocean.
 !!   - TG3  .. Deep soil temperature (at 500cm)
 !!   - STC  .. Soil temperature (LSOIL layrs)
 !!   - SMC  .. Total soil moisture (LSOIL layrs)
 !!   - AIS  .. Sea ice mask (0 or 1)
 !!   - CNP  .. Canopy water content
 !!   - CV   .. Convective cloud cover
 !!   - CVB  .. Convective cloud base
 !!   - CVT  .. Convective cloud top
 !!   - SLI  .. LAND/SEA/SEA-ICE mask. (1/0/2 respectively)
 !!   - VEG  .. Vegetation cover
 !!   - SOT  .. Soil type
 !!   - SIH  .. Sea ice thickness
 !!   - SIC  .. Sea ice concentration
 !!   - SND  .. Snow depth
 !!   - SLC  .. Liquid soil moisture (LSOIL layers)
 !!   - VMN  .. Vegetation cover minimum
 !!   - VMX  .. Vegetation cover maximum
 !!   - SLP  .. Slope type
 !!   - ABS  .. Maximum snow albedo
 !!   - T2M  .. 2m Temperature
 !!   - Q2M  .. 2m Specific Humidity
 !!   - TICE .. Ice Temperature
 !!   - OROG_UF .. Orography unfiltered
 !!
 !! Most fields have a blending coefficient. This controls
 !! the blending of the forecast (first guess) and interpolated
 !! climatology or analyzed fields. When it is equal to 1.0, the
 !! pure forecast is used. When the coefficient is equal to 0,
 !! the pure climatology or analysis is used. The default values
 !! are set as follows:
 !!
 !!   Variables |  Land  |  Sea
 !!   ----------|--------|-------------------------------------
 !!   Surface temperature    |   Forecast          |  Analysis
 !!   Albedo                 |   Analysis          |  Analysis
 !!   Sea-ice                |   Analysis          |  Analysis
 !!   Snow                   |   Analysis          |  Forecast (over sea ice)
 !!   Roughness              |   Analysis          |  Forecast
 !!   Plant resistance       |   Analysis          |  Analysis
 !!   Soil moisture          |   Weighted average  |  Analysis
 !!   Soil temperature       |   Forecast          |  Analysis
 !!   Canopy waver content   |   Forecast          |  Forecast
 !!   Convective cloud cover |   Forecast          |  Forecast
 !!   Convective cloud bottm |   Forecast          |  Forecast
 !!   Convective cloud top   |   Forecast          |  Forecast
 !!   Vegetation greenness   |   Analysis          |  Analysis
 !!   Vegetation type        |   Analysis          |  Analysis
 !!   Soil type              |   Analysis          |  Analysis
 !!
 !! @param[in] LUGB Fortran unit number uses in sfccycle subprogram to
 !!            read input datasets.
 !! @param[in] IDIM 'i' dimension of the cubed-sphere tile
 !! @param[in] JDIM 'j' dimension of the cubed-sphere tile
 !! @param[in] LENSFC Total numberof points for the cubed-sphere tile
 !! @param[in] LSOIL Number of soil layers
 !! @param[in] DELTSFC Cycling frequency in hours
 !! @param[in] IY Year of initial state
 !! @param[in] IM Month of initial state
 !! @param[in] ID Day of initial state
 !! @param[in] IH Hour of initial state
 !! @param[in] FH Forecast hour
 !! @param[in] IALB Use modis albedo when '1'. Use brigleb when '0'.
 !! @param[in] USE_UFO When true, adjust SST and soil temperature for
 !!            differences between the filtered and unfiltered terrain.
 !! @param[in] DO_NSST When true, process NSST records.
 !! @param[in] DO_SFCCYCLE Call sfccycle routine to update surface fields
 !! @param[in] DO_LNDINC Read in land increment files, and add increments to
 !!            requested states.
 !! @param[in] ZSEA1 When running NSST model, this is the lower bound
 !!            of depth of sea temperature.  In whole mm.
 !! @param[in] ZSEA2 When running NSST model, this is the upper bound
 !!            of depth of sea temperature.  In whole mm.
 !! @param[in] ISOT Use STATSGO soil type when '1'.  Use Zobler when '0'.
 !! @param[in] IVEGSRC Use IGBP vegetation type when '1'.  Use SIB when '2'.
 !! @param[in] MYRANK MPI rank number
 !! @author Mark Iredell, George Gayno
 SUBROUTINE SFCDRV(LUGB, IDIM,JDIM,LENSFC,LSOIL,DELTSFC,  &
                   IY,IM,ID,IH,FH,IALB,                  &
                   USE_UFO,DO_NSST,DO_SFCCYCLE,DO_LNDINC,&
                   ZSEA1,ZSEA2,ISOT,IVEGSRC,MYRANK)
!
 USE READ_WRITE_DATA
 use machine
 USE MPI
 USE LAND_INCREMENTS, ONLY: ADD_INCREMENT_SOIL,     &
                            ADD_INCREMENT_SNOW,     &
                            CALCULATE_LANDINC_MASK, &
                            APPLY_LAND_DA_ADJUSTMENTS_SOIL, &
                            APPLY_LAND_DA_ADJUSTMENTS_SND, &
                            LSM_NOAH, LSM_NOAHMP

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: IDIM, JDIM, LENSFC, LSOIL, IALB
 INTEGER, INTENT(IN) :: LUGB, IY, IM, ID, IH
 INTEGER, INTENT(IN) :: ISOT, IVEGSRC, MYRANK

 LOGICAL, INTENT(IN) :: USE_UFO, DO_NSST,DO_SFCCYCLE
 LOGICAL, INTENT(IN) :: DO_LNDINC
 
 REAL, INTENT(IN)    :: FH, DELTSFC, ZSEA1, ZSEA2

 INTEGER, PARAMETER  :: NLUNIT=35
 INTEGER, PARAMETER  :: SZ_NML=1

 CHARACTER(LEN=5)    :: TILE_NUM
 CHARACTER(LEN=500)  :: NST_FILE
 CHARACTER(LEN=500)  :: LND_SOI_FILE
 CHARACTER(LEN=4)    :: INPUT_NML_FILE(SZ_NML)

 INTEGER             :: I, IERR
 INTEGER             :: I_INDEX(LENSFC), J_INDEX(LENSFC)
 INTEGER             :: IDUM(IDIM,JDIM)
 integer             :: num_parthds, num_threads

 LOGICAL             :: IS_NOAHMP
 INTEGER             :: LSM

 real(kind=kind_io8) :: min_ice(lensfc)

 REAL                :: SLMASK(LENSFC), OROG(LENSFC)
 REAL                :: SIHFCS(LENSFC), SICFCS(LENSFC)
 REAL                :: SITFCS(LENSFC), TSFFCS(LENSFC)
 REAL                :: SWEFCS(LENSFC), ZORFCS(LENSFC)
 REAL                :: ALBFCS(LENSFC,4), TG3FCS(LENSFC)
 REAL                :: CNPFCS(LENSFC), SMCFCS(LENSFC,LSOIL)
 REAL                :: STCFCS(LENSFC,LSOIL), SLIFCS(LENSFC)
 REAL                :: AISFCS(LENSFC), F10M(LENSFC)
 REAL                :: VEGFCS(LENSFC), VETFCS(LENSFC)
 REAL                :: SOTFCS(LENSFC), ALFFCS(LENSFC,2)
 REAL                :: CVFCS(LENSFC), CVTFCS(LENSFC)
 REAL                :: CVBFCS(LENSFC), TPRCP(LENSFC)
 REAL                :: SRFLAG(LENSFC), SNDFCS(LENSFC)
 REAL                :: SLCFCS(LENSFC,LSOIL), VMXFCS(LENSFC)
 REAL                :: VMNFCS(LENSFC), T2M(LENSFC)
 REAL                :: Q2M(LENSFC), SLPFCS(LENSFC)
 REAL                :: ABSFCS(LENSFC), OROG_UF(LENSFC)
 REAL                :: USTAR(LENSFC)
 REAL                :: FMM(LENSFC), FHH(LENSFC)
 REAL                :: RLA(LENSFC), RLO(LENSFC)
 REAL(KIND=4)        :: ZSOIL(LSOIL)
 REAL                :: SIG1T(LENSFC) !< The sigma level 1 temperature for a
                                      !! dead start. Set to zero for non-dead
                                      !! start.
 REAL, ALLOCATABLE   :: STC_BCK(:,:), SMC_BCK(:,:), SLC_BCK(:,:)
 REAL, ALLOCATABLE   :: SLIFCS_FG(:)
 INTEGER, ALLOCATABLE :: LANDINC_MASK_FG(:), LANDINC_MASK(:)
 REAL, ALLOCATABLE   :: SND_BCK(:), SND_INC(:), SWE_BCK(:)
 REAL(KIND=KIND_IO8), ALLOCATABLE :: SLMASKL(:), SLMASKW(:)

 TYPE(NSST_DATA)     :: NSST
 real, dimension(idim,jdim) :: tf_clm,tf_trd,sal_clm
 real, dimension(lensfc)    :: tf_clm_tile,tf_trd_tile,sal_clm_tile
 INTEGER             :: veg_type_landice
 INTEGER, DIMENSION(LENSFC) :: STC_UPDATED, SLC_UPDATED

 LOGICAL :: FILE_EXISTS, DO_SOI_INC, DO_SNO_INC
 CHARACTER(LEN=3)       :: RANKCH

!--------------------------------------------------------------------------------
! NST_FILE is the path/name of the gaussian GSI file which contains NSST
! increments.
!--------------------------------------------------------------------------------
 
 NAMELIST/NAMSFCD/ NST_FILE, LND_SOI_FILE, DO_SNO_INC

 DATA NST_FILE/'NULL'/
 DATA LND_SOI_FILE/'NULL'/

 DO_SNO_INC = .FALSE.
 DO_SOI_INC = .FALSE.
 

 SIG1T = 0.0            ! Not a dead start!

 INPUT_NML_FILE = "NULL"

 CALL BAOPENR(37, "fort.37", IERR)
 READ (37, NML=NAMSFCD)

 PRINT*
 PRINT*,'IN ROUTINE SFCDRV,IDIM=',IDIM,'JDIM=',JDIM,'FH=',FH

!--------------------------------------------------------------------------------
! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
!--------------------------------------------------------------------------------

 CALL READ_LAT_LON_OROG(RLA,RLO,OROG,OROG_UF,TILE_NUM,IDIM,JDIM,LENSFC)

 DO I = 1, IDIM
   IDUM(I,:) = I
 ENDDO

 I_INDEX = RESHAPE(IDUM, (/LENSFC/))

 DO I = 1, JDIM
   IDUM(:,I) = I
 ENDDO

 J_INDEX = RESHAPE(IDUM, (/LENSFC/) )

 IF (DO_NSST) THEN
   PRINT*
   PRINT*,"WILL PROCESS NSST RECORDS."
   ALLOCATE(NSST%C_0(LENSFC))
   ALLOCATE(NSST%C_D(LENSFC))
   ALLOCATE(NSST%D_CONV(LENSFC))
   ALLOCATE(NSST%DT_COOL(LENSFC))
   ALLOCATE(NSST%IFD(LENSFC))
   ALLOCATE(NSST%QRAIN(LENSFC))
   ALLOCATE(NSST%TREF(LENSFC))
   ALLOCATE(NSST%TFINC(LENSFC))
   ALLOCATE(NSST%W_0(LENSFC))
   ALLOCATE(NSST%W_D(LENSFC))
   ALLOCATE(NSST%XS(LENSFC))
   ALLOCATE(NSST%XT(LENSFC))
   ALLOCATE(NSST%XTTS(LENSFC))
   ALLOCATE(NSST%XU(LENSFC))
   ALLOCATE(NSST%XV(LENSFC))
   ALLOCATE(NSST%XZ(LENSFC))
   ALLOCATE(NSST%XZTS(LENSFC))
   ALLOCATE(NSST%Z_C(LENSFC))
   ALLOCATE(NSST%ZM(LENSFC))
   ALLOCATE(SLIFCS_FG(LENSFC))
 ENDIF

IF (DO_LNDINC) THEN
   ! identify variables to be updated, and allocate arrays.
   IF  (TRIM(LND_SOI_FILE) .NE. "NULL") THEN
       DO_SOI_INC = .TRUE.
       PRINT*
       PRINT*," APPLYING SOIL INCREMENTS FROM THE GSI"
       ALLOCATE(STC_BCK(LENSFC, LSOIL), SMC_BCK(LENSFC, LSOIL), SLC_BCK(LENSFC,LSOIL))
       ALLOCATE(LANDINC_MASK_FG(LENSFC))
   ENDIF
   ! FOR NOW, CODE SO CAN DO BOTH, BUT MIGHT NEED TO THINK ABOUT THIS SOME MORE.
   IF  (DO_SNO_INC) THEN
       ! ideally, would check here that sfcsub snow DA update is not also requested
       ! but latter is controlled by fnsol, which is read in within that routine.
       ! should be done at script level.
       PRINT*
       PRINT*," APPLYING SNOW INCREMENTS FROM JEDI"
       ALLOCATE(SND_BCK(LENSFC), SND_INC(LENSFC), SWE_BCK(LENSFC))
   ENDIF
   ! set-up land mask info
   ALLOCATE(LANDINC_MASK(LENSFC))
   if (ivegsrc == 2) then   ! sib
      veg_type_landice=13
   else
      veg_type_landice=15
   endif
ENDIF

!--------------------------------------------------------------------------------
! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE.
!--------------------------------------------------------------------------------

 CALL READ_DATA(LSOIL,LENSFC,DO_NSST,.false.,IS_NOAHMP=IS_NOAHMP, &
                TSFFCS=TSFFCS,SMCFCS=SMCFCS,   &
                SWEFCS=SWEFCS,STCFCS=STCFCS,TG3FCS=TG3FCS,ZORFCS=ZORFCS,  &
                CVFCS=CVFCS,  CVBFCS=CVBFCS,CVTFCS=CVTFCS,ALBFCS=ALBFCS,  &
                VEGFCS=VEGFCS,SLIFCS=SLIFCS,CNPFCS=CNPFCS,F10M=F10M    ,  &
                VETFCS=VETFCS,SOTFCS=SOTFCS,ALFFCS=ALFFCS,USTAR=USTAR  ,  &
                FMM=FMM      ,FHH=FHH      ,SIHFCS=SIHFCS,SICFCS=SICFCS,  &
                SITFCS=SITFCS,TPRCP=TPRCP  ,SRFLAG=SRFLAG,SNDFCS=SNDFCS,  &
                VMNFCS=VMNFCS,VMXFCS=VMXFCS,SLCFCS=SLCFCS,SLPFCS=SLPFCS,  &
                ABSFCS=ABSFCS,T2M=T2M      ,Q2M=Q2M      ,SLMASK=SLMASK,  &
                ZSOIL=ZSOIL,   NSST=NSST)

 IF (IS_NOAHMP) THEN 
        LSM=LSM_NOAHMP
 ELSE
        LSM=LSM_NOAH
 ENDIF

 IF (USE_UFO) THEN
   PRINT*
   PRINT*,'USE UNFILTERED OROGRAPHY.'
 ELSE
   OROG_UF = 0.0
 ENDIF
 
 DO I=1,LENSFC
   AISFCS(I) = 0.
   IF(NINT(SLIFCS(I)).EQ.2) AISFCS(I) = 1.
 ENDDO

 IF (DO_NSST) THEN
   IF (.NOT. DO_SFCCYCLE ) THEN
     PRINT*
     PRINT*,"FIRST GUESS MASK ADJUSTED BY IFD RECORD"
     SLIFCS_FG = SLIFCS
     WHERE(NINT(NSST%IFD) == 3) SLIFCS_FG = 2.0
   ELSE
     PRINT*
     PRINT*,"SAVE FIRST GUESS MASK"
     SLIFCS_FG = SLIFCS
   ENDIF
 ENDIF
 
 ! CALCULATE MASK FOR LAND INCREMENTS
 IF (DO_LNDINC)  &
    CALL CALCULATE_LANDINC_MASK(SMCFCS(:,1),SWEFCS, VETFCS,  &
                    LENSFC,VEG_TYPE_LANDICE,  LANDINC_MASK)

!--------------------------------------------------------------------------------
! UPDATE SURFACE FIELDS.
!
! FIRST, SET WATER AND LAND MASKS - SLMASKW/SLMASKL. FOR UNCOUPLED
! (NON-FRACTIONAL) MODE, THESE ARE IDENTICAL TO THE CURRENT
! MASK - '0' WATER; '1' LAND.
!--------------------------------------------------------------------------------

 IF (DO_SFCCYCLE) THEN
   ALLOCATE(SLMASKL(LENSFC), SLMASKW(LENSFC))
! for running uncoupled (non-fractional grid)
   DO I=1,LENSFC
     IF(NINT(SLMASK(I)) == 1) THEN
       SLMASKL(I) = 1.0_KIND_io8
       SLMASKW(I) = 1.0_KIND_io8
     ELSE
       SLMASKL(I) = 0.0_KIND_io8
       SLMASKW(I) = 0.0_KIND_io8
     ENDIF
     if(nint(slmask(i)) == 0) then
       min_ice(i) = 0.15_KIND_io8
     else
       min_ice(i) = 0.0_KIND_io8
     endif
   ENDDO  
   num_threads = num_parthds()
   PRINT*
   PRINT*,"CALL SFCCYCLE TO UPDATE SURFACE FIELDS."
   CALL SFCCYCLE(LUGB,LENSFC,LSOIL,SIG1T,DELTSFC,          &
               IY,IM,ID,IH,FH,RLA,RLO,                   &
               SLMASKL,SLMASKW, OROG, OROG_UF, USE_UFO, DO_NSST,   &
               SIHFCS,SICFCS,SITFCS,SNDFCS,SLCFCS,       &
               VMNFCS,VMXFCS,SLPFCS,ABSFCS,              &
               TSFFCS,SWEFCS,ZORFCS,ALBFCS,TG3FCS,       &
               CNPFCS,SMCFCS,STCFCS,SLIFCS,AISFCS,       &
               VEGFCS,VETFCS,SOTFCS,ALFFCS,              &
               CVFCS,CVBFCS,CVTFCS,MYRANK,num_threads, NLUNIT,        &
               SZ_NML, INPUT_NML_FILE,                   &
               min_ice, &
               IALB,ISOT,IVEGSRC,TILE_NUM,I_INDEX,J_INDEX)
   DEALLOCATE(SLMASKL, SLMASKW)
 ENDIF

!--------------------------------------------------------------------------------
! IF RUNNING WITH NSST, READ IN GSI FILE WITH THE UPDATED INCREMENTS (ON THE
! GAUSSIAN GRID), INTERPOLATE INCREMENTS TO THE CUBED-SPHERE TILE, AND PERFORM
! REQUIRED ADJUSTMENTS AND QC.
!--------------------------------------------------------------------------------

 IF (DO_NSST) THEN
   IF (NST_FILE == "NULL") THEN
     PRINT*
     PRINT*,"NO GSI FILE.  ADJUST IFD FOR FORMER ICE POINTS."
     DO I = 1, LENSFC
       IF (NINT(SLIFCS_FG(I)) == 2 .AND. NINT(SLIFCS(I)) == 0) THEN
         NSST%IFD(I) = 3.0
       ENDIF
     ENDDO
     NSST%TFINC = 0.0
   ELSE
     PRINT*
     PRINT*,"ADJUST TREF FROM GSI INCREMENT"
!
!    Get tf climatology at the time
!
     call get_tf_clm(rla,rlo,jdim,idim,iy,im,id,ih,tf_clm,tf_trd)
     tf_clm_tile(:) = reshape(tf_clm, (/lensfc/) )
     tf_trd_tile(:) = reshape(tf_trd, (/lensfc/) )
!
!    Get salinity climatology at the time
!
     call get_sal_clm(rla,rlo,jdim,idim,iy,im,id,ih,sal_clm)
     sal_clm_tile(:) = reshape(sal_clm, (/lensfc/) )
!
!    read tf analysis increment generated by GSI
!
     CALL READ_GSI_DATA(NST_FILE, 'NST')
!
!    update foundation & surface temperature for NSST
!
     CALL ADJUST_NSST(RLA,RLO,SLIFCS,SLIFCS_FG,TSFFCS,SITFCS,SICFCS,STCFCS, &
                    NSST,LENSFC,LSOIL,IDIM,JDIM,ZSEA1,ZSEA2,IM,ID,DELTSFC,  &
                    tf_clm_tile,tf_trd_tile,sal_clm_tile)
   ENDIF
 ENDIF

!--------------------------------------------------------------------------------
! READ IN AND APPLY LAND INCREMENTS FROM THE GSI
!--------------------------------------------------------------------------------

 IF (DO_LNDINC) THEN

    ! SNOW INCREMENTS
    ! do snow first, as temperature updates will use snow analaysis
    IF (DO_SNO_INC) THEN
    ! updates are made to snow depth only over land (and not-land ice).
    ! SWE is then updated from the snow depth analysis, using the model
    ! forecast density

       !--------------------------------------------------------------------------------
       ! read increments in
       !--------------------------------------------------------------------------------

       ! Only coded for DA on native model grid (would always be the case for cycling DA)
       CALL READ_DATA(LSOIL,LENSFC,.false.,.true.,SNDFCS=SND_INC)

       !--------------------------------------------------------------------------------
       ! add increments to state vars
       !--------------------------------------------------------------------------------

       ! store background states
       SND_BCK = SNDFCS
       SWE_BCK = SWEFCS

       CALL ADD_INCREMENT_SNOW(SND_INC,LANDINC_MASK,LENSFC,SNDFCS)

       !--------------------------------------------------------------------------------
       ! make any necessary adjustments to dependent variables
       !--------------------------------------------------------------------------------

       CALL APPLY_LAND_DA_ADJUSTMENTS_SND(LSM, LENSFC, LANDINC_MASK, SWE_BCK, SND_BCK, &
                        SNDFCS, SWEFCS)

    ENDIF

    ! SOIL INCREMENTS
    IF (DO_SOI_INC) THEN

       !--------------------------------------------------------------------------------
       ! re-calculate soilsnow mask if snow has been updated.
       !--------------------------------------------------------------------------------

        LANDINC_MASK_FG = LANDINC_MASK

        IF (DO_SFCCYCLE .OR. DO_SNO_INC)  THEN
            CALL CALCULATE_LANDINC_MASK(SMCFCS(:,1),SWEFCS, VETFCS, LENSFC, &
                                                        VEG_TYPE_LANDICE, LANDINC_MASK )
        ENDIF

       !--------------------------------------------------------------------------------
       ! read increments in
       !--------------------------------------------------------------------------------

        WRITE(RANKCH, '(I3.3)') (MYRANK+1)

        LND_SOI_FILE = trim(LND_SOI_FILE) // "." //  RANKCH

        INQUIRE(FILE=trim(LND_SOI_FILE), EXIST=file_exists)
        IF (.not. file_exists) then
            print *, 'FATAL ERROR: land increment update requested, but file does not exist: ', &
                    trim(lnd_soi_file)
            call MPI_ABORT(MPI_COMM_WORLD, 10, IERR)
        ENDIF

        CALL READ_GSI_DATA(LND_SOI_FILE, 'LND', LSOIL=LSOIL)

        !--------------------------------------------------------------------------------
        ! add increments to state vars
        !--------------------------------------------------------------------------------
        ! when applying increments, will often need to adjust other land states in response
        ! to the changes made. Need to store bacground, apply the increments, then make
        ! secondart adjustments. When updating more than one state, be careful about the
        ! order if increments and secondary adjustments.

        ! store background states
        STC_BCK = STCFCS
        SMC_BCK = SMCFCS
        SLC_BCK = SLCFCS

        ! below updates [STC/SMC/STC]FCS to hold the analysis
        CALL ADD_INCREMENT_SOIL(RLA,RLO,STCFCS,SMCFCS,SLCFCS,STC_UPDATED, SLC_UPDATED, &
                LANDINC_MASK,LANDINC_MASK_FG,LENSFC,LSOIL,IDIM,JDIM,LSM,MYRANK)

        !--------------------------------------------------------------------------------
        ! make any necessary adjustments to dependent variables
        !--------------------------------------------------------------------------------


        CALL APPLY_LAND_DA_ADJUSTMENTS_SOIL(LSM, ISOT, IVEGSRC,LENSFC, LSOIL, &
            SOTFCS, LANDINC_MASK_FG, STC_BCK, STCFCS, SMCFCS, SLCFCS, STC_UPDATED, &
            SLC_UPDATED,ZSOIL)

   ENDIF ! soil increments

!--------------------------------------------------------------------------------
! clean up
!--------------------------------------------------------------------------------

   ! to do - save and write out  STC_INC? (soil temperature increments)
   IF(ALLOCATED(LANDINC_MASK_FG)) DEALLOCATE(LANDINC_MASK_FG)
   IF(ALLOCATED(LANDINC_MASK)) DEALLOCATE(LANDINC_MASK)
   IF(ALLOCATED(STC_BCK)) DEALLOCATE(STC_BCK)
   IF(ALLOCATED(SMC_BCK)) DEALLOCATE(SMC_BCK)
   IF(ALLOCATED(SLC_BCK)) DEALLOCATE(SLC_BCK)
   IF(ALLOCATED(SND_BCK)) DEALLOCATE(SND_BCK)
   IF(ALLOCATED(SWE_BCK)) DEALLOCATE(SWE_BCK)
   IF(ALLOCATED(SND_INC)) DEALLOCATE(SND_INC)

 ENDIF
!--------------------------------------------------------------------------------
! WRITE OUT UPDATED SURFACE DATA ON THE CUBED-SPHERE TILE.
!--------------------------------------------------------------------------------

 IF (LSM==LSM_NOAHMP) THEN

   CALL WRITE_DATA(LENSFC,IDIM,JDIM,LSOIL,DO_NSST,NSST,VEGFCS=VEGFCS, &
                   SLCFCS=SLCFCS,SMCFCS=SMCFCS,STCFCS=STCFCS)

 ELSEIF (LSM==LSM_NOAH) THEN

   CALL WRITE_DATA(LENSFC,IDIM,JDIM,LSOIL, &
                   DO_NSST,NSST,SLIFCS=SLIFCS,TSFFCS=TSFFCS,VEGFCS=VEGFCS, &
                   SWEFCS=SWEFCS,TG3FCS=TG3FCS,ZORFCS=ZORFCS, &
                   ALBFCS=ALBFCS,ALFFCS=ALFFCS,CNPFCS=CNPFCS, &
                   F10M=F10M,T2M=T2M,Q2M=Q2M,VETFCS=VETFCS, &
                   SOTFCS=SOTFCS,USTAR=USTAR,FMM=FMM,FHH=FHH, &
                   SICFCS=SICFCS,SIHFCS=SIHFCS,SITFCS=SITFCS,TPRCP=TPRCP, &
                   SRFLAG=SRFLAG,SWDFCS=SNDFCS,VMNFCS=VMNFCS, &
                   VMXFCS=VMXFCS,SLPFCS=SLPFCS,ABSFCS=ABSFCS, &
                   SLCFCS=SLCFCS,SMCFCS=SMCFCS,STCFCS=STCFCS)

 ENDIF

 IF (DO_NSST) THEN
   DEALLOCATE(NSST%C_0)
   DEALLOCATE(NSST%C_D)
   DEALLOCATE(NSST%D_CONV)
   DEALLOCATE(NSST%DT_COOL)
   DEALLOCATE(NSST%IFD)
   DEALLOCATE(NSST%QRAIN)
   DEALLOCATE(NSST%TREF)
   DEALLOCATE(NSST%TFINC)
   DEALLOCATE(NSST%W_0)
   DEALLOCATE(NSST%W_D)
   DEALLOCATE(NSST%XS)
   DEALLOCATE(NSST%XT)
   DEALLOCATE(NSST%XTTS)
   DEALLOCATE(NSST%XU)
   DEALLOCATE(NSST%XV)
   DEALLOCATE(NSST%XZ)
   DEALLOCATE(NSST%XZTS)
   DEALLOCATE(NSST%Z_C)
   DEALLOCATE(NSST%ZM)
   DEALLOCATE(SLIFCS_FG)
 ENDIF

 RETURN

 END SUBROUTINE SFCDRV
 
 !> Read in gsi file with the updated reference temperature increments (on the gaussian
 !! grid), interpolate increments to the cubed-sphere tile, and
 !! perform required nsst adjustments and qc.
 !!
 !! @param[inout] RLA Latitude on the cubed-sphere tile
 !! @param[inout] RLO Longitude on the cubed-sphere tile
 !! @param[in] SLMSK_TILE Land-sea mask on the cubed-sphere tile
 !! @param[in] SLMSK_FG_TILE First guess land-sea mask on the cubed-sphere tile
 !! @param[inout] SKINT_TILE Skin temperature on the cubed-sphere tile
 !! @param[inout] SICET_TILE Ice temperature on the cubed-sphere tile
 !! @param[inout] sice_tile Ice concentration on the cubed-sphere tile
 !! @param[inout] SOILT_TILE Soil temperature on the cubed-sphere tile
 !! @param[in] NSST Data structure holding nsst fields
 !! @param[in] LENSFC Number of points on a tile
 !! @param[in] LSOIL Number of soil layers
 !! @param[in] IDIM 'I' dimension of a tile
 !! @param[in] JDIM 'J' dimension of a tile
 !! @param[in] ZSEA1 When running nsst model, this is the lower bound of
 !! depth of sea temperature. In whole mm.
 !! @param[in] ZSEA2 When running nsst model, this is the upper bound of
 !! depth of sea temperature. In whole mm.
 !! @param[in] MON Month
 !! @param[in] DAY Day
 !! @param[in] DELTSFC Cycling frequency in hours
 !! @param[in] tf_clm_tile Climatological reference temperature on the
 !! cubed-sphere tile.
 !! @param[in] tf_trd_tile Climatolocial reference temperature trend on the
 !! cubed-sphere tile.
 !! @param[in] sal_clm_tile Climatological salinity on the cubed-sphere tile.
 !!
 !! @author Xu Li, George Gayno
 SUBROUTINE ADJUST_NSST(RLA,RLO,SLMSK_TILE,SLMSK_FG_TILE,SKINT_TILE,&
                        SICET_TILE,sice_tile,SOILT_TILE,NSST,LENSFC,LSOIL,    &
                        IDIM,JDIM,ZSEA1,ZSEA2,MON,DAY,DELTSFC, &
                        tf_clm_tile,tf_trd_tile,sal_clm_tile)

 USE UTILS
 USE GDSWZD_MOD
 USE READ_WRITE_DATA, ONLY : IDIM_GAUS, JDIM_GAUS, &
                             SLMSK_GAUS, DTREF_GAUS, &
                             NSST_DATA

 USE MPI

 IMPLICIT NONE

 INTEGER, INTENT(IN)      :: LENSFC, LSOIL, IDIM, JDIM, MON, DAY

 REAL, INTENT(IN)         :: SLMSK_TILE(LENSFC), SLMSK_FG_TILE(LENSFC)
 real, intent(in)         :: tf_clm_tile(lensfc),tf_trd_tile(lensfc),sal_clm_tile(lensfc)
 REAL, INTENT(IN)         :: ZSEA1, ZSEA2, DELTSFC
 REAL, INTENT(INOUT)      :: RLA(LENSFC), RLO(LENSFC), SKINT_TILE(LENSFC)
 REAL, INTENT(INOUT)      :: SICET_TILE(LENSFC),sice_tile(lensfc),SOILT_TILE(LENSFC,LSOIL)

 TYPE(NSST_DATA)          :: NSST

 REAL, PARAMETER          :: TMAX=313.0,tzero=273.16

 INTEGER                  :: IOPT, NRET, KGDS_GAUS(200)
 INTEGER                  :: IGAUS, JGAUS, IJ, II, JJ, III, JJJ, KRAD
 INTEGER                  :: ISTART, IEND, JSTART, JEND
 INTEGER                  :: MASK_TILE, MASK_FG_TILE
 INTEGER                  :: ITILE, JTILE
 INTEGER                  :: MAX_SEARCH, J, IERR
 INTEGER                  :: IGAUSP1, JGAUSP1
 integer                  :: nintp,nsearched,nice,nland
 integer                  :: nfill,nfill_tice,nfill_clm
 integer                  :: nset_thaw,nset_thaw_s,nset_thaw_i,nset_thaw_c

 INTEGER, ALLOCATABLE     :: ID1(:,:), ID2(:,:), JDC(:,:)

 LOGICAL                  :: IS_ICE
 
 real                     :: tfreez
 REAL                     :: TREF_SAVE,WSUM,tf_ice,tf_thaw
 REAL                     :: FILL, DTZM, GAUS_RES_KM, DTREF
 REAL, ALLOCATABLE        :: XPTS(:), YPTS(:), LATS(:), LONS(:)
 REAL, ALLOCATABLE        :: DUM2D(:,:), LATS_RAD(:), LONS_RAD(:)
 REAL, ALLOCATABLE        :: AGRID(:,:,:), S2C(:,:,:)

 KGDS_GAUS     = 0
 KGDS_GAUS(1)  = 4          ! OCT 6 - TYPE OF GRID (GAUSSIAN)
 KGDS_GAUS(2)  = IDIM_GAUS  ! OCT 7-8 - # PTS ON LATITUDE CIRCLE
 KGDS_GAUS(3)  = JDIM_GAUS
 KGDS_GAUS(4)  = 90000      ! OCT 11-13 - LAT OF ORIGIN
 KGDS_GAUS(5)  = 0          ! OCT 14-16 - LON OF ORIGIN
 KGDS_GAUS(6)  = 128        ! OCT 17 - RESOLUTION FLAG
 KGDS_GAUS(7)  = -90000     ! OCT 18-20 - LAT OF EXTREME POINT
 KGDS_GAUS(8)  = NINT(-360000./FLOAT(IDIM_GAUS))  ! OCT 21-23 - LON OF EXTREME POINT
 KGDS_GAUS(9)  = NINT((360.0 / FLOAT(IDIM_GAUS))*1000.0)
                            ! OCT 24-25 - LONGITUDE DIRECTION INCR.
 KGDS_GAUS(10) = JDIM_GAUS/2     ! OCT 26-27 - NUMBER OF CIRCLES POLE TO EQUATOR
 KGDS_GAUS(12) = 255        ! OCT 29 - RESERVED
 KGDS_GAUS(20) = 255        ! OCT 5  - NOT USED, SET TO 255

 PRINT*
 PRINT*,'ADJUST NSST USING GSI INCREMENTS ON GAUSSIAN GRID'

!----------------------------------------------------------------------
! CALL GDSWZD TO COMPUTE THE LAT/LON OF EACH GSI GAUSSIAN GRID POINT.
!----------------------------------------------------------------------

 IOPT = 0
 FILL = -9999.
 ALLOCATE(XPTS(IDIM_GAUS*JDIM_GAUS))
 ALLOCATE(YPTS(IDIM_GAUS*JDIM_GAUS))
 ALLOCATE(LATS(IDIM_GAUS*JDIM_GAUS))
 ALLOCATE(LONS(IDIM_GAUS*JDIM_GAUS))
 XPTS = FILL
 YPTS = FILL
 LATS = FILL
 LONS = FILL

 CALL GDSWZD(KGDS_GAUS,IOPT,(IDIM_GAUS*JDIM_GAUS),FILL,XPTS,YPTS,LONS,LATS,NRET)

 IF (NRET /= (IDIM_GAUS*JDIM_GAUS)) THEN
   PRINT*,'FATAL ERROR: PROBLEM IN GDSWZD. STOP.'
   CALL MPI_ABORT(MPI_COMM_WORLD, 12, IERR)
 ENDIF

 DEALLOCATE (XPTS, YPTS)

 ALLOCATE(DUM2D(IDIM_GAUS,JDIM_GAUS))
 DUM2D = RESHAPE(LATS, (/IDIM_GAUS,JDIM_GAUS/) )
 DEALLOCATE(LATS)

 ALLOCATE(LATS_RAD(JDIM_GAUS))
 DO J = 1, JDIM_GAUS
   LATS_RAD(J) = DUM2D(1,JDIM_GAUS-J+1) * 3.1415926 / 180.0
 ENDDO

 DUM2D = RESHAPE(LONS, (/IDIM_GAUS,JDIM_GAUS/) )
 DEALLOCATE(LONS)
 ALLOCATE(LONS_RAD(IDIM_GAUS))
 LONS_RAD = DUM2D(:,1) * 3.1415926 / 180.0

 DEALLOCATE(DUM2D)

 ALLOCATE(AGRID(IDIM,JDIM,2))
 AGRID(:,:,1) = RESHAPE (RLO, (/IDIM,JDIM/) )
 AGRID(:,:,2) = RESHAPE (RLA, (/IDIM,JDIM/) )
 AGRID        = AGRID * 3.1415926 / 180.0

 ALLOCATE(ID1(IDIM,JDIM))
 ALLOCATE(ID2(IDIM,JDIM))
 ALLOCATE(JDC(IDIM,JDIM))
 ALLOCATE(S2C(IDIM,JDIM,4))

!----------------------------------------------------------------------
! COMPUTE BILINEAR WEIGHTS FOR EACH MODEL POINT FROM THE NEAREST
! FOUR GSI/GAUSSIAN POINTS.  DOES NOT ACCOUNT FOR MASK.  THAT
! HAPPENS LATER.
!----------------------------------------------------------------------

 CALL REMAP_COEF( 1, IDIM, 1, JDIM, IDIM_GAUS, JDIM_GAUS, &
                  LONS_RAD, LATS_RAD, ID1, ID2, JDC, S2C, AGRID )

 DEALLOCATE(LONS_RAD, LATS_RAD, AGRID)

!----------------------------------------------------------------------
! THE MAXIMUM DISTANCE TO SEARCH IS 500 KM. HOW MANY GAUSSIAN
! GRID LENGTHS IS THAT?
!----------------------------------------------------------------------

 GAUS_RES_KM = 360.0 / IDIM_GAUS * 111.0
 MAX_SEARCH  = CEILING(500.0/GAUS_RES_KM)

 PRINT*
 PRINT*,'MAXIMUM SEARCH IS ',MAX_SEARCH, ' GAUSSIAN POINTS.'
 PRINT*

!
! Initialize variables for counts statitics to be zeros
!
 nintp = 0
 nset_thaw = 0
 nset_thaw_s = 0
 nset_thaw_i = 0
 nset_thaw_c = 0
 nsearched = 0
 nfill = 0
 nfill_tice = 0
 nfill_clm = 0
 nice = 0
 nland = 0
!----------------------------------------------------------------------
! TREF INCREMENT WILL BE OUTPUT.  INITIALIZE TO ZERO.
!----------------------------------------------------------------------

 NSST%TFINC = 0.0

 IJ_LOOP : DO IJ = 1, LENSFC

   MASK_TILE    = NINT(SLMSK_TILE(IJ))
   MASK_FG_TILE = NINT(SLMSK_FG_TILE(IJ))

!
!  when sea ice exists, get salinity dependent water temperature
!
   tf_ice = tfreez(sal_clm_tile(ij)) + tzero
!----------------------------------------------------------------------
! SKIP LAND POINTS.  NSST NOT APPLIED AT LAND.
!----------------------------------------------------------------------

   IF (MASK_TILE == 1) THEN
     nland = nland + 1
     CYCLE IJ_LOOP  
   ENDIF

!
! these are ice points.  set tref to tf_ice and update tmpsfc.
!
   if (mask_tile == 2) then
     nsst%tref(ij)=tf_ice      ! water part tmp set
     skint_tile(ij)=(1.0-sice_tile(ij))*nsst%tref(ij)+sice_tile(ij)*sicet_tile(ij)
     nice = nice + 1
     cycle ij_loop
   endif

!
!  Get i,j index on array of (idim,jdim) from known ij
!
   JTILE = (IJ-1) / IDIM + 1
   ITILE = MOD(IJ,IDIM)
   IF (ITILE==0) ITILE = IDIM

!----------------------------------------------------------------------
! IF THE MODEL POINT WAS ICE COVERED, BUT IS NOW OPEN WATER, SET
! TREF TO searched adjascent open water onea, if failed the search, set to
! weighted average of tf_ice and tf_clm. For NSST vars, set xz TO '30' AND ALL OTHER FIELDS TO ZERO.
!----------------------------------------------------------------------

   IF (MASK_FG_TILE == 2 .AND. MASK_TILE == 0) THEN
!
!    set background for the thaw (just melted water) situation
!
     call tf_thaw_set(nsst%tref,nint(slmsk_fg_tile),itile,jtile,tf_ice,tf_clm_tile(ij),tf_thaw,idim,jdim, &
                      nset_thaw_s,nset_thaw_i,nset_thaw_c)
     call nsst_water_reset(nsst,ij,tf_thaw)
     nset_thaw = nset_thaw + 1
   ENDIF

!----------------------------------------------------------------------
! THESE ARE POINTS THAT ARE OPEN WATER AND WERE OPEN WATER PRIOR
! TO ANY ICE UPDATE BY SFCCYCLE. UPDATE TREF AND SKIN TEMP.
! AT OPEN WATER POINTS, THE SEA ICE TEMPERATURE (SICET_TILE) AND
! SOIL COLUMN TEMPERATURE (SOILT_TILE) ARE SET TO THE SKIN TEMP.
! IT IS SIMPLY A FILLER VALUE.  THESE FIELDS ARE NOT USED AT
! OPEN WATER POINTS.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! SEE IF ANY OF THE NEAREST GSI POINTS MASK AREA OPEN WATER.  
! IF SO, APPLY NSST INCREMENT USING BILINEAR INTERPOLATION.
!----------------------------------------------------------------------

   IGAUS   = ID1(ITILE,JTILE)
   JGAUS   = JDC(ITILE,JTILE)
   IGAUSP1 = ID2(ITILE,JTILE)
   JGAUSP1 = JDC(ITILE,JTILE)+1

   IF (SLMSK_GAUS(IGAUS,JGAUS)     == 0 .OR. &
       SLMSK_GAUS(IGAUSP1,JGAUS)   == 0 .OR. &
       SLMSK_GAUS(IGAUSP1,JGAUSP1) == 0 .OR. &
       SLMSK_GAUS(IGAUS,JGAUSP1)   == 0) THEN

     DTREF = 0.0
     WSUM  = 0.0

     IF (SLMSK_GAUS(IGAUS,JGAUS) == 0) THEN
       DTREF = DTREF + (S2C(ITILE,JTILE,1) * DTREF_GAUS(IGAUS,JGAUS))
       WSUM  = WSUM + S2C(ITILE,JTILE,1)
     ENDIF

     IF (SLMSK_GAUS(IGAUSP1,JGAUS) == 0) THEN
       DTREF = DTREF + (S2C(ITILE,JTILE,2) * DTREF_GAUS(IGAUSP1,JGAUS))
       WSUM  = WSUM + S2C(ITILE,JTILE,2)
     ENDIF

     IF (SLMSK_GAUS(IGAUSP1,JGAUSP1) == 0) THEN
       DTREF = DTREF + (S2C(ITILE,JTILE,3) * DTREF_GAUS(IGAUSP1,JGAUSP1))
       WSUM  = WSUM + S2C(ITILE,JTILE,3)
     ENDIF

     IF (SLMSK_GAUS(IGAUS,JGAUSP1) == 0) THEN
       DTREF = DTREF + (S2C(ITILE,JTILE,4) * DTREF_GAUS(IGAUS,JGAUSP1))
       WSUM  = WSUM + S2C(ITILE,JTILE,4)
     ENDIF

     nintp = nintp + 1
     DTREF = DTREF / WSUM

     TREF_SAVE      = NSST%TREF(IJ)
     NSST%TREF(IJ)  = NSST%TREF(IJ) + DTREF
     NSST%TREF(IJ)  = MAX(NSST%TREF(IJ), tf_ice)
     NSST%TREF(IJ)  = MIN(NSST%TREF(IJ), TMAX)
     NSST%TFINC(IJ) = NSST%TREF(IJ) - TREF_SAVE

     CALL DTZM_POINT(NSST%XT(IJ),NSST%XZ(IJ),NSST%DT_COOL(IJ),  &
                     NSST%Z_C(IJ),ZSEA1,ZSEA2,DTZM)

     SKINT_TILE(IJ) = NSST%TREF(IJ) + DTZM
     SKINT_TILE(IJ) = MAX(SKINT_TILE(IJ), tf_ice)
     SKINT_TILE(IJ) = MIN(SKINT_TILE(IJ), TMAX)

     SICET_TILE(IJ)   = SKINT_TILE(IJ)
     SOILT_TILE(IJ,:) = SKINT_TILE(IJ)

!----------------------------------------------------------------------
! NO NEARBY GSI/GAUSSIAN OPEN WATER POINTS. PERFORM A SPIRAL SEARCH TO
! FIND NEAREST NON-LAND POINT ON GSI/GAUSSIAN GRID.
!----------------------------------------------------------------------

   ELSE  

     IS_ICE = .FALSE.

     DO KRAD = 1, MAX_SEARCH

       ISTART = IGAUS - KRAD
       IEND   = IGAUS + KRAD
       JSTART = JGAUS - KRAD
       JEND   = JGAUS + KRAD

       DO JJ = JSTART, JEND
       DO II = ISTART, IEND

         IF((JJ == JSTART) .OR. (JJ == JEND) .OR.   &
            (II == ISTART) .OR. (II == IEND))  THEN

           IF ((JJ >= 1) .AND. (JJ <= JDIM_GAUS)) THEN

             JJJ = JJ
             IF (II <= 0) THEN
               III = IDIM_GAUS + II
             ELSE IF (II >= (IDIM_GAUS+1)) THEN
               III = II - IDIM_GAUS
             ELSE
               III = II
             END IF

!----------------------------------------------------------------------
! SEE IF NEARBY POINTS ARE SEA ICE.  IF THEY ARE, AND THE SEARCH FOR
! A GAUSSIAN GRID OPEN WATER POINT FAILS, THEN TREF WILL BE SET TO
! FREEZING BELOW.
!----------------------------------------------------------------------

             IF (KRAD <= 2 .AND. SLMSK_GAUS(III,JJJ) == 2) IS_ICE = .TRUE.

             IF (SLMSK_GAUS(III,JJJ) == 0) THEN

!              PRINT*,'MISMATCH AT TILE POINT  ',ITILE,JTILE
!              PRINT*,'UPDATE TREF USING GSI INCREMENT AT ',III,JJJ,DTREF_GAUS(III,JJJ)
               nsearched = nsearched + 1

               TREF_SAVE      = NSST%TREF(IJ)
               NSST%TREF(IJ ) = NSST%TREF(IJ) + DTREF_GAUS(III,JJJ)
               NSST%TREF(IJ)  = MAX(NSST%TREF(IJ), tf_ice)
               NSST%TREF(IJ)  = MIN(NSST%TREF(IJ), TMAX)
               NSST%TFINC(IJ) = NSST%TREF(IJ) - TREF_SAVE

               CALL DTZM_POINT(NSST%XT(IJ),NSST%XZ(IJ),NSST%DT_COOL(IJ),  &
                     NSST%Z_C(IJ),ZSEA1,ZSEA2,DTZM)

               SKINT_TILE(IJ) = NSST%TREF(IJ) + DTZM
               SKINT_TILE(IJ) = MAX(SKINT_TILE(IJ), tf_ice)
               SKINT_TILE(IJ) = MIN(SKINT_TILE(IJ), TMAX)

               SICET_TILE(IJ)   = SKINT_TILE(IJ)
               SOILT_TILE(IJ,:) = SKINT_TILE(IJ)
               CYCLE IJ_LOOP

             ENDIF ! GSI/Gaussian mask is open water

           ENDIF

         ENDIF

       ENDDO
       ENDDO

     ENDDO ! KRAD LOOP

!----------------------------------------------------------------------
! THE SEARCH FAILED.  IF THERE IS NEARBY ICE, SET TREF TO FREEZING.
! ELSE UPDATE TREF BASED ON THE ANNUAL SST CYCLE.
!----------------------------------------------------------------------

!    PRINT*,'WARNING !!!!!! SEARCH FAILED AT TILE POINT ',ITILE,JTILE

     nfill = nfill  + 1
     IF (IS_ICE) THEN
       NSST%TREF(IJ) = tf_ice
!      PRINT*,"NEARBY ICE.  SET TREF TO FREEZING"
       nfill_tice = nfill_tice + 1
     ELSE
       TREF_SAVE      = NSST%TREF(IJ)
       NSST%TREF(IJ)  = NSST%TREF(IJ) + TF_TRD_TILE(IJ)
       NSST%TREF(IJ)  = MAX(NSST%TREF(IJ), tf_ice)
       NSST%TREF(IJ)  = MIN(NSST%TREF(IJ), TMAX)
       NSST%TFINC(IJ) = NSST%TREF(IJ) - TREF_SAVE
!      PRINT*,'UPDATE TREF FROM SST CLIMO ',DTREF
       nfill_clm = nfill_clm + 1
     ENDIF

     CALL DTZM_POINT(NSST%XT(IJ),NSST%XZ(IJ),NSST%DT_COOL(IJ),  &
                     NSST%Z_C(IJ),ZSEA1,ZSEA2,DTZM)

     SKINT_TILE(IJ) = NSST%TREF(IJ) + DTZM
     SKINT_TILE(IJ) = MAX(SKINT_TILE(IJ), tf_ice)
     SKINT_TILE(IJ) = MIN(SKINT_TILE(IJ), TMAX)

     SICET_TILE(IJ)   = SKINT_TILE(IJ)
     SOILT_TILE(IJ,:) = SKINT_TILE(IJ)

   ENDIF  ! NEARBY GAUSSIAN POINTS ARE OPEN WATER?

 ENDDO IJ_LOOP

 write(*,'(a)') 'statistics of grids number processed for tile : '
 write(*,'(a,I8)') ' nintp = ',nintp
 write(*,'(a,4I8)') 'nset_thaw,nset_thaw_s,nset_thaw_i,nset_thaw_c =',nset_thaw,nset_thaw_s,nset_thaw_i,nset_thaw_c
 write(*,'(a,I8)') ' nsearched = ',nsearched
 write(*,'(a,3I6)') ' nfill,nfill_tice,nfill_clm = ',nfill,nfill_tice,nfill_clm
 write(*,'(a,I8)') ' nice = ',nice
 write(*,'(a,I8)') ' nland = ',nland

 DEALLOCATE(ID1, ID2, JDC, S2C)

 END SUBROUTINE ADJUST_NSST

 !> If the tile point is an isolated water point that has no
 !! corresponding gsi water point, then tref is updated using the rtg
 !! sst climo trend. This monthly trend is sorted by latitude band.
 !!
 !! @param[in] LATITUDE Latitude of tile point
 !! @param[in] MON Month
 !! @param[in] DAY Day
 !! @param[in] DELTSFC Cycling frequency in hours
 !! @param[out] DTREF Monthly trend of reference temperature
 !! @author Xu Li, George Gayno
 SUBROUTINE CLIMO_TREND(LATITUDE, MON, DAY, DELTSFC, DTREF)
 IMPLICIT NONE

 INTEGER, INTENT(IN)    :: MON, DAY

 REAL, INTENT(IN)       :: LATITUDE, DELTSFC
 REAL, INTENT(OUT)      :: DTREF

 INTEGER                :: NUM_DAYS(12), MON2, MON1

 REAL, TARGET           :: SST_80_90(12)
 REAL, TARGET           :: SST_70_80(12)
 REAL, TARGET           :: SST_60_70(12)
 REAL, TARGET           :: SST_50_60(12)
 REAL, TARGET           :: SST_40_50(12)
 REAL, TARGET           :: SST_30_40(12)
 REAL, TARGET           :: SST_20_30(12)
 REAL, TARGET           :: SST_10_20(12)
 REAL, TARGET           :: SST_00_10(12)
 REAL, TARGET           :: SST_M10_00(12)
 REAL, TARGET           :: SST_M20_M10(12)
 REAL, TARGET           :: SST_M30_M20(12)
 REAL, TARGET           :: SST_M40_M30(12)
 REAL, TARGET           :: SST_M50_M40(12)
 REAL, TARGET           :: SST_M60_M50(12)
 REAL, TARGET           :: SST_M70_M60(12)
 REAL, TARGET           :: SST_M80_M70(12)
 REAL, TARGET           :: SST_M90_M80(12)

 REAL, POINTER          :: SST(:)

 DATA NUM_DAYS  /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

 DATA SST_80_90 /271.466, 271.458, 271.448, 271.445, 271.519, 271.636, & 
                 272.023, 272.066, 272.001, 271.698, 271.510, 271.472/

 DATA SST_70_80 /272.149, 272.103, 272.095, 272.126, 272.360, 272.988, &
                 274.061, 274.868, 274.415, 273.201, 272.468, 272.268/

 DATA SST_60_70 /274.240, 274.019, 273.988, 274.185, 275.104, 276.875, &
                 279.005, 280.172, 279.396, 277.586, 275.818, 274.803/

 DATA SST_50_60 /277.277, 276.935, 277.021, 277.531, 279.100, 281.357, &
                 283.735, 285.171, 284.399, 282.328, 279.918, 278.199/

 DATA SST_40_50 /281.321, 280.721, 280.850, 281.820, 283.958, 286.588, &
                 289.195, 290.873, 290.014, 287.652, 284.898, 282.735/

 DATA SST_30_40 /289.189, 288.519, 288.687, 289.648, 291.547, 293.904, &
                 296.110, 297.319, 296.816, 295.225, 292.908, 290.743/

 DATA SST_20_30 /294.807, 294.348, 294.710, 295.714, 297.224, 298.703, &
                 299.682, 300.127, 300.099, 299.455, 297.953, 296.177/

 DATA SST_10_20 /298.878, 298.720, 299.033, 299.707, 300.431, 300.709, &
                 300.814, 300.976, 301.174, 301.145, 300.587, 299.694/

 DATA SST_00_10 /300.415, 300.548, 300.939, 301.365, 301.505, 301.141, &
                 300.779, 300.660, 300.818, 300.994, 300.941, 300.675/

 DATA SST_M10_00 /300.226, 300.558, 300.914, 301.047, 300.645, 299.870, &
                  299.114, 298.751, 298.875, 299.294, 299.721, 299.989/

 DATA SST_M20_M10 /299.547, 299.985, 300.056, 299.676, 298.841, 297.788, &
                   296.893, 296.491, 296.687, 297.355, 298.220, 298.964/

 DATA SST_M30_M20 /297.524, 298.073, 297.897, 297.088, 295.846, 294.520, &
                   293.525, 293.087, 293.217, 293.951, 295.047, 296.363/

 DATA SST_M40_M30 /293.054, 293.765, 293.468, 292.447, 291.128, 289.781, &
                   288.773, 288.239, 288.203, 288.794, 289.947, 291.553/

 DATA SST_M50_M40 /285.052, 285.599, 285.426, 284.681, 283.761, 282.826, &
                   282.138, 281.730, 281.659, 281.965, 282.768, 283.961/

 DATA SST_M60_M50 /277.818, 278.174, 277.991, 277.455, 276.824, 276.229, &
                   275.817, 275.585, 275.560, 275.687, 276.142, 276.968/

 DATA SST_M70_M60 /273.436, 273.793, 273.451, 272.813, 272.349, 272.048, &
                   271.901, 271.838, 271.845, 271.889, 272.080, 272.607/

 DATA SST_M80_M70 /271.579, 271.578, 271.471, 271.407, 271.392, 271.391, &
                   271.390, 271.391, 271.394, 271.401, 271.422, 271.486/

 DATA SST_M90_M80 /271.350, 271.350, 271.350, 271.350, 271.350, 271.350, &
                   271.350, 271.350, 271.350, 271.350, 271.350, 271.350/

 NULLIFY(SST)
 IF (LATITUDE > 80.0) THEN
   SST => SST_80_90
 ELSEIF (LATITUDE > 70.0) THEN
   SST => SST_70_80
 ELSEIF (LATITUDE > 60.0) THEN
   SST => SST_60_70
 ELSEIF (LATITUDE > 50.0) THEN
   SST => SST_50_60
 ELSEIF (LATITUDE > 40.0) THEN
   SST => SST_40_50
 ELSEIF (LATITUDE > 30.0) THEN
   SST => SST_30_40
 ELSEIF (LATITUDE > 20.0) THEN
   SST => SST_20_30
 ELSEIF (LATITUDE > 10.0) THEN
   SST => SST_10_20
 ELSEIF (LATITUDE > 0.0) THEN
   SST => SST_00_10
 ELSEIF (LATITUDE > -10.0) THEN
   SST => SST_M10_00
 ELSEIF (LATITUDE > -20.0) THEN
   SST => SST_M20_M10
 ELSEIF (LATITUDE > -30.0) THEN
   SST => SST_M30_M20
 ELSEIF (LATITUDE > -40.0) THEN
   SST => SST_M40_M30
 ELSEIF (LATITUDE > -50.0) THEN
   SST => SST_M50_M40
 ELSEIF (LATITUDE > -60.0) THEN
   SST => SST_M60_M50
 ELSEIF (LATITUDE > -70.0) THEN
   SST => SST_M70_M60
 ELSEIF (LATITUDE > -80.0) THEN
   SST => SST_M80_M70
 ELSE
   SST => SST_M90_M80
 END IF

 IF (DAY >= 15) THEN
   MON2 = MON+1
   IF(MON2 == 13) MON2 = 1
   MON1 = MON
   DTREF = (SST(MON2) - SST(MON1)) / NUM_DAYS(MON1)
 ELSE
   MON1 = MON - 1
   IF (MON1 == 0) MON1=12
   MON2 = MON
   DTREF = (SST(MON2) - SST(MON1)) / NUM_DAYS(MON1)
 ENDIF

 DTREF = DTREF * (DELTSFC / 24.0)

 END SUBROUTINE CLIMO_TREND

 !> Compute the vertical mean of the NSST t-profile.
 !!                                                                 
 !! @param[in] xt Heat content in the diurnal thermocline layer.
 !! @param[in] xz Thickness of the diurnal thermocline layer.
 !! @param[in] dt_cool Skin-layer cooling amount.
 !! @param[in] zc Thickness of skin-layer.
 !! @param[in] z1 Lower bound of depth of sea temperature.
 !! @param[in] z2 Upper bound of depth of sea temperature.
 !! @param[out] dtzm Mean of the NSST t-profile from z1 to z2.
 !!
 !! @author Xu Li @date 2015
 SUBROUTINE DTZM_POINT(XT,XZ,DT_COOL,ZC,Z1,Z2,DTZM)
  implicit none

  real, intent(in)  :: xt,xz,dt_cool,zc,z1,z2
  real, intent(out) :: dtzm

  real, parameter   :: zero = 0.0
  real, parameter   :: one  = 1.0
  real, parameter   :: half = 0.5
  real              :: dt_warm,dtw,dtc
!
! get the mean warming in the range of z=z1 to z=z2
!
  dtw = zero
  if ( xt > zero ) then
    dt_warm = (xt+xt)/xz      ! Tw(0)
    if ( z1 < z2) then
      if ( z2 < xz ) then
        dtw = dt_warm*(one-(z1+z2)/(xz+xz))
      elseif ( z1 < xz .and. z2 >= xz ) then
        dtw = half*(one-z1/xz)*dt_warm*(xz-z1)/(z2-z1)
      endif
    elseif ( z1 == z2 ) then
      if ( z1 < xz ) then
        dtw = dt_warm*(one-z1/xz)
      endif
    endif
  endif
!
! get the mean cooling in the range of z=z1 to z=z2
!
  dtc = zero
  if ( zc > zero ) then
    if ( z1 < z2) then
      if ( z2 < zc ) then
        dtc = dt_cool*(one-(z1+z2)/(zc+zc))
      elseif ( z1 < zc .and. z2 >= zc ) then
        dtc = half*(one-z1/zc)*dt_cool*(zc-z1)/(z2-z1)
      endif
    elseif ( z1 == z2 ) then
      if ( z1 < zc ) then
        dtc = dt_cool*(one-z1/zc)
      endif
    endif
  endif

!
! get the mean T departure from Tf in the range of z=z1 to z=z2
!
  dtzm = dtw - dtc

 END SUBROUTINE DTZM_POINT

 !> Set the background reference temperature (tf) for points where
 !! the ice has just melted.
 !!
 !! @param[in] tf_ij Foundation temperature background on FV3 native grids.
 !! @param[in] mask_ij Mask of the tile (FV3 native grids).
 !! @param[in] itile Location index in the 'i' direction.
 !! @param[in] jtile Location index in the 'j' direction.
 !! @param[in] tice Water temperature (calulated with a salinity formula).
 !! @param[in] tclm SST climatology valid at the analysis time.
 !! @param [out] tf_thaw Foundation temperature of thawed points.
 !! @param[inout] nx 'i' dimension of tf_ij
 !! @param[inout] ny 'j' dimension of tf_ij
 !! @param[inout] nset_thaw_s Number of foundation temperature points filled
 !! via a search.
 !! @param[inout] nset_thaw_i Number of ice points filled with a calculated
 !! tice.
 !! @param[inout] nset_thaw_c Number of points filled with a weighted 
 !! average of tice and tclm.
 !! @author Xu Li
 subroutine tf_thaw_set(tf_ij,mask_ij,itile,jtile,tice,tclm,tf_thaw,nx,ny, &
                        nset_thaw_s,nset_thaw_i,nset_thaw_c)

 real,    dimension(nx*ny), intent(in)    :: tf_ij
 integer, dimension(nx*ny), intent(in)    :: mask_ij
 real,                      intent(in)    :: tice,tclm
 integer,                   intent(in)    :: itile,jtile,nx,ny
 real,                      intent(out)   :: tf_thaw
 integer,                   intent(inout) :: nset_thaw_s,nset_thaw_i,nset_thaw_c
! Local
 real, parameter :: bmiss = -999.0
 real,    dimension(nx,ny) :: tf
 integer, dimension(nx,ny) :: mask
 integer :: krad,max_search,istart,iend,jstart,jend
 integer :: ii,jj,iii,jjj
 logical :: is_ice

 max_search = 2

 mask(:,:) = reshape(mask_ij,(/nx,ny/) )
 tf(:,:)   = reshape(tf_ij,(/nx,ny/) )

 tf_thaw = bmiss

 do krad = 1, max_search

    istart = itile - krad
    iend   = itile + krad
    jstart = jtile - krad
    jend   = jtile + krad

    do jj = jstart, jend
       do ii = istart, iend


          if ((jj == jstart) .or. (jj == jend) .or.   &
              (ii == istart) .or. (ii == iend))  then

             if ((jj >= 1) .and. (jj <= ny)) then
                jjj = jj
                if (ii <= 0) then
                   iii = nx + ii
                else if (ii >= (nx+1)) then
                   iii = ii - nx
                else
                   iii = ii
                endif

!----------------------------------------------------------------------
! SEE IF NEARBY POINTS ARE SEA ICE.  IF THEY ARE, AND THE SEARCH FOR
! A GAUSSIAN GRID OPEN WATER POINT FAILS, THEN TREF WILL BE SET TO
! FREEZING BELOW.
!----------------------------------------------------------------------
                if (krad <= 2 .and. mask(iii,jjj) == 2) is_ice = .true.

                if (mask(iii,jjj) == 0) then
                   tf_thaw = tf(iii,jjj)
                   nset_thaw_s = nset_thaw_s + 1
                   write(*,'(a,I4,2F9.3)') 'nset_thaw_s,tf(iii,jjj),tclm : ',nset_thaw_s,tf(iii,jjj),tclm
                   go to 100
                endif ! tile mask is open water

             endif
          endif
       enddo
    enddo
 enddo  ! krad loop

 100 continue

 if ( tf_thaw == bmiss ) then
    if (is_ice) then
       tf_thaw = tice
       nset_thaw_i = nset_thaw_i + 1
       write(*,'(a,I4,F9.3)') 'nset_thaw_i,tf_ice : ',nset_thaw_i,tice
    else
       tf_thaw = 0.8*tice+0.2*tclm
       nset_thaw_c = nset_thaw_c + 1
       write(*,'(a,I4,2F9.3)') 'nset_thaw_c,tf_ice,tclm : ',nset_thaw_c,tice,tclm
    endif
 endif

 end subroutine tf_thaw_set

 !> If the first guess was sea ice, but the analysis is open water,
 !! reset all nsst variables.
 !! 
 !! @param[inout] nsst Data structure that holds the nsst fields
 !! @param[in] ij Index of point to be updated
 !! @param[in] tf_thaw Reference temperature for former ice points
 !! @author Xu Li
 subroutine nsst_water_reset(nsst,ij,tf_thaw)
 use read_write_data, only : nsst_data
 implicit none

 integer, intent(in)  :: ij

 real, intent(in)     :: tf_thaw

 type(nsst_data), intent(inout)      :: nsst

 nsst%c_0(ij)     = 0.0
 nsst%c_d(ij)     = 0.0
 nsst%d_conv(ij)  = 0.0
 nsst%dt_cool(ij) = 0.0
 nsst%ifd(ij)     = 0.0
 nsst%qrain(ij)   = 0.0
 nsst%tref(ij)    = tf_thaw
 nsst%w_0(ij)     = 0.0
 nsst%w_d(ij)     = 0.0
 nsst%xs(ij)      = 0.0
 nsst%xt(ij)      = 0.0
 nsst%xtts(ij)    = 0.0
 nsst%xu(ij)      = 0.0
 nsst%xv(ij)      = 0.0
 nsst%xz(ij)      = 30.0
 nsst%xzts(ij)    = 0.0
 nsst%z_c(ij)     = 0.0
 nsst%zm(ij)      = 0.0

 end subroutine nsst_water_reset

 !> Get the sst climatology at the valid time and on the target
 !! grid.
 !!
 !! @param[in] xlats_ij latitude of target grid
 !! @param[in] xlons_ij longitude of target grid
 !! @param[in] ny 'j' dimension of target grid
 !! @param[in] nx 'i' dimension of target grid
 !! @param[in] iy Year 
 !! @param[in] im Month
 !! @param[in] id Day
 !! @param[in] ih Hour
 !! @param[out] tf_clm sst climatology at the valid time and on the target grid
 !! @param[out] tf_trd 6-hourly sst climatology tendency at the valid time
 !! and on the target grid.
 !! @author Xu Li
subroutine get_tf_clm(xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,tf_clm,tf_trd)
 use read_write_data, only : get_tf_clm_dim

 implicit none

 real,    dimension(nx*ny), intent(in)  :: xlats_ij 
 real,    dimension(nx*ny), intent(in)  :: xlons_ij 
 real,    dimension(nx,ny), intent(out) :: tf_clm   
 real,    dimension(nx,ny), intent(out) :: tf_trd   
 integer, intent(in) :: iy,im,id,ih,nx,ny
! local declare
 real,    allocatable, dimension(:,:)   :: tf_clm0    ! sst climatology at the valid time (nxc,nyc)
 real,    allocatable, dimension(:,:)   :: tf_trd0    ! 6-hourly sst climatology tendency valid at atime (nxc,nyc)
 real,    allocatable, dimension(:)     :: cxlats     ! latitudes of sst climatology
 real,    allocatable, dimension(:)     :: cxlons     ! longitudes of sst climatology

 real,    dimension(nx*ny)  :: tf_clm_ij  ! sst climatology at target grids (nx*ny)
 real,    dimension(nx*ny)  :: tf_trd_ij  ! 6-hourly sst climatology tendency 
 real :: wei1,wei2
 integer :: nxc,nyc,mon1,mon2,i,j
 character (len=6), parameter :: fin_tf_clm='sstclm' ! sst climatology file name
!
! get which two months used and their weights from atime
!
 call get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
!
! get the dimensions of the sst climatology & allocate the related arrays
!
 call get_tf_clm_dim(fin_tf_clm,nyc,nxc)
 allocate( tf_clm0(nxc,nyc),tf_trd0(nxc,nyc),cxlats(nyc),cxlons(nxc) )
!
! get tf_clm at the analysis time from monthly climatology & cxlats, cxlons
!
 call get_tf_clm_ta(tf_clm0,tf_trd0,cxlats,cxlons,nyc,nxc,mon1,mon2,wei1,wei2)
!
! get tf_clm (nx by ny lat/lon) valid at atime
!
 if ( nx == nxc .and. ny == nyc ) then
    tf_clm(:,:) = tf_clm0(:,:)
    tf_trd(:,:) = tf_trd0(:,:)
!   write(*,'(a,2f9.3)') 'same dimensions, tf_clm, min : ',minval(tf_clm),maxval(tf_clm)
 else
!   write(*,'(a,4i8)') 'different dimensions,nx,ny,nxc,nyc : ',nx,ny,nxc,nyc
    call intp_tile(tf_clm0,  cxlats,  cxlons,  nyc, nxc, &
                   tf_clm_ij,xlats_ij,xlons_ij,ny,  nx)
    call intp_tile(tf_trd0,  cxlats,  cxlons,  nyc, nxc, &
                   tf_trd_ij,xlats_ij,xlons_ij,ny,  nx)
!   write(*,'(a,2f9.3)') 'tf_clm0, min, max                        : ',minval(tf_clm0),maxval(tf_clm0)

    tf_clm(:,:) = reshape (tf_clm_ij, (/nx,ny/) )
    tf_trd(:,:) = reshape (tf_trd_ij, (/nx,ny/) )
 endif

end subroutine get_tf_clm

!> Get the reference temperature/sst climatology and its trend at analysis time.
!! The data is time interpolated between two bounding months.
!!
!! @param[out] tf_clm_ta Climatological tf/sst at analysis time
!! @param[out] tf_clm_trend  Climatological tf/sst trend at analysis time
!! @param[out] xlats Latitudes on the climatological data grid
!! @param[out] xlons Longitudes on the climatological data grid
!! @param[in] nlat 'j' dimension on the climatological grid
!! @param[in] nlon 'i' dimension on the climatological grid
!! @param[in] mon1 First bounding month
!! @param[in] mon2 Second bounding month
!! @param[in] wei1 Weighting of first bounding month
!! @param[in] wei2 Weighting of second bounding month
!! @author Xu Li @date March 2019 
subroutine get_tf_clm_ta(tf_clm_ta,tf_clm_trend,xlats,xlons,nlat,nlon,mon1,mon2,wei1,wei2)
 use read_write_data, only : read_tf_clim_grb
 implicit none

! input
 integer, intent(in) :: nlat,nlon,mon1,mon2
 real,    intent(in) :: wei1,wei2
! output
 real, dimension(nlon,nlat), intent(out)   :: tf_clm_ta,tf_clm_trend
 real, dimension(nlat),      intent(out)   :: xlats
 real, dimension(nlon),      intent(out)   :: xlons

!input/output data file names
 character (len=6),  parameter :: fin_tf_clm='sstclm'

! local declare
 real, dimension(nlon,nlat) :: tf_clm1,tf_clm2

!
! read in rtg sst climatology without bitmap (surface mask) for mon1 and mon2
!
  call read_tf_clim_grb(trim(fin_tf_clm),tf_clm1,xlats,xlons,nlat,nlon,mon1)
  call read_tf_clim_grb(trim(fin_tf_clm),tf_clm2,xlats,xlons,nlat,nlon,mon2)
!
!  tf_clim at the analysis time
!
   tf_clm_ta(:,:) = wei1*tf_clm1(:,:)+wei2*tf_clm2(:,:)
!
!  tf_clim trend at the analysis time (month to month)
!
   tf_clm_trend(:,:) = (tf_clm2(:,:)-tf_clm1(:,:))/120.0

   write(*,'(a,2f9.3)') 'tf_clm_ta, min, max : ',minval(tf_clm_ta),maxval(tf_clm_ta)
   write(*,'(a,2f9.3)') 'tf_clm_trend, min, max : ',minval(tf_clm_trend),maxval(tf_clm_trend)
 end subroutine get_tf_clm_ta

 !> Get salinity climatology at the valid time on the target grid.
 !!
 !! @param[in] xlats_ij Latitudes of target grid
 !! @param[in] xlons_ij Longitudes of target grid
 !! @param[in] ny 'j' dimension of target grid
 !! @param[in] nx 'i' dimension of target grid
 !! @param[in] iy Year
 !! @param[in] im Month
 !! @param[in] id Day
 !! @param[in] ih Hour
 !! @param[out] sal_clm Salinity climatology on the target grid at the valid time
 !! @author Xu Li
subroutine get_sal_clm(xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,sal_clm)
 use read_write_data, only : get_dim_nc
 implicit none

 real,    dimension(nx*ny), intent(in)  :: xlats_ij   ! 
 real,    dimension(nx*ny), intent(in)  :: xlons_ij   ! 
 real,    dimension(nx,ny), intent(out) :: sal_clm    ! 
 integer, intent(in) :: iy,im,id,ih,nx,ny
! local declare
 real,    allocatable, dimension(:,:)   :: sal_clm0   ! salinity climatology at the valid time
 real,    allocatable, dimension(:)     :: cxlats     ! latitudes of sst climatology
 real,    allocatable, dimension(:)     :: cxlons     ! longitudes of sst climatology

 real,    dimension(nx*ny)  :: sal_clm_ij  ! salinity climatology at target grids (nx*ny)
 real :: wei1,wei2
 integer :: nxc,nyc,mon1,mon2,i,j
 character (len=6), parameter :: fin_sal_clm='salclm' ! salinity climatology file name
!
! get which two months used and their weights from atime
!
 call get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
!
! get the dimensions of the sst climatology & allocate the related arrays
!
 call get_dim_nc(fin_sal_clm,nyc,nxc)
 allocate( sal_clm0(nxc,nyc),cxlats(nyc),cxlons(nxc) )
!
! get sal_clm at the analysis time from monthly climatology & cxlats, cxlons
!
 call get_sal_clm_ta(sal_clm0,cxlats,cxlons,nyc,nxc,mon1,mon2,wei1,wei2)
!
! get sal_clm (nx by ny lat/lon) valid at atime
!
 if ( nx == nxc .and. ny == nyc ) then
    sal_clm(:,:) = sal_clm0(:,:)
!   write(*,'(a,2f9.3)') 'same dimensions, sal_clm, min : ',minval(sal_clm),maxval(sal_clm)
 else
!   write(*,'(a,4i8)') 'different dimensions,nx,ny,nxc,nyc : ',nx,ny,nxc,nyc
    call intp_tile(sal_clm0,  cxlats,  cxlons,  nyc, nxc, &
                   sal_clm_ij,xlats_ij,xlons_ij,ny,  nx)
!   write(*,'(a,2f9.3)') 'sal_clm0, min, max                        : ',minval(sal_clm0),maxval(sal_clm0)
!   write(*,'(a,2f9.3)') 'done with intp_tile for sal_clm, min, max : ',minval(sal_clm_ij),maxval(sal_clm_ij)

    sal_clm(:,:) = reshape (sal_clm_ij, (/nx,ny/) )
 endif

end subroutine get_sal_clm

!> Get climatological salinity at the analysis time.
!!
!! @param[in] nlat 'j' dimension of climatological data
!! @param[in] nlon 'i' dimension of climatological data
!! @param[in] mon1 First bounding month
!! @param[in] mon2 Second bounding month
!! @param[in] wei1 Weight of first bounding month
!! @param[in] wei2 Weight of second bounding month
!! @param[out] sal_clm_ta Climatological salinity at the analysis time
!! @param[out] xlats Latitudes on the climatological grid
!! @param[out] xlons Longitudes on the climatological grid
!! @author Xu Li @date March 2019 
subroutine get_sal_clm_ta(sal_clm_ta,xlats,xlons,nlat,nlon,mon1,mon2,wei1,wei2)

 use read_write_data, only : read_salclm_gfs_nc
 implicit none

! input
 integer, intent(in) :: nlat,nlon,mon1,mon2
 real,    intent(in) :: wei1,wei2
! output
 real, dimension(nlon,nlat), intent(out)   :: sal_clm_ta
 real, dimension(nlat),      intent(out)   :: xlats
 real, dimension(nlon),      intent(out)   :: xlons

!input/output data file names
 character (len=6),  parameter :: fin_sal_clm='salclm'

! local declare
 real, dimension(nlon,nlat) :: sal_clm1,sal_clm2

!
! read in rtg sst climatology without bitmap (surface mask) for mon1 and mon2
!
  call read_salclm_gfs_nc(trim(fin_sal_clm),sal_clm1,xlats,xlons,nlat,nlon,mon1)
  call read_salclm_gfs_nc(trim(fin_sal_clm),sal_clm2,xlats,xlons,nlat,nlon,mon2)
!
!  sal_clim at the analysis time
!
   sal_clm_ta(:,:) = wei1*sal_clm1(:,:)+wei2*sal_clm2(:,:)
   write(*,'(a,2f9.3)') 'sal_clm_ta, min, max : ',minval(sal_clm_ta),maxval(sal_clm_ta)
 end subroutine get_sal_clm_ta

 !> Interpolate lon/lat grid data to the fv3 native grid (tf_lalo => tf_tile). Does not
 !! account for a mask.
 !!
 !! @param[in] tf_lalo (idim_lalo,idim_lalo) field on the lat/lon regular grid.
 !! @param[in] dlats_lalo (jdim_lalo) latitudes along y direction of lat/lon regular grid points.
 !! @param[in] dlons_lalo (idim_lalo) longitudes along x direction of lat/lon regular grid points.
 !! @param[in] jdim_lalo number of y dimension of tf_lalo.
 !! @param[in] idim_lalo number of x dimension of tf_lalo.
 !! @param[in] xlats_tile (jdim_tile*idim_tile) latitudes of all tile grid points.
 !! @param[in] xlons_tile (jdim_tile*idim_tile) longitudes of all tile grid points.
 !! @param[in] jdim_tile number of y dimension of tf_tile.
 !! @param[in] idim_tile number of x dimension of tf_tile.
 !! @param[out] tf_tile (jdim_tile*idim_tile) field on the cubed sphere grid.
 !! @author Xu Li
subroutine intp_tile(tf_lalo,dlats_lalo,dlons_lalo,jdim_lalo,idim_lalo, &
                     tf_tile,xlats_tile,xlons_tile,jdim_tile,idim_tile)

 use utils

 implicit none

! input/output
 real,    dimension(idim_lalo,jdim_lalo), intent(in)  :: tf_lalo
 real,    dimension(jdim_lalo),           intent(in)  :: dlats_lalo
 real,    dimension(idim_lalo),           intent(in)  :: dlons_lalo
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlats_tile
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlons_tile
 integer,                                 intent(in)  :: jdim_lalo,idim_lalo,jdim_tile,idim_tile
 real,    dimension(jdim_tile*idim_tile), intent(out) :: tf_tile

! local
 real, parameter :: deg2rad=3.1415926/180.0
 real,    dimension(jdim_lalo) :: xlats_lalo
 real,    dimension(idim_lalo) :: xlons_lalo
 real    :: tf,wsum,res_km
 integer :: itile,jtile
 integer :: ii,jj,ij,iii,jjj
 integer :: ilalo,jlalo,ilalop1,jlalop1
 integer :: istart,iend,jstart,jend,krad

 integer, allocatable, dimension(:,:)   :: id1,id2,jdc
 real,    allocatable, dimension(:,:,:) :: agrid,s2c

 print*
 print*,'interpolate from lat/lon grids to any one grid with known lat/lon'

 xlats_lalo = dlats_lalo*deg2rad
 xlons_lalo = dlons_lalo*deg2rad

 allocate(agrid(idim_tile,jdim_tile,2))
 agrid(:,:,1) = reshape (xlons_tile, (/idim_tile,jdim_tile/) )
 agrid(:,:,2) = reshape (xlats_tile, (/idim_tile,jdim_tile/) )
 agrid        = agrid*deg2rad

 allocate(id1(idim_tile,jdim_tile))
 allocate(id2(idim_tile,jdim_tile))
 allocate(jdc(idim_tile,jdim_tile))
 allocate(s2c(idim_tile,jdim_tile,4))

!----------------------------------------------------------------------
! compute bilinear weights for each model point from the nearest
! four lalo points. does not account for mask.  that
! happens later.
!----------------------------------------------------------------------

 call remap_coef( 1, idim_tile, 1, jdim_tile, idim_lalo, jdim_lalo, &
                  xlons_lalo, xlats_lalo, id1, id2, jdc, s2c, agrid )

 do ij = 1, jdim_tile*idim_tile

    jtile = (ij-1)/idim_tile + 1
    itile = mod(ij,idim_tile)
    if (itile==0) itile = idim_tile

    ilalo   = id1(itile,jtile)
    ilalop1 = id2(itile,jtile)
    jlalo   = jdc(itile,jtile)
    jlalop1 = jdc(itile,jtile) + 1

    wsum = s2c(itile,jtile,1) + s2c(itile,jtile,2) + &
           s2c(itile,jtile,3) + s2c(itile,jtile,4)

    tf_tile(ij)  = ( s2c(itile,jtile,1)*tf_lalo(ilalo,jlalo)     + &
                     s2c(itile,jtile,2)*tf_lalo(ilalop1,jlalo)   + &
                     s2c(itile,jtile,3)*tf_lalo(ilalop1,jlalop1) + &
                     s2c(itile,jtile,4)*tf_lalo(ilalo,jlalop1) )/wsum
 enddo

 deallocate(id1, id2, jdc, s2c)

end subroutine intp_tile

!> For a given date, determine the bounding months and the linear
!! time interpolation weights.
!!
!! @param[in] iy The year
!! @param[in] im The month
!! @param[in] id The day
!! @param[in] ih The hour
!! @param[out] mon1 First bounding month
!! @param[out] mon2 Second bounding month
!! @param[out] wei1 Weighting of first bounding month
!! @param[out] wei2 Weighting of second bounding month
!! @author Xu Li @date March 2019 
subroutine get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
 implicit none

! input
 integer, intent(in) :: iy,im,id,ih
! output
 integer, intent(out) :: mon1,mon2
 real,    intent(out) :: wei1,wei2

! local declare
 real :: rjday
 integer :: mon,monend,monm,monp,jdow,jdoy,jday
 integer :: jda(8)
!
!dayhf : julian day of the middle of each month
!
 real, dimension(13) ::  dayhf
 data dayhf/15.5,45.0,74.5,105.0,135.5,166.0,196.5,227.5,258.0,288.5,319.0,349.5,380.5/

! 15, 44, 73.5, 104, 134.5, 165, 195.5, 226.5, 257, 287.5, 318.5, 349  ! from
! woa05

 jda=0
 jda(1)=iy
 jda(2)=im
 jda(3)=id
 jda(5)=ih
 jdow = 0
 jdoy = 0
 jday = 0
 call w3doxdat(jda,jdow,jdoy,jday)
 rjday=jdoy+jda(5)/24.
 if(rjday.lt.dayhf(1)) rjday=rjday+365.

 monend = 12
 do mon = 1, monend
    monm = mon
    monp = mon + 1
    if( rjday >= dayhf(monm) .and. rjday < dayhf(monp) ) then
       mon1 = monm
       mon2 = monp
       go to 10
    endif
 enddo

 print *,'FATAL ERROR in get_tim_wei, wrong rjday',rjday
 call abort
 10     continue

 wei1 = (dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
 wei2 = (rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))

 if( mon2 == 13 ) mon2=1

 write(*,'(a,2i4,3f9.3)') 'mon1,mon2,rjday,wei1,wei2=',mon1,mon2,rjday,wei1,wei2

 end subroutine get_tim_wei

 !> Compute the freezing point of water as a function of salinity.
 !!
 !! Constants taken from Gill, 1982.
 !!
 !! @date 21 September 1994.
 !! @author Robert Grumbine
 !!
 !! @param [in] salinity The salinity.
 !! @return tfreez The freezing point of water.
 real function tfreez(salinity)

 implicit none

 REAL salinity,sal
 REAL a1, a2, a3
 PARAMETER (a1 = -0.0575)
 PARAMETER (a2 =  1.710523E-3)
 PARAMETER (a3 = -2.154996E-4)

 IF (salinity .LT. 0.) THEN
   sal = 0.
 ELSE
   sal = salinity
 ENDIF
 tfreez = sal*(a1+a2*SQRT(sal)+a3*sal)

 return
 end
