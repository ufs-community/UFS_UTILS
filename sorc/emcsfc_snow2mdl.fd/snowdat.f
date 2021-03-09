!> @file
!! @brief Read and qc afwa, nesdis/ims and autosnow snow data.
!! @author gayno org: w/np2 @date 2005-dec-16 

!> Read and qc afwa, nesdis/ims and autosnow snow data.
!!
!! program history log:
!! -  2005-dec-16  gayno   - initial version
!! -  2007-aug-10  gayno   - Allow program to run with no nesdis/ims data.
!!                          Add 16th mesh afwa grib data.
!! -  2008-feb-04  gayno   - Add autosnow cover data for sh.
!! -  2009-jun-03  gayno   - Add qc check for nesdis/ims and afwa data.
!! -  2014-feb-07  gayno   - Read nesdis/ims data in grib1 or grib 2
!!                          format.
!! -  2014-sep-30  gayno   - Convert weekly nh snow climatology - used to 
!!                          qc input data - to grib 2.
!! variable definitions:
!! -  bad_afwa_Xh        -  is afwa data corrupt?
!! -  bitmap_afwa_Xh     -  bitmap of afwa grid (false-non land, true-land)
!! -  kgds_afwa_Xh       - afwa grid description section (grib section 2)
!! -  nesdis_res         - resolution of nesdis/ims data in km
!! -  snow_dep_afwa_Xh   - afwa snow depth data (inches*10 on input, 
!!                                              meters on output)
!! -  use_xh_afwa        - true if afwa data to be used
!!
 module snowdat

 use program_setup, only  : autosnow_file,       &
                            nesdis_snow_file,    &
                            nesdis_lsmask_file,  &
                            afwa_snow_global_file, &
                            afwa_snow_nh_file,   &
                            afwa_snow_sh_file,   &
                            afwa_lsmask_nh_file, &
                            afwa_lsmask_sh_file

 use model_grid, only     : imdl,                &
                            jmdl

 integer                 :: iafwa  !< i-dimension of afwa grid
 integer                 :: iautosnow  !< i-dimension of autosnow grid
 integer                 :: inesdis   !< i-dimension of nesdis grid
 integer                 :: jafwa  !< j-dimension of afwa grid
 integer                 :: jautosnow  !< j-dimension of autosnow grid
 integer                 :: jnesdis  !< j-dimension of nesdis grid
 integer                 :: kgds_afwa_global(200) !< grib1 grid description section for
                                                  !! global afwa data.
 integer                 :: kgds_afwa_nh(200) !< grib1 grid description section for northern
                                              !! hemisphere 16th mesh afwa data.
 integer                 :: kgds_afwa_nh_8th(200) !< grib1 grid description section for
                                                  !! northern hemisphere 8th mesh afwa data.
 integer                 :: kgds_afwa_sh(200) !< grib1 grid description section for southern
                                              !! hemisphere 16th mesh afwa data.
 integer                 :: kgds_afwa_sh_8th(200) !< grib1 grid description section for
                                                  !! southern hemisphere 8th mesh afwa data.
 integer                 :: kgds_autosnow(200)  !< autosnow grid description section (grib section 2)
 integer                 :: kgds_nesdis(200)  !< nesdis/ims grid description section (grib section 2)
 integer                 :: mesh_nesdis !< nesdis/ims data is 96th mesh (or bediant)
 integer*1, allocatable  :: sea_ice_nesdis(:,:) !< nesdis/ims sea ice flag (0-open water, 1-ice) 
 logical                 :: bad_afwa_nh !< When true, the northern hemisphere afwa data failed
                                        !! its quality control check.
 logical                 :: bad_afwa_sh !< When true, the southern hemisphere afwa data failed
                                        !! its quality control check.
 logical                 :: bad_nesdis   !< When true, the nesdis data failed its quality 
                                         !! control check.
 logical                 :: bad_afwa_global !< When true, the global afwa data failed its quality
                                            !! control check.
 logical*1, allocatable  :: bitmap_afwa_global(:,:) !< The global afwa data grib bitmap.
                                                    !!(false-non land, true-land).
 logical*1, allocatable  :: bitmap_afwa_nh(:,:) !< The northern hemisphere afwa data grib bitmap.
                                                !! (false-non land, true-land).
 logical*1, allocatable  :: bitmap_afwa_sh(:,:) !< The southern hemisphere afwa data grib bitmap.
                                                !! (false-non land, true-land).
 logical*1, allocatable  :: bitmap_nesdis(:,:) !< nesdis data grib bitmap (false-non land, true-land).
 logical*1, allocatable  :: bitmap_autosnow(:,:) !< autosnow data grib bitmap (false-non land,
                                                 !! true-land).
 logical                 :: use_nh_afwa !< True if northern hemisphere afwa data to be used.
 logical                 :: use_sh_afwa !< True if southern hemisphere afwa data to be used.
 logical                 :: use_global_afwa !< True if global hemisphere afwa data to be used.
 logical                 :: use_autosnow !< True if autosnow data to be used
 logical                 :: use_nesdis  !< True if nesdis/ims data to be used

 real                    :: autosnow_res  !< Resolution of autosnow in km 
 real                    :: afwa_res   !<  Resolution of afwa data in km
 real                    :: nesdis_res !< Resolution of the nesdis data in km.
 real, allocatable       :: snow_cvr_nesdis(:,:)   !< nesdis/ims snow cover flag (0-no, 100-yes)
 real, allocatable       :: snow_cvr_autosnow(:,:)  !< autosnow snow cover flag (0-no, 100-yes)
 real, allocatable       :: snow_dep_afwa_global(:,:) !< The global afwa snow depth.
 real, allocatable       :: snow_dep_afwa_nh(:,:)  !< Northern hemisphere afwa snow depth.
 real, allocatable       :: snow_dep_afwa_sh(:,:)  !< Southern hemisphere afwa snow depth.

! the afwa 8th mesh binary data has no grib header, so set it from these
! data statements. needed for ipolates routines.

 data kgds_afwa_nh_8th/5,2*512,-20826,145000,8,-80000,2*47625,0,  &
                       9*0,255,180*0/
 data kgds_afwa_sh_8th/5,2*512,20826,-125000,8,-80000,2*47625,128, &
                       9*0,255,180*0/
 contains
!> Read autosnow snow cover.
!!
!! program history log:
!! 2008-feb-04  gayno    - initial version
!!
!! files:
!!   input:
!!     - autosnow data, grib 2, unit=lugb
!!
!! condition codes:  all fatal
!!   74 - bad open of autosnow file
!!   75 - bad read of autosnow file
!!
!! @note    Autosnow data is available only for southern hemis.
!!          Autosnow data is in grib 2.          
!!
!! @author  George Gayno   org: w/np2   @date  2008-Feb-04
 subroutine readautosnow
 use grib_mod  ! grib 2 libraries

 implicit none

 type(gribfield)            :: gfld

 integer                    :: iret, j, k, lugb, lugi
 integer                    :: jdisc, jgdtn, jpdtn
 integer                    :: jids(200), jgdt(200), jpdt(200)

 logical                    :: unpack

 use_autosnow = .true.

 if ( len_trim(autosnow_file) == 0 ) then
   print*,"- WILL NOT USE AUTOSNOW DATA."
   use_autosnow = .false.
   return
 end if

 print*,"- OPEN AND READ AUTOSNOW FILE ", trim(autosnow_file)

 lugb=12
 call baopenr(lugb,autosnow_file,iret)

 if (iret /= 0) then
   print*,'- FATAL ERROR: BAD OPEN OF FILE, IRET IS ', iret
   call w3tage('SNOW2MDL')
   call errexit(74)
 endif

 call grib2_null(gfld)

 j       = 0      ! search at beginning of file
 lugi    = 0      ! no grib index file
 jdisc   = 0      ! search for discipline; 0 - meteorological products
 jpdtn   = 30     ! search for product definition template number; 30 - satellite product
 jgdtn   = 0      ! search for grid definition template number; 0 - lat/lon grid
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definiation template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 jpdt(1) = 1      ! search for parameter category - moisture
 jpdt(2) = 42     ! search for parameter number - snow cover in percent.
 unpack  = .true. ! unpack data

 call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

 if (iret /=0) then
  print*,'- FATAL ERROR: BAD DEGRIB OF FILE, IRET IS ', iret
  call w3tage('SNOW2MDL')
  call errexit(75)
 endif

 print*,"- DATA VALID AT (YYYYMMDDHH): ", gfld%idsect(6),gfld%idsect(7), &
                                          gfld%idsect(8),gfld%idsect(9)

 call baclose (lugb, iret)

!-----------------------------------------------------------------------
! set the grib1 kgds array from the g2 grid definition template array. 
! the kgds array is used by ipolates.
!-----------------------------------------------------------------------

 call gdt_to_gds(gfld%igdtnum, gfld%igdtmpl, gfld%igdtlen, kgds_autosnow, &
                 iautosnow, jautosnow, autosnow_res)
 
 allocate (bitmap_autosnow(iautosnow,jautosnow))
 bitmap_autosnow = reshape (gfld%bmap , (/iautosnow,jautosnow/) )
 
 allocate (snow_cvr_autosnow(iautosnow,jautosnow))
 snow_cvr_autosnow = reshape (gfld%fld , (/iautosnow,jautosnow/) )

 call grib2_free(gfld)

 end subroutine readautosnow

!> Read nesdis/ims snow cover/ice data.
!!   
!! program history log:
!! 2005-dec-16  gayno    - initial version
!! 2014-feb-07  gayno    - Read 4km ims data in either
!!                         grib1 or grib 2 format.
!! files:
!!   input:
!!      - ims snow cover and ice file, grib 1 or grib 2
!!      - 16th-mesh ims land mask, binary
!!
!! condition codes: all fatal
!!   41 - ims file not grib 1 or grib 2
!!   53 - ims data failed quality check
!!   70 - bad read of ims snow cover data
!!   71 - bad read of ims ice data
!!   72 - bad read of ims grib 1 header
!!   73 - bad open of ims file
!!   87 - bad open ims land mask file
!!   88 - bad read ims land mask file
!!   
!! @note    Nesdis/ims data available only for n hemis.  Ims data used
!!          to be created by nesdis,  hence the references to "nesdis"
!!          in this routine.  Ims data is now created by the national
!!          ice center.
!!
!! @author  George Gayno org: w/np2 @date 2005-Dec-16
 subroutine readnesdis
 use grib_mod

 implicit none

 integer, parameter         :: iunit = 13  ! input grib file unit number
 integer, parameter         :: iunit2 = 43  ! input landmask file unit number

 integer*4, allocatable     :: dummy4(:,:)
 integer                    :: i, j
 integer                    :: iret
 integer                    :: jgds(200)
 integer                    :: jpds(200)
 integer                    :: lskip
 integer, parameter         :: lugi = 0    ! grib index file unit number - not used
 integer                    :: jdisc, jgdtn, jpdtn, k
 integer                    :: jids(200), jgdt(200), jpdt(200)
 integer                    :: kgds(200)
 integer                    :: kpds(200)
 integer                    :: message_num
 integer                    :: numbytes
 integer                    :: numpts
 integer                    :: isgrib

 logical                    :: unpack

 real, allocatable          :: dummy(:,:)
 real                       :: dum
 
 type(gribfield)            :: gfld

 use_nesdis = .true.

 if ( len_trim(nesdis_snow_file) == 0 ) then
   print*,"- WILL NOT USE NESDIS/IMS DATA."
   use_nesdis = .false.
   return
 end if

 print*,"- OPEN AND READ NESDIS/IMS SNOW FILE ", trim(nesdis_snow_file)

 call grib_check(nesdis_snow_file, isgrib)

 if (isgrib==0) then
   print*,'- FATAL ERROR: IMS FILE MUST BE GRIB 1 OR GRIB2 FORMAT'
   call w3tage('SNOW2MDL')
   call errexit(41)
 end if

 call baopenr (iunit, nesdis_snow_file, iret)

 if (iret /= 0) then
   print*,'- FATAL ERROR: BAD OPEN OF FILE, IRET IS ', iret
   call w3tage('SNOW2MDL')
   call errexit(73)
 end if

 if (isgrib==1) then  ! grib 1 format

!-----------------------------------------------------------------------
! tell degribber to look for requested data.
!-----------------------------------------------------------------------

   lskip    = -1
   jpds     = -1
   jgds     = -1
   jpds(5)  = 91     ! ice cover
   kpds     = jpds
   kgds     = jgds

   print*,"- GET GRIB HEADER"

   call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
               numpts, message_num, kpds, kgds, iret)

   if (iret /= 0) then
     print*,"- FATAL ERROR: BAD DEGRIB OF HEADER. IRET IS ", iret
     call w3tage('SNOW2MDL')
     call errexit(72)
   end if

   kgds_nesdis = kgds
   inesdis     = kgds(2)
   jnesdis     = kgds(3)

   mesh_nesdis = inesdis / 64
   nesdis_res  = 381. / float(mesh_nesdis)   ! in km

   print*,"- DATA VALID AT (YYMMDDHH): ", kpds(8:11)
 
   allocate (dummy(inesdis,jnesdis))
   allocate (sea_ice_nesdis(inesdis,jnesdis))
   allocate (bitmap_nesdis(inesdis,jnesdis))

   print*,"- DEGRIB SEA ICE."

   call getgb(iunit, lugi, (inesdis*jnesdis), lskip, jpds, jgds, &
              numpts, lskip, kpds, kgds, bitmap_nesdis, dummy, iret)

   if (iret /= 0) then
     print*,"- FATAL ERROR: BAD DEGRIB OF DATA. IRET IS ", iret
     call w3tage('SNOW2MDL')
     call errexit(71)
   end if

   sea_ice_nesdis = nint(dummy)  ! only needed as yes/no flag
   deallocate (dummy)

   lskip    = -1
   jpds     = -1
   jgds     = -1
   jpds(5)  = 238     ! snow cover
   kpds     = jpds
   kgds     = jgds

   allocate (snow_cvr_nesdis(inesdis,jnesdis))

   print*,"- DEGRIB SNOW COVER."

   call getgb(iunit, lugi, (inesdis*jnesdis), lskip, jpds, jgds, &
              numpts, lskip, kpds, kgds, bitmap_nesdis, snow_cvr_nesdis, iret)

   if (iret /= 0) then
     print*,"- FATAL ERROR: BAD DEGRIB OF DATA. IRET IS ", iret
     call w3tage('SNOW2MDL')
     call errexit(70)
   end if

 elseif (isgrib==2) then  ! grib 2 format

   print*,"- DEGRIB SNOW COVER."

   j       = 0      ! search at beginning of file
   jdisc   = 0      ! search for discipline; 0 - meteorological products
   jpdtn   = 0      ! search for product definition template number; 0 - analysis at one level
   jgdtn   = 20     ! search for grid definition template number; 20 - polar stereographic grid
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template 3.m
   jpdt    = -9999  ! array of values in product definition template 4.n
   jpdt(1) = 1      ! search for parameter category - moisture
   jpdt(2) = 201    ! search for parameter number - snow cover in percent.
   unpack  = .true. ! unpack data

   call grib2_null(gfld)

   call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, iret)

   if (iret /=0) then
    print*,'- FATAL ERROR: BAD DEGRIB OF FILE, IRET IS ', iret
    call w3tage('SNOW2MDL')
    call errexit(70)
   endif

   print*,"- DATA VALID AT (YYYYMMDDHH): ", gfld%idsect(6),gfld%idsect(7), &
                                            gfld%idsect(8),gfld%idsect(9)

!-----------------------------------------------------------------------
! set the grib1 kgds array from the g2 grid definition template array. 
! the kgds array is used by ipolates.
!-----------------------------------------------------------------------

   call gdt_to_gds(gfld%igdtnum, gfld%igdtmpl, gfld%igdtlen, kgds_nesdis, &
                   inesdis, jnesdis, dum)

   mesh_nesdis = inesdis / 64
   nesdis_res  = 381. / float(mesh_nesdis)   ! in km

   if (mesh_nesdis==16) kgds_nesdis(6)=136  ! the ims 16th mesh grib2 data
                                            ! is gribbed with an elliptical
                                            ! earth.  that is wrong. hardwire
                                            ! a fix here.

   allocate (snow_cvr_nesdis(inesdis,jnesdis))
   allocate (sea_ice_nesdis(inesdis,jnesdis))
   allocate (bitmap_nesdis(inesdis,jnesdis))

   bitmap_nesdis = reshape (gfld%bmap , (/inesdis,jnesdis/) )
   snow_cvr_nesdis = reshape (gfld%fld , (/inesdis,jnesdis/) )
   
   call grib2_free(gfld)

   print*,"- DEGRIB SEA ICE."

   j       = 0      ! search at beginning of file
   jdisc   = 10     ! search for discipline; 10 - ocean products
   jpdtn   = 0      ! search for product definition template number; 0 - analysis at one level
   jgdtn   = 20     ! search for grid definition template number; 20 - polar stereographic grid
   jids    = -9999  ! array of values in identification section, set to wildcard
   jgdt    = -9999  ! array of values in grid definition template 3.m
   jpdt    = -9999  ! array of values in product definition template 4.n
   jpdt(1) = 2      ! search for parameter category - ice
   jpdt(2) = 0      ! search for parameter number - ice cover in percent.
   unpack  = .true. ! unpack data

   call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld, iret)

   if (iret /=0) then
    print*,'- FATAL ERROR: BAD DEGRIB OF FILE, IRET IS ', iret
    call w3tage('SNOW2MDL')
    call errexit(71)
   endif

   sea_ice_nesdis = reshape (gfld%fld , (/inesdis,jnesdis/) )

   call grib2_free(gfld)

 end if

 call baclose(iunit,iret)

!-----------------------------------------------------------------------
! the 16th mesh nesdis/ims grib data does not have a proper
! bitmap section.  therefore, need to read in the mask
! from file.  but the 96th mesh data has a proper bitmap, so use it.
!-----------------------------------------------------------------------

 if (mesh_nesdis == 16) then

   print*,"- OPEN NESDIS/IMS 16TH MESH LAND MASK: ", trim(nesdis_lsmask_file)

   open(iunit2, file=trim(nesdis_lsmask_file), form="formatted", &
        iostat = iret)

   if (iret /= 0) then
     print*,"- FATAL ERROR OPENING NESDIS/IMS LAND MASK FILE. ISTAT IS: ", iret
     call errexit(87)
   end if

   print*,"- READ NESDIS/IMS 16TH MESH LAND MASK."

   allocate (dummy4(inesdis,jnesdis))
   
   do j = 1, 1024
     read(iunit2, 123, iostat=iret) (dummy4(i,j),i=1,1024)
     if (iret /= 0) then
       print*,"- FATAL ERROR READING NESDIS/IMS LAND MASK FILE. ISTAT IS: ", iret
       call errexit(88)
     end if
   enddo

   close (iunit2)

!-----------------------------------------------------------------------
! the file has 0-sea, 1-land, 9-off hemi.  this code expects
! 0-non-land (or don't use data), 1-land (use data).
!-----------------------------------------------------------------------

   bitmap_nesdis=.false.
   do j = 1, 1024
     do i = 1, 1024
       if (dummy4(i,j) == 1) bitmap_nesdis(i,j) = .true.
     enddo
   enddo
 
   deallocate(dummy4)

123 FORMAT(80I1)

 endif  ! is nesdis/ims data 16th mesh?

 bad_nesdis=.false.
 call nh_climo_check(kgds_nesdis,snow_cvr_nesdis,bitmap_nesdis,inesdis,jnesdis,2,bad_nesdis)

!-----------------------------------------------------------------------
! for the 2009 nmm-b implementation, it was decided to not run with
! afwa only.  so even if afwa data is current and not corrupt, 
! but the ims is bad, then abort program. exception, if ims is very old
! (there is a catastropic outage) then program will run with afwa
! only.  this is done by setting the nesdis_snow_file variable to
! a zero length string (i.e., ims data not selected).  this variable
! setting is accomplished in the run script. 
!-----------------------------------------------------------------------

 if (bad_nesdis) then
   print*,'- FATAL ERROR: NESDIS/IMS DATA BAD, DO NOT USE.'
   print*,'- DONT RUN PROGRAM.'
   use_nesdis=.false.
   call w3tage('SNOW2MDL')
   call errexit(53)
   stop
 endif

 return

 end subroutine readnesdis

!>  Read snow depth data and masks.
!!
!! @note Read nh and sh afwa snow depth data and
!!   land sea mask. 
!!
!! program history log:
!!
!! 2005-dec-16  gayno    - initial version
!! 2007-nov-28  gayno    - read 16th mesh afwa data in grib format
!!
!! files:
!!   input:
!!     - global afwa data in grib 1 (if selected)
!!     - nh afwa data in grib 1 (if selected)
!!     - sh afwa data in grib 1 (if selected)
!!
!! condition codes:
!!   60 - bad open afwa file
!!   61 - bad degrib of afwa file
!!
!! @author  George Gayno org: w/np2 @date 2005-Dec-16
 subroutine readafwa
 implicit none

 integer, parameter            :: iunit=11
 integer                       :: jgds(200), jpds(200), kgds(200), kpds(200)
 integer                       :: istat
 integer                       :: lugi, lskip, numbytes, numpts, message_num
 integer                       :: isgrib

 bad_afwa_nh=.false.
 bad_afwa_sh=.false.
 bad_afwa_global=.false.

 use_global_afwa=.true.
 use_nh_afwa = .true.
 use_sh_afwa = .true.

 if (len_trim(afwa_snow_nh_file) == 0 .and.   &
     len_trim(afwa_snow_sh_file) == 0 .and.   &
     len_trim(afwa_snow_global_file) == 0) then
   print*,"- WILL NOT USE AFWA DATA."
   use_nh_afwa = .false.
   use_sh_afwa = .false.
   use_global_afwa = .false.
   return
 end if

 if ( len_trim(afwa_snow_global_file) > 0 ) then

   print*,"- OPEN AND READ AFWA SNOW FILE ", trim(afwa_snow_global_file)
   call baopenr (iunit, afwa_snow_global_file, istat)
   if (istat /= 0) then
     print*,'- FATAL ERROR: BAD OPEN OF FILE, ISTAT IS ', istat
     call w3tage('SNOW2MDL')
     call errexit(60)
   end if

!-----------------------------------------------------------------------
! tell degribber to look for requested data.
!-----------------------------------------------------------------------

   lugi     = 0
   lskip    = -1
   jpds     = -1
   jgds     = -1
   jpds(5)  = 66     ! snow depth
   kpds     = jpds
   kgds     = jgds

   print*,"- GET GRIB HEADER"
   call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
               numpts, message_num, kpds, kgds, istat)

   if (istat /= 0) then
     print*,"- FATAL ERROR: BAD DEGRIB OF HEADER. ISTAT IS ", istat
     call w3tage('SNOW2MDL')
     call errexit(61)
   end if

   iafwa = kgds(2)
   jafwa = kgds(3)
   afwa_res = float(kgds(10))*0.001*111.0  ! in km.  

   print*,"- DATA VALID AT (YYMMDDHH): ", kpds(8:11)
   print*,"- DEGRIB SNOW DEPTH."

   allocate(bitmap_afwa_global(iafwa,jafwa))
   allocate(snow_dep_afwa_global(iafwa,jafwa))

   call getgb(iunit, lugi, (iafwa*jafwa), lskip, jpds, jgds, &
              numpts, lskip, kpds, kgds, bitmap_afwa_global, snow_dep_afwa_global, istat)

   if (istat /= 0) then
     print*,"- FATAL ERROR: BAD DEGRIB OF DATA. ISTAT IS ", istat
     call w3tage('SNOW2MDL')
     call errexit(61)
   end if

   kgds_afwa_global = kgds

   call baclose(iunit, istat) 

   call nh_climo_check(kgds_afwa_global,snow_dep_afwa_global,bitmap_afwa_global,iafwa,jafwa,1,bad_afwa_global)

   if (bad_afwa_global) then
     print*,'- WARNING: AFWA DATA BAD, DO NOT USE.'
     use_global_afwa = .false.
   endif

   use_nh_afwa=.false.   ! use global or hemispheric files. not both.
   use_sh_afwa=.false.

   return  ! use global or hemispheric files. not both.

 else

   use_global_afwa=.false.

 endif
   
 if ( len_trim(afwa_snow_nh_file) > 0 ) then  ! afwa nh data selected

   call grib_check(afwa_snow_nh_file, isgrib)

   if (isgrib==0) then ! old ncep binary format

     iafwa = 512
     jafwa = 512
     afwa_res = 47.625   ! in kilometers
     kgds_afwa_nh = kgds_afwa_nh_8th

     allocate (snow_dep_afwa_nh(iafwa,jafwa))
     call read_afwa_binary(afwa_snow_nh_file, snow_dep_afwa_nh)

     allocate (bitmap_afwa_nh(iafwa,jafwa))
     call read_afwa_mask(afwa_lsmask_nh_file, bitmap_afwa_nh) 

   else ! afwa data is grib

     print*,"- OPEN AND READ AFWA SNOW FILE ", trim(afwa_snow_nh_file)

     call baopenr (iunit, afwa_snow_nh_file, istat)

     if (istat /= 0) then
       print*,'- FATAL ERROR: BAD OPEN OF FILE, ISTAT IS ', istat
       call w3tage('SNOW2MDL')
       call errexit(60)
     end if

!-----------------------------------------------------------------------
! tell degribber to look for requested data.
!-----------------------------------------------------------------------

     lugi     = 0
     lskip    = -1
     jpds     = -1
     jgds     = -1
     jpds(5)  = 66     ! snow depth
     kpds     = jpds
     kgds     = jgds

     print*,"- GET GRIB HEADER"
     call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
                 numpts, message_num, kpds, kgds, istat)

     if (istat /= 0) then
       print*,"- FATAL ERROR: BAD DEGRIB OF HEADER. ISTAT IS ", istat
       call w3tage('SNOW2MDL')
       call errexit(61)
     end if

     iafwa = kgds(2)
     jafwa = kgds(3)
     afwa_res = float(kgds(8))*0.001  ! in km.  

     print*,"- DATA VALID AT (YYMMDDHH): ", kpds(8:11)

     print*,"- DEGRIB SNOW DEPTH."

     allocate(bitmap_afwa_nh(iafwa,jafwa))
     allocate(snow_dep_afwa_nh(iafwa,jafwa))

     call getgb(iunit, lugi, (iafwa*jafwa), lskip, jpds, jgds, &
                numpts, lskip, kpds, kgds, bitmap_afwa_nh, snow_dep_afwa_nh, istat)

     if (istat /= 0) then
       print*,"- FATAL ERROR: BAD DEGRIB OF DATA. ISTAT IS ", istat
       call w3tage('SNOW2MDL')
       call errexit(61)
     end if

     kgds_afwa_nh = kgds

     kgds_afwa_nh(7) = -80000  ! ipolates definition of orientation angle is
                               ! 180 degrees off from grib standard.

     call baclose(iunit, istat) 

   endif ! is nh afwa data grib?

   call nh_climo_check(kgds_afwa_nh,snow_dep_afwa_nh,bitmap_afwa_nh,iafwa,jafwa,1,bad_afwa_nh)

 else

   use_nh_afwa=.false.

 endif

!-----------------------------------------------------------------------
! now, read southern hemisphere data.
!-----------------------------------------------------------------------

 if ( len_trim(afwa_snow_sh_file) > 0 ) then

   call grib_check(afwa_snow_sh_file, isgrib)

   if (isgrib==0) then ! old ncep binary format

     iafwa = 512
     jafwa = 512
     afwa_res = 47.625
     kgds_afwa_sh = kgds_afwa_sh_8th

     allocate (snow_dep_afwa_sh(iafwa,jafwa))
     call read_afwa_binary(afwa_snow_sh_file, snow_dep_afwa_sh)

     allocate (bitmap_afwa_sh(iafwa,jafwa))
     call read_afwa_mask(afwa_lsmask_sh_file, bitmap_afwa_sh) 

   else   ! sh afwa data is grib

     print*,"- OPEN AND READ AFWA SNOW FILE ", trim(afwa_snow_sh_file)

     call baopenr (iunit, afwa_snow_sh_file, istat)

     if (istat /= 0) then
       print*,'- FATAL ERROR: BAD OPEN OF FILE, ISTAT IS ', istat
       call w3tage('SNOW2MDL')
       call errexit(60)
     end if

!-----------------------------------------------------------------------
! tell degribber to look for requested data.
!-----------------------------------------------------------------------

     lugi     = 0
     lskip    = -1
     jpds     = -1
     jgds     = -1
     jpds(5)  = 66     ! snow cover
     kpds     = jpds
     kgds     = jgds

     print*,"- GET GRIB HEADER"
     call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
                 numpts, message_num, kpds, kgds, istat)

     if (istat /= 0) then
       print*,"- FATAL ERROR: BAD DEGRIB OF HEADER. ISTAT IS ", istat
       call w3tage('SNOW2MDL')
       call errexit(61)
     end if

     iafwa = kgds(2)
     jafwa = kgds(3)
     afwa_res = float(kgds(8))*0.001  ! in km.  

     print*,"- DATA VALID AT (YYMMDDHH): ", kpds(8:11)

     print*,"- DEGRIB SNOW DEPTH."

     allocate(bitmap_afwa_sh(iafwa,jafwa))
     allocate(snow_dep_afwa_sh(iafwa,jafwa))

     call getgb(iunit, lugi, (iafwa*jafwa), lskip, jpds, jgds, &
                numpts, lskip, kpds, kgds, bitmap_afwa_sh, snow_dep_afwa_sh, istat)

     if (istat /= 0) then
       print*,"- FATAL ERROR: BAD DEGRIB OF DATA. ISTAT IS ", istat
       call w3tage('SNOW2MDL')
       call errexit(61)
     end if

     kgds_afwa_sh = kgds

     kgds_afwa_sh(7) = -80000  ! ipolates definition of orientation angle is
                               ! 180 degrees off from grib standard.

     call baclose(iunit, istat)

   endif  ! is sh afwa data grib or not?

   call afwa_check(2)

 else

   use_sh_afwa = .false.

 endif

!-------------------------------------------------------------------
!if either hemisphere is bad, don't trust all hemispheres
!-------------------------------------------------------------------

 if (bad_afwa_nh .or. bad_afwa_sh) then
   print*,'- WARNING: AFWA DATA BAD, DO NOT USE.'
   use_nh_afwa = .false.
   use_sh_afwa = .false.
 endif

 return

 end subroutine readafwa 

!>  Check for corrupt nh snow cover data.
!!
!! @note  Check for corrupt nh data by comparing it
!!            to climatology.
!!  
!! program history log:
!! 2009-jun-3   gayno    - initial version
!! 2011-apr-26  gayno    - Perform gross check first,
!!                         then check against climo.
!! 2014-sep-30  gayno    - Weekly climo file converted
!!                         to grib 2.
!!
!! @param[in] kgds_data Grib 1 grid description sect of data to be qcd.
!! @param[in] snow_data Snow cover to be qcd.
!! @param[in] bitmap_data bitmap of data to be qcd.
!! @param[in] idata I dimension of data to be qcd.
!! @param[in] jdata J dimension of data to be qcd.
!! @param[in] isrc Flag indicating data source; 1- afwa depth, 2-ims cover.
!! @param[out] bad When true, data failed check.
!!
!! files:
!!   input:
!!     - NH weekly climatological snow cover file (grib 2).
!!
!! @author  George Gayno org: w/np2 @date  2009-Jun-3
 subroutine nh_climo_check(kgds_data,snow_data,bitmap_data,idata,jdata,isrc,bad)
 use gdswzd_mod

 use program_setup, only    : climo_qc_file,  &
                              grib_year, grib_month, grib_day, &
                              grib_century

 use grib_mod   ! for grib2 library

 implicit none

! describes the climo data grid.
 integer, parameter        :: iclim = 1080
 integer, parameter        :: jclim = 270
 real,  parameter          :: lat11_clim = 90.0
 real,  parameter          :: lon11_clim = -180.0
 real,  parameter          :: dx_clim = 1./3.
 real,  parameter          :: dy_clim = 1./3.

 integer, intent(in)       :: idata, jdata, kgds_data(200), isrc
 logical*1, intent(in)     :: bitmap_data(idata,jdata)
 logical, intent(out)      :: bad
 real, intent(in)          :: snow_data(idata,jdata)

! local variables
 integer                   :: idat(8), jdow, jdoy, jday
 integer                   :: century, year, week, iret, lugb, i, j, ii, jj
 integer                   :: lugi, jdisc, jpdtn, jgdtn, k, nret
 integer                   :: jids(200), jgdt(200), jpdt(200)
 integer                   :: count_nosnow_climo, count_nosnow_data
 integer                   :: count_snow_climo, count_snow_data, count_grosschk_data

 logical*1, allocatable    :: bitmap_clim(:,:)
 logical                   :: unpack

 real, allocatable         :: climo(:,:)
 real                      :: fill, percent, x, y
 real, allocatable         :: xpts(:,:),ypts(:,:),rlon_data(:,:),rlat_data(:,:)
 real                      :: thresh_gross, thresh

 type(gribfield)           :: gfld

 bad=.false.
 if (len_trim(climo_qc_file)==0) return

 print*,"- QC SNOW DATA IN NH."

 if (isrc==1) then
   thresh_gross=50.0   ! afwa data is depth in meters
 elseif (isrc==2) then
   thresh_gross=100.0  ! nesdis/ims data is coverage in percent
 endif

 fill=999.
 allocate(xpts(idata,jdata))
 allocate(ypts(idata,jdata))
 allocate(rlon_data(idata,jdata))
 allocate(rlat_data(idata,jdata))
 do j=1,jdata
 do i=1,idata
   xpts(i,j)=i
   ypts(i,j)=j
 enddo
 enddo

 print*,"- CALC LAT/LONS OF SOURCE POINTS."
 call gdswzd(kgds_data,1,(idata*jdata),fill,xpts,ypts,rlon_data,rlat_data,nret)

 deallocate(xpts,ypts)

 if (nret /= (idata*jdata)) then
   print*,"- WARNING: CALC FAILED. WILL NOT PERFORM QC."
   deallocate (rlon_data,rlat_data)
   return
 endif

 count_grosschk_data=0
 do j=1,jdata
 do i=1,idata
   if (rlat_data(i,j)>0.0 .and. bitmap_data(i,j)) then
     if (snow_data(i,j) < 0.0 .or. snow_data(i,j) > thresh_gross) then
       count_grosschk_data=count_grosschk_data+1
     endif
   endif
 enddo
 enddo

 if (count_grosschk_data > 1) then
   print*,'- NUMBER OF DATA POINTS THAT FAIL GROSS CHECK ',count_grosschk_data
   deallocate (rlon_data,rlat_data)
   bad=.true.
   return
 endif

 print*,"- QC DATA SOURCE AGAINST CLIMO."
 print*,"- OPEN CLIMO SNOW COVER FILE ",trim(climo_qc_file)
 lugb=11
 call baopenr(lugb,climo_qc_file,iret)

 if (iret /= 0) then
   print*,"- WARNING: BAD OPEN, WILL NOT PERFORM QC ", iret
   deallocate (rlon_data,rlat_data)
   return
 endif

!---------------------------------------------------------------
! climo file is weekly. so calculate the current week
! then read that record from the climo file.
!---------------------------------------------------------------

 if (grib_year == 100) then
   century = grib_century
 else
   century = grib_century-1
 endif

 year = century*100 + grib_year

 idat=0
 idat(1)=year
 idat(2)=grib_month
 idat(3)=grib_day

 call w3doxdat(idat,jdow,jdoy,jday)

! the climo file date is the beginning of the 7 day period

 week = nint((jdoy+3.)/7.)
 if (week==0) week=52
 if (week==53) week=1

 print*,"- READ CLIMO FOR WEEK ",week

 call grib2_null(gfld)

 j        = week-1 ! search for specific week (# records to skip)
 lugi     = 0      ! no grib index file
 jdisc    = 0      ! search for discipline; 0 - meteorological products
 jpdtn    = 8      ! search for product definition template number; 8 - average
 jgdtn    = 0      ! search for grid definition template number; 0 - lat/lon grid
 jids     = -9999  ! array of values in identification section, set to wildcard

 jgdt     = -9999  ! array of values in grid definition template 3.m
 jgdt(8)  = iclim  ! search for assumed grid specs - i/j dimensions and corner
                   ! point lat/lons must match.
 jgdt(9)  = jclim
 jgdt(12) = nint(lat11_clim * 1e6)
 jgdt(13) = nint(abs(lon11_clim) * 1e6)

 jpdt     = -9999  ! array of values in product definition template 4.n
 jpdt(1)  = 1      ! search for parameter category - moisture
 jpdt(2)  = 201    ! search for parameter number - snow cover in percent.
 unpack   = .true. ! unpack data

 call getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

 if (iret /= 0) then 
   print*,"- WARNING: PROBLEM READING GRIB FILE ", iret
   print*,"- WILL NOT PERFORM QC."
   deallocate(rlon_data,rlat_data)
   deallocate(climo, bitmap_clim)
   call baclose(lugb,iret)
   return
 endif

 call baclose(lugb,iret)

 allocate(climo(iclim,jclim))
 climo = reshape (gfld%fld , (/iclim,jclim/) )
 allocate(bitmap_clim(iclim,jclim))
 bitmap_clim = reshape (gfld%bmap , (/iclim,jclim/) )

 call grib2_free(gfld)

!---------------------------------------------------------------
! loop over all data points in nh.  gross check data.
! afwa is a depth in meters, ims is % coverage.  there should be 
! no neg values or very large values. if point passes gross check,
! then check against climatology.  find the
! nearest point on the climo snow cover grid.  if
! climo indicates snow is likely (100% coverage), then
! check if afwa/ims has snow.  if climo indicates snow is
! impossible (0% coverage), then check if afwa/ims has no snow.  if
! afwa/ims differs from climo too much, then afwa/ims is
! considered suspect and will not be used.
!---------------------------------------------------------------

 count_nosnow_climo=0
 count_nosnow_data=0
 count_snow_data=0
 count_snow_climo=0

 if (isrc==1) then
   thresh=.005
 elseif (isrc==2) then
   thresh=50.0
 endif

 do j=1,jdata
 do i=1,idata
   if (rlat_data(i,j)>0.0 .and. bitmap_data(i,j)) then
     y = (lat11_clim-rlat_data(i,j))/dy_clim + 1.0
     if (rlon_data(i,j)>180.0) rlon_data(i,j)=rlon_data(i,j)-360.0
     x = (rlon_data(i,j)-lon11_clim)/dx_clim + 1.0
     jj=nint(y)
     if (jj<1) jj=1
     if (jj>jclim) jj=jclim
     ii=nint(x)
     if (ii<1) ii=ii+iclim
     if (ii>iclim) ii=ii-iclim
     if (bitmap_clim(ii,jj)) then  ! climo point is land
       if (climo(ii,jj) <1.0) then ! climo point is snow impossible
         count_nosnow_climo=count_nosnow_climo+1
         if (snow_data(i,j) == 0.0) then
           count_nosnow_data=count_nosnow_data+1
         endif
       endif
       if (climo(ii,jj) > 99.) then  ! climo point is snow likely
         count_snow_climo=count_snow_climo+1
         if (snow_data(i,j) >thresh) then
           count_snow_data=count_snow_data+1
         endif
       endif
     endif
   endif
 enddo
 enddo

 percent = float(count_snow_climo-count_snow_data) / float(count_snow_climo)
 percent = percent*100.
 write(6,200) '- NUMBER OF DATA POINTS THAT SHOULD HAVE SNOW',count_snow_climo
 write(6,201) '- NUMBER OF THESE POINTS THAT ARE BARE GROUND',(count_snow_climo-count_snow_data), &
        'OR', percent, '%'

 200 format(1x,a45,1x,i10)
 201 format(1x,a45,1x,i10,1x,a2,1x,f6.2,a1)

 if (percent>50.0) then
   print*,"- WARNING: PERCENTAGE OF BARE GROUND POINTS EXCEEDS ACCEPTABLE LEVEL."
   print*,"- WILL NOT USE SOURCE DATA." 
   bad=.true.
 endif
   
 percent = float(count_nosnow_climo-count_nosnow_data) / float(count_nosnow_climo)
 percent = percent*100.
 write(6,202) '- NUMBER OF DATA POINTS THAT SHOULD *NOT* HAVE SNOW',count_nosnow_climo
 write(6,203) '- NUMBER OF THESE POINTS WITH SNOW',(count_nosnow_climo-count_nosnow_data), &
        'OR', percent, '%'

 202 format(1x,a51,1x,i10)
 203 format(1x,a34,1x,i10,1x,a2,1x,f6.2,a1)

 if (percent>20.0) then
   print*,"- WARNING: PERCENTAGE OF POINTS WITH SNOW EXCEEDS ACCEPTABLE LEVEL."
   print*,"- WILL NOT USE SOURCE DATA." 
   bad=.true.
 endif

 if (allocated(rlat_data)) deallocate (rlat_data)
 if (allocated(rlon_data)) deallocate (rlon_data)
 if (allocated(climo)) deallocate (climo)
 if (allocated(bitmap_clim)) deallocate (bitmap_clim)

 return

 end subroutine nh_climo_check

!> Check for corrupt afwa data.
!!
!! @param[in] hemi (1-nh, 2-sh)
!! @author  George Gayno  org: w/np2 @date 2009-Jun-3
 subroutine afwa_check(hemi)
  use gdswzd_mod

 implicit none

 integer, intent(in) :: hemi
 integer             :: kgds(200), nret
 integer, parameter  :: npts=1

 real                :: fill, xpts(npts), ypts(npts)
 real                :: rlon(npts), rlat(npts)

 kgds=0
 fill=9999.

 if (hemi==1) then
   print*,'- QC DATA IN NH.'
   kgds=kgds_afwa_nh
   rlat=75.0
   rlon=-40.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_nh(nint(xpts(1)),nint(ypts(1))) < 0.001) then
     print*,'- WARNING: NO SNOW IN GREENLAND: ',snow_dep_afwa_nh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_nh=.true.
   endif
   rlat=3.0
   rlon=-60.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_nh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN S AMERICA: ',snow_dep_afwa_nh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_nh=.true.
   endif
   rlat=23.0
   rlon=10.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_nh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN SAHARA: ',snow_dep_afwa_nh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_nh=.true.
   endif
   rlat=15.0
   rlon=10.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_nh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN S INDIA: ',snow_dep_afwa_nh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_nh=.true.
   endif
 endif

 if (hemi==2) then
   print*,'- QC DATA IN SH.'
   kgds=kgds_afwa_sh
   rlat=-88.0
   rlon=0.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_sh(nint(xpts(1)),nint(ypts(1))) < 0.001) then
     print*,'- WARNING: NO SNOW IN ANTARCTICA: ',snow_dep_afwa_sh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_sh=.true.
   endif
   rlat=-10.
   rlon=-45.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_sh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN SOUTH AMERICA: ',snow_dep_afwa_sh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_sh=.true.
   endif
   rlat=-20.0
   rlon=130.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_sh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN AUSTRALIA: ',snow_dep_afwa_sh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_sh=.true.
   endif
   rlat=-9.0
   rlon=25.
   call gdswzd(kgds,(-1),npts,fill,xpts,ypts,rlon,rlat,nret)
   if (snow_dep_afwa_sh(nint(xpts(1)),nint(ypts(1))) > 0.0) then
     print*,'- WARNING: SNOW IN AFRICA: ',snow_dep_afwa_sh(nint(xpts),nint(ypts))
     print*,'- DONT USE AFWA DATA.'
     bad_afwa_sh=.true.
   endif
 endif

 end subroutine afwa_check

!> Read afwa binary snow depth file.
!!
!! @param[in] file_name file name
!! @param[out] snow_dep_afwa snow depth in meters
!!
!! files:
!!   input:
!!     - nh/sh afwa data in simple binary format
!!
!! condition codes: all fatal
!!   60 - bad open of afwa file
!!   61 - bad read of afwa file
!!
!! @note Read logic for binary data is  taken from hua-lu's code,
!! /nwprod/sorc/grib_snowgrib.fd.
!!
!! @author  George Gayno org: w/np2 @date  2007-Nov-28
 subroutine read_afwa_binary(file_name, snow_dep_afwa) 

 implicit none

 character*8                            :: afwa_file_info(2)
 character*(*), intent(in)              :: file_name

 integer*2, allocatable                 :: dummy(:,:)
 integer                                :: i,j, istat
 integer, parameter                     :: iafwa = 512
 integer, parameter                     :: jafwa = 512
 integer, parameter                     :: iunit=11  ! input afwa data file

 real, intent(out)                      :: snow_dep_afwa(iafwa,jafwa)

 print*,"- OPEN AFWA BINARY FILE ", trim(file_name)
 open (iunit, file=trim(file_name), access="direct", recl=iafwa*2, iostat=istat)

 if (istat /= 0) then
   print*,'- FATAL ERROR: BAD OPEN.  ISTAT IS ',istat
   call w3tage('SNOW2MDL')
   call errexit(60)
 end if

 print*,"- READ AFWA BINARY FILE ", trim(file_name)
 read(iunit, rec=2, iostat = istat) afwa_file_info

 if (istat /= 0) then
   print*,'- FATAL ERROR: BAD READ.  ISTAT IS ',istat
   call w3tage('SNOW2MDL')
   call errexit(61)
 end if

 print*,"- AFWA DATA IS ", afwa_file_info(1), " AT TIME ", afwa_file_info(2)(2:7)

 allocate(dummy(iafwa,jafwa))
 
 do j = 1, jafwa
   read(iunit, rec=j+2, iostat=istat) (dummy(i,j),i=1,iafwa)
   if (istat /= 0) then
     print*,'- FATAL ERROR: BAD READ.  ISTAT IS ',istat
     call w3tage('SNOW2MDL')
     call errexit(61)
   end if
 enddo

 close(iunit)

!-----------------------------------------------------------------------
! "4090" is the sea ice flag.  we don't use the afwa sea ice.
!-----------------------------------------------------------------------

 where (dummy == 4090) dummy = 0

 snow_dep_afwa = float(dummy)  

!---------------------------------------------------------------------
! afwa data is a snow depth in units of tenths of inches.
! convert this to meters.
!---------------------------------------------------------------------

 snow_dep_afwa = snow_dep_afwa * 2.54 / 1000.0

 deallocate (dummy)

 return

 end subroutine read_afwa_binary

!> Read afwa land mask file to get a bitmap.
!!  
!! @param[in]   file_name      land mask file name
!! @param[out]  bitmap_afwa   .true. if land
!!
!! files:
!!   input:
!!    - afwa landmask in simple binary format
!!
!! condition codes: all fatal
!!   62 - bad open of afwa landmask file
!!   63 - bad read of afwa landmask file
!!
!! @author  George Gayno org: w/np2 @date 2007-Nov-28
 subroutine read_afwa_mask(file_name, bitmap_afwa)
 implicit none

 character*(*), intent(in)         :: file_name

 integer, parameter                :: iunit=11 ! input mask file
 integer, parameter                :: iafwa = 512
 integer, parameter                :: jafwa = 512
 integer                           :: i, j, istat
 integer*4, allocatable            :: dummy4(:,:)

 logical*1, intent(out)            :: bitmap_afwa(iafwa,jafwa)

 allocate (dummy4(iafwa,jafwa))

 print*,'- OPEN AFWA MASK FILE ', trim(file_name)
 open(iunit, file=trim(file_name), access='direct', &
      recl=iafwa*jafwa*4, iostat=istat)

 if (istat /= 0) then
   print*,'- FATAL ERROR: BAD OPEN. ISTAT IS ', istat
   call w3tage('SNOW2MDL')
   call errexit(62)
 end if

 print*,'- READ AFWA MASK FILE ', trim(file_name)
 read(iunit, rec=1, iostat=istat) dummy4

 if (istat /= 0) then
   print*,'- FATAL ERROR: BAD READ. ISTAT IS ', istat
   call w3tage('SNOW2MDL')
   call errexit(63)
 end if
 
 close(iunit)

!-----------------------------------------------------------------------
! here -1-offhemi, 1-ocean, 2-land, 4-coast.
!-----------------------------------------------------------------------

 bitmap_afwa = .false.

 do j = 1, jafwa
   do i = 1, iafwa
     if (dummy4(i,j) > 1) then
       bitmap_afwa(i,j) = .true.
     endif
   enddo
 enddo

 deallocate (dummy4)

 end subroutine read_afwa_mask


 end module snowdat
