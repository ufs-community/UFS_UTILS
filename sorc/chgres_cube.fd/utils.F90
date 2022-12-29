module utilities

contains
!> @file
!! @brief Contains utility routines.
!!

!> General error handler.
!!
!! @param[in] string  error message
!! @param[in] rc      error status code
 subroutine error_handler(string, rc)

 use mpi

 implicit none

 character(len=*), intent(in)    :: string
 
 integer,          intent(in)    :: rc

 integer :: ierr

 print*,"- FATAL ERROR: ", string
 print*,"- IOSTAT IS: ", rc
 call mpi_abort(mpi_comm_world, 999, ierr)

 end subroutine error_handler

!> Error handler for netcdf
!!
!! @param[in] err     error status code
!! @param[in] string  error message
 subroutine netcdf_err( err, string )

 use mpi
 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg
 integer :: iret

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 call mpi_abort(mpi_comm_world, 999, iret)

 return
 end subroutine netcdf_err
 
!> Convert string from lower to uppercase.
!! @author Clive Page
!!
!! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
!!
!! @param[in] strIn   string to convert
!! @return strOut string in uppercase
function to_upper(strIn) result(strOut)

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

!> Convert from upper to lowercase
!! @author Clive Page
!!
!! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
!!
!! @param[in,out] strIn   string to convert
subroutine to_lower(strIn)

     implicit none

     character(len=*), intent(inout) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do
  strIn(:) = strOut(:)
end subroutine to_lower

!> Handle GRIB2 read error based on the user selected
!! method in the varmap file.
!!
!! @param [in] vname  grib2 variable name
!! @param [in] lev    grib2 variable level
!! @param [in] method  how missing data is handled
!! @param [in] value   fill value for missing data
!! @param [in] varnum  grib2 variable number
!! @param [inout] iret  return status code
!! @param [inout] var   4-byte array of corrected data
!! @param [inout] var8  8-byte array of corrected data
!! @param [inout] var3d 3-d array of corrected data
!! @param [inout] read_from_input  logical array indicating if variable was read in 
!! @author Larissa Reames
subroutine handle_grib_error(vname,lev,method,value,varnum,read_from_input, iret,var,var8,var3d)

  use, intrinsic :: ieee_arithmetic
  use esmf

  implicit none
   
  real(esmf_kind_r4), intent(in)    :: value
  logical, intent(inout)            :: read_from_input(:)
  real(esmf_kind_r4), intent(inout), optional :: var(:,:)
  real(esmf_kind_r8), intent(inout), optional :: var8(:,:)
  real(esmf_kind_r8), intent(inout), optional  :: var3d(:,:,:)
  
  character(len=20), intent(in)     :: vname, lev, method
  
  integer, intent(in)               :: varnum 
  integer, intent(inout)            :: iret
  
  iret = 0
  if (varnum == 9999) then
    print*, "WARNING: ", trim(vname), " NOT FOUND AT LEVEL ", lev, " IN EXTERNAL FILE ", &
            "AND NO ENTRY EXISTS IN VARMAP TABLE. VARIABLE WILL NOT BE USED."
    iret = 1

    return
  endif

  if (trim(method) == "skip" ) then
    print*, "WARNING: SKIPPING ", trim(vname), " IN FILE"
    read_from_input(varnum) = .false.
    iret = 1
  elseif (trim(method) == "set_to_fill") then
    print*, "WARNING: ,", trim(vname), " NOT AVAILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO FILL VALUE OF ", value
    if(present(var)) var(:,:) = value
    if(present(var8)) var8(:,:) = value
    if(present(var3d)) var3d(:,:,:) = value
  elseif (trim(method) == "set_to_NaN") then
    print*, "WARNING: ,", trim(vname), " NOT AVAILABLE AT LEVEL ", trim(lev), &
           ". SETTING EQUAL TO NaNs"
    if(present(var)) var(:,:) = ieee_value(var,IEEE_QUIET_NAN)
    if(present(var8)) var8(:,:) = ieee_value(var8,IEEE_QUIET_NAN)
    if(present(var3d)) var3d(:,:,:) = ieee_value(var3d,IEEE_QUIET_NAN)
  elseif (trim(method) == "stop") then
    call error_handler("READING "//trim(vname)// " at level "//lev//". TO MAKE THIS NON- &
                        FATAL, CHANGE STOP TO SKIP FOR THIS VARIABLE IN YOUR VARMAP &
                        FILE.", iret)
  elseif (trim(method) == "intrp") then
    print*, "WARNING: ,"//trim(vname)//" NOT AVAILABLE AT LEVEL "//trim(lev)// &
          ". WILL INTERPOLATE INTERSPERSED MISSING LEVELS AND/OR FILL MISSING"//&
          " LEVELS AT EDGES."
  else
    call error_handler("ERROR USING MISSING_VAR_METHOD. PLEASE SET VALUES IN" // &
                       " VARMAP TABLE TO ONE OF: set_to_fill, set_to_NaN,"// &
                       " , intrp, skip, or stop.", 1)
  endif

end subroutine handle_grib_error

!> Sort an array of values.
!!
!! @param a      the sorted array
!! @param first  the first value of sorted array
!! @param last   the last value of sorted array
!! @author Jili Dong NOAA/EMC
recursive subroutine quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

!> Check for and replace certain values in soil temperature.
!> At open water points (landmask=0) use skin temperature as
!> a filler value. At land points (landmask=1) with excessive
!> soil temperature, replace soil temperature with skin temperature. 
!> In GEFSv12.0 data there are some erroneous missing values at
!> land points which this corrects. At sea ice points (landmask=2),
!> store a default ice column temperature because grib2 files do not 
!> have ice column temperature which FV3 expects at these points.
!!
!! @param soilt    [inout] 3-dimensional soil temperature arrray
!! @param landmask [in]    landmask of the input grid
!! @param skint    [in]    2-dimensional skin temperature array
!! @param ICET_DEFAULT [in] Default temperature to apply at ice points
!! @param i_input  [in]    i-dimension of input grid
!! @param j_input  [in]    j-dimension of input grid
!! @param lsoil_input  [in]  soil layers dimension of input grid 
!! @author Larissa Reames CIMMS/NSSL

subroutine check_soilt(soilt, landmask, skint,ICET_DEFAULT,i_input,j_input,lsoil_input)
  use esmf
  implicit none
  integer, intent(in)               :: i_input, j_input, lsoil_input
  real(esmf_kind_r8), intent(inout) ::  soilt(i_input,j_input,lsoil_input)
  real(esmf_kind_r8), intent(in)    ::  skint(i_input,j_input)
  real,  intent(in) :: ICET_DEFAULT
  integer(esmf_kind_i4), intent(in)    ::  landmask(i_input,j_input)
  
  integer                           :: i, j, k

  do k=1,lsoil_input
    do j = 1, j_input
      do i = 1, i_input
        if (landmask(i,j) == 0_esmf_kind_i4 ) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 1_esmf_kind_i4 .and. soilt(i,j,k) > 350.0_esmf_kind_r8) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 2_esmf_kind_i4 ) then 
          soilt(i,j,k) = ICET_DEFAULT
        endif
      enddo
    enddo
  enddo
end subroutine check_soilt

!> When using GEFS data, some points on the target grid have 
!> unreasonable canpy moisture content, so zero out any 
!> locations with unrealistic canopy moisture values (>0.5).
!!
!! @param cnwat [input] 2-dimensional canopy moisture content
!! @param i_input [in]  i-dimension of input grid
!! @param j_input  [in] j-dimension of input grid
!! @author Larissa Reames CIMMS/NSSL

subroutine check_cnwat(cnwat,i_input,j_input)
  use esmf
  implicit none 
  integer, intent(in)               :: i_input, j_input
  real(esmf_kind_r8), intent(inout) :: cnwat(i_input,j_input)
  
  real(esmf_kind_r8)                :: max_cnwat = 0.5
  
  integer :: i, j

  do i = 1,i_input
    do j = 1,j_input
      if (cnwat(i,j) .gt. max_cnwat) cnwat(i,j) = 0.0_esmf_kind_r8
    enddo
  enddo
end subroutine check_cnwat

!> Pressure to presure vertical interpolation for tracers with linear or lnP
!> interpolation. Input tracers on pres levels are interpolated 
!> to the target output pressure levels. The matching levels of input and 
!> output will keep the same. Extrapolation is also allowed but needs
!> caution. The routine is mostly for GFSV16 combined grib2 input when spfh has
!> missing levels in low and mid troposphere from U/T/HGT/DZDT. 
!!
!! @param [in] ppin  1d input pres levs
!! @param [in] xxin  1d input tracer
!! @param [in] npin  number of input levs 
!! @param [in] ppout 1d target pres levs 
!! @param [out] xxout 1d interpolated tracer
!! @param [in] npout number of target levs 
!! @param [in] linlog interp method.1:linear;not 1:log;neg:extrp allowed
!! @param [in] xmsg  fill values of missing levels (-999.0)
!! @param [out] ier  error status. non 0: failed interpolation
!! @author Jili Dong NCEP/EMC  
!! @date 2021/07/30
SUBROUTINE DINT2P(PPIN,XXIN,NPIN,PPOUT,XXOUT,NPOUT   &                   
   ,LINLOG,XMSG,IER)
      IMPLICIT NONE

! NCL code for pressure level interpolation
!
! This code was designed for one simple task. It has since
! been mangled and abused for assorted reasons. For example,
! early gfortran compilers had some issues with automatic arrays.
! Hence, the C-Wrapper was used to create 'work' arrays which
! were then passed to this code.  The original focused (non-NCL) 
! task was to handle PPIN & PPOUT that had the same 'monotonicity.' 
! Extra code was added to handle the more general case. 
! Blah-Blah:  Punch line: it is embarrassingly convoluted!!!
!
!                                                ! input types
      INTEGER NPIN,NPOUT,LINLOG,IER
      real*8 PPIN(NPIN),XXIN(NPIN),PPOUT(NPOUT),XMSG
                                                ! output
      real*8 XXOUT(NPOUT)
                                                ! work
      real*8 PIN(NPIN),XIN(NPIN),P(NPIN),X(NPIN)
      real*8 POUT(NPOUT),XOUT(NPOUT)

! local
      INTEGER J1,NP,NL,NIN,NLMAX,NPLVL,NLSAVE,NP1,NO1,N1,N2,LOGLIN,   &
             NLSTRT
      real*8 SLOPE,PA,PB,PC

      LOGLIN = ABS(LINLOG)

! error check: enough points: pressures consistency?

      IER = 0
      IF (NPOUT.GT.0) THEN
          DO NP = 1,NPOUT
              XXOUT(NP) = XMSG
          END DO
      END IF
! Jili Dong input levels have to be the same as output levels:
! we only interpolate for levels with missing variables
!      IF (.not. all(PPIN .eq. PPOUT)) IER = IER+1

      IF (NPIN.LT.2 .OR. NPOUT.LT.1) IER = IER + 1

      IF (IER.NE.0) THEN
!          PRINT *,'INT2P: error exit: ier=',IER
          RETURN
      END IF

! should *input arrays* be reordered? want p(1) > p(2) > p(3) etc
! so that it will match order for which code was originally designed
! copy to 'work'  arrays

      NP1 = 0
      NO1 = 0
      IF (PPIN(1).LT.PPIN(2)) THEN
          NP1 = NPIN + 1
      END IF
      IF (PPOUT(1).LT.PPOUT(2)) THEN
          NO1 = NPOUT + 1
      END IF

      DO NP = 1,NPIN
          PIN(NP) = PPIN(ABS(NP1-NP))
          XIN(NP) = XXIN(ABS(NP1-NP))
      END DO

      DO NP = 1,NPOUT
          POUT(NP) = PPOUT(ABS(NO1-NP))
      END DO

! eliminate XIN levels with missing data. 
! .   This can happen with observational data.

      NL = 0
      DO NP = 1,NPIN
          IF (XIN(NP).NE.XMSG .AND. PIN(NP).NE.XMSG) THEN
              NL = NL + 1
              P(NL) = PIN(NP)
              X(NL) = XIN(NP)
          END IF
      END DO
      NLMAX = NL

                                                ! all missing data
      IF (NLMAX.LT.2) THEN
          IER = IER + 1000
          PRINT *,'INT2P: ier=',IER
          RETURN
      END IF

! ===============> pressure in decreasing order <================
! perform the interpolation  [pin(1)>pin(2)>...>pin(npin)]
!                                                      ( p ,x)
! ------------------------- p(nl+1), x(nl+1)   example (200,5)
! .
! ------------------------- pout(np), xout(np)         (250,?)
! .
! ------------------------- p(nl)  , x(nl)             (300,10)


! exact p-level matches
      NLSTRT = 1
      NLSAVE = 1
      DO NP = 1,NPOUT
          XOUT(NP) = XMSG
          DO NL = NLSTRT,NLMAX
              IF (POUT(NP).EQ.P(NL)) THEN
                  XOUT(NP) = X(NL)
                  NLSAVE = NL + 1
                  GO TO 10
              END IF
          END DO
   10     NLSTRT = NLSAVE
      END DO

      IF (LOGLIN.EQ.1) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      SLOPE = (X(NL)-X(NL+1))/ (P(NL)-P(NL+1))
                      XOUT(NP) = X(NL+1) + SLOPE* (POUT(NP)-P(NL+1))
                  END IF
              END DO
          END DO
      ELSE
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      PA = LOG(P(NL))
                      PB = LOG(POUT(NP))
! special case: In case someome inadvertently enter p=0.
                      if (p(nl+1).gt.0.d0) then
                          PC = LOG(P(NL+1))
                      else
                          PC = LOG(1.d-4)
                      end if

                      SLOPE = (X(NL)-X(NL+1))/ (PA-PC)
                      XOUT(NP) = X(NL+1) + SLOPE* (PB-PC)
                  END IF
              END DO
          END DO
      END IF

! extrapolate?
! . use the 'last' valid slope for extrapolating

      IF (LINLOG.LT.0) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX
                  IF (POUT(NP).GT.P(1)) THEN
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(2)-X(1))/ (P(2)-P(1))
                          XOUT(NP) = X(1) + SLOPE* (POUT(NP)-P(1))
                      ELSE
                          PA = LOG(P(2))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(1))
                          SLOPE = (X(2)-X(1))/ (PA-PC)
                          XOUT(NP) = X(1) + SLOPE* (PB-PC)
                      END IF
                  ELSE IF (POUT(NP).LT.P(NLMAX)) THEN
                      N1 = NLMAX
                      N2 = NLMAX - 1
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(N1)-X(N2))/ (P(N1)-P(N2))
                          XOUT(NP) = X(N1) + SLOPE* (POUT(NP)-P(N1))
                      ELSE
                          PA = LOG(P(N1))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(N2))
                          SLOPE = (X(N1)-X(N2))/ (PA-PC)
                          !XOUT(NP) = X(N1) + SLOPE* (PB-PC) !bug fixed below
                          XOUT(NP) = X(N1) + SLOPE* (PB-PA)
                      END IF
                  END IF
              END DO
          END DO
      END IF

! place results in the return array;
! .   possibly .... reverse to original order

      if (NO1.GT.0) THEN
          DO NP = 1,NPOUT
             n1 = ABS(NO1-NP)
             PPOUT(NP) = POUT(n1)
             XXOUT(NP) = XOUT(n1)
          END DO
      ELSE
          DO NP = 1,NPOUT
             PPOUT(NP) = POUT(NP)
             XXOUT(NP) = XOUT(NP)
          END DO
      END IF


      RETURN
      END SUBROUTINE DINT2P
end module utilities
