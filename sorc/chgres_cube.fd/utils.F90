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

subroutine handle_grib_error(vname,lev,method,value,varnum, iret,var,var8,var3d)

  use, intrinsic :: ieee_arithmetic

  implicit none
  
  real(esmf_kind_r4), intent(in)    :: value
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
!! @author Larissa Reames CIMMS/NSSL

subroutine check_soilt(soilt, landmask, skint,ICET_DEFAULT)
  implicit none
  real(esmf_kind_r8), intent(inout) ::  soilt(i_input,j_input,lsoil_input)
  real(esmf_kind_r8), intent(in)    ::  skint(i_input,j_input)
  real(esmf_kind_r8), parameter, intent(in) :: ICET_DEFAULT
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
!! @author Larissa Reames CIMMS/NSSL

subroutine check_cnwat(cnwat)
  implicit none 
  real(esmf_kind_r8), intent(inout) :: cnwat(i_input,j_input)
  
  real(esmf_kind_r8)                :: max_cnwat = 0.5
  
  integer :: i, j

  do i = 1,i_input
    do j = 1,j_input
      if (cnwat(i,j) .gt. max_cnwat) cnwat(i,j) = 0.0_esmf_kind_r8
    enddo
  enddo
end subroutine check_cnwat

