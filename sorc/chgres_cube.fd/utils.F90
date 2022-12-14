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

 print*,"- FATAL ERROR: ", trim(string)
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
