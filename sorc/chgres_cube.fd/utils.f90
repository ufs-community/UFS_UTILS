 subroutine error_handler(string, rc)

 implicit none

 character(len=*), intent(in)    :: string
 
 integer,          intent(in)    :: rc

 print*,"- FATAL ERROR: ", string
 print*,"- IOSTAT IS: ", rc
 call mpi_abort

 end subroutine error_handler

 subroutine netcdf_err( err, string )

 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg

 include "mpif.h"

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 call mpi_abort(mpi_comm_world, 999)

 return
 end subroutine netcdf_err
 
 subroutine to_upper(strIn)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(inout) :: strIn
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
  strIn(:) = strOut(:)
end subroutine to_upper

subroutine to_lower(strIn)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

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
