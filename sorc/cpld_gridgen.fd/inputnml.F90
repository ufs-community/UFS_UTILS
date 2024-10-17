!> @file
!! @brief Define the input namelist variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module contains the namelist variables
!! @author Denise.Worthen@noaa.gov

module inputnml

  use grdvars,     only : nx,ny,ni,nj,npx
  use grdvars,     only : editmask, debug, do_postwgts
  use charstrings, only : dirsrc, dirout, fv3dir, res, topofile, editsfile

  implicit none

contains

  !>  Read input namelist file
  !!
  !! @param[in]  fname  the file name to read
  !!
  !! @author Denise.Worthen@noaa.gov

  subroutine read_inputnml(fname)

    character(len=*),   intent(in) :: fname

    ! local variables
    integer :: iounit, rc
    character(len=200) :: tmpstr

    namelist /grid_nml/ ni, nj, dirsrc, dirout, fv3dir,  topofile, editsfile, &
         res, editmask, debug, do_postwgts

    ! Check whether file exists.
    inquire (file=trim(fname), iostat=rc)
    if (rc /= 0) then
       write (0, '(3a)') 'Error: input file "', trim(fname), '" does not exist.'
       stop 1
    end if

    ! Open and read Namelist file.
    open (action='read', file=trim(fname), iostat=rc, newunit=iounit)
    read (nml=grid_nml, iostat=rc, unit=iounit)
    if (rc /= 0) then
       backspace(iounit)
       read(iounit,'(a)')tmpstr
       write (6, '(a)') 'Error: invalid Namelist format '//trim(tmpstr)
       stop 1
    end if
    close(iounit)

    ! set supergrid dimensions
    nx = ni*2
    ny = nj*2

  end subroutine read_inputnml
end module inputnml
