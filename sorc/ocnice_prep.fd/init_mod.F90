!> @file
!! @brief Define the input namelist variables
!! @author Denise.Worthen@noaa.gov
!!
!> This module contains the namelist variables
!! @author Denise.Worthen@noaa.gov
module init_mod

  implicit none

  public

  integer, parameter :: maxvars = 60         !< The maximum number of fields expected in a source file
  character(len=10)  :: maskvar = 'h'        !< The variable in the ocean source file used to create
                                             !< the interpolation mask with dynamic masking
  type :: vardefs
     character(len= 20)   :: var_name        !< A variable's variable name
     character(len=120)   :: long_name       !< A variable's long name
     character(len= 20)   :: units           !< A variable's unit
     character(len= 20)   :: var_remapmethod !< A variable's mapping method
     integer              :: var_dimen       !< A variable's dimensionality
     character(len=  4)   :: var_grid        !< A variable's input grid location
     character(len= 20)   :: var_pair        !< A variable's pair
     character(len=  4)   :: var_pair_grid   !< A pair variable grid
  end type vardefs

  type(vardefs) :: outvars(maxvars)          !< An empty structure filled by reading a csv file
                                             !< describing the fields

  character(len=10)  :: ftype      !< The type of tripole grid (ocean or ice)
  character(len=10)  :: fsrc       !< A character string for tripole grid
  character(len=10)  :: fdst       !< A character string for the destination grid
  character(len=120) :: wgtsdir    !< The directory containing the regridding weights
  character(len=120) :: griddir    !< The directory containing the master tripole grid file
  character(len=20)  :: input_file !< The input file name

  integer :: nxt        !< The x-dimension of the source tripole grid
  integer :: nyt        !< The y-dimension of the source tripole grid
  integer :: nlevs      !< The vertical or category dimension of the source tripole grid

  integer :: nxr        !< The x-dimension of the destination tripole grid
  integer :: nyr        !< The y-dimension of the destination tripole grid

  integer :: logunit    !< The log unit
  logical :: debug      !< If true, print debug messages and intermediate files
  logical :: do_ocnprep !< If true, the source file is ocean, otherwise ice

contains
  !>  Read input namelist file
  !!
  !! @param[in]   fname     namelist file
  !! @param[out]  errmsg    return error message
  !! @param[out]  rc        return error code
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine readnml(fname,errmsg,rc)

    character(len=*), intent(in)  :: fname
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: rc

    ! local variable
    logical :: fexist
    integer :: ierr, iounit
    integer :: srcdims(2), dstdims(2)
    !----------------------------------------------------------------------------

    namelist /ocniceprep_nml/ ftype, wgtsdir, griddir, srcdims, dstdims, debug

    srcdims = 0; dstdims = 0
    errmsg='' ! for successful return
    rc = 0    ! for successful retun

    inquire(file=trim(fname), exist=fexist)
    if (.not. fexist) then
       write (errmsg, '(a)') 'FATAL ERROR: input file '//trim(fname)//' does not exist.'
       rc = 1
       return
    else
       ! Open and read namelist file.
       open (action='read', file=trim(fname), iostat=ierr, newunit=iounit)
       read (nml=ocniceprep_nml, iostat=ierr, unit=iounit)
       if (ierr /= 0) then
          rc = 1
          write (errmsg, '(a)') 'FATAL ERROR: invalid namelist format.'
          return
       end if
       close (iounit)
    end if

    ! check that model is either ocean or ice
    if (trim(ftype) /= 'ocean' .and. trim(ftype) /= 'ice') then
       rc = 1
       write (errmsg, '(a)') 'FATAL ERROR: ftype must be ocean or ice'
       return
    end if

    ! set grid dimensions and names
    nxt = srcdims(1); nyt = srcdims(2)
    nxr = dstdims(1); nyr = dstdims(2)
    fsrc = '' ; fdst = ''
    if (nxt == 1440 .and. nyt == 1080) fsrc = 'mx025'    ! 1/4deg tripole
    if (len_trim(fsrc) == 0) then
       rc = 1
       write(errmsg,'(a)')'FATAL ERROR: source grid dimensions incorrect'
       return
    end if

    if (nxr == 720  .and. nyr == 576) fdst = 'mx050'     ! 1/2deg tripole
    if (nxr == 360  .and. nyr == 320) fdst = 'mx100'     ! 1deg tripole
    if (nxr == 72   .and. nyr == 35)  fdst = 'mx500'     ! 5deg tripole
    if (len_trim(fdst) == 0) then
       rc = 1
       write(errmsg,'(a)')'FATAL ERROR: destination grid dimensions incorrect'
       return
    end if

    ! initialize the source file types
    if (trim(ftype) == 'ocean') then
       do_ocnprep = .true.
    else
       do_ocnprep = .false.
    end if

    input_file = trim(ftype)//'.nc'
    inquire (file=trim(input_file), exist=fexist)
    if (.not. fexist) then
       write (errmsg, '(a)') 'FATAL ERROR: input file '//trim(input_file)//' does not exist.'
       rc=1
       return
    end if

    ! log file
    open(newunit=logunit, file=trim(ftype)//'.prep.log',form='formatted')
    if (debug) write(logunit, '(a)')'input file: '//trim(input_file)

    ! all checks pass, continue
    write(errmsg,'(a)') 'Namelist successfully read, continue'
    rc = 0

  end subroutine readnml

  !> Read the input csv file and fill the vardefs type
  !!
  !! @param[in]    fname     input csv file
  !! @param[out]   errmsg    return error message
  !! @param[out]   rc        return error code
  !! @param[out]  nvalid    the number of variables in the csv file
  !!
  !! @author Denise.Worthen@noaa.gov
  subroutine readcsv(fname,errmsg,rc,nvalid)

    character(len=*), intent(in)  :: fname
    character(len=*), intent(out) :: errmsg
    integer, intent(out)          :: rc
    integer, intent(out)          :: nvalid

    ! local variables
    character(len=100) :: chead
    character(len= 20) :: c1,c3,c4,c5,c6
    integer :: i2, idx1,idx2
    integer :: nn,n,ierr,iounit
    !----------------------------------------------------------------------------

    errmsg='' ! for successful return
    rc = 0    ! for successful retun

    open(newunit=iounit, file=trim(fname), status='old', iostat=ierr)
    if (ierr /= 0) then
       rc = 1
       write (errmsg, '(a)') 'FATAL ERROR: input file '//trim(fname)//' does not exist.'
       return
    end if

    read(iounit,*)chead
    nn=0
    do n = 1,maxvars
       read(iounit,*,iostat=ierr)c1,i2,c3,c4,c5,c6
       if (ierr .ne. 0) exit
       if (len_trim(c1) > 0) then
          nn = nn+1
          outvars(nn)%var_name = trim(c1)
          outvars(nn)%var_dimen = i2
          outvars(nn)%var_grid = trim(c3)
          outvars(nn)%var_remapmethod = trim(c4)
          outvars(nn)%var_pair = trim(c5)
          outvars(nn)%var_pair_grid = trim(c6)
       end if
    end do
    close(iounit)
    nvalid = nn

    ! check for u,v pairs, these should be listed in csv file in ordered pairs
    idx1 = 0; idx2 = 0
    do n = 1,nvalid
       if (len_trim(outvars(n)%var_pair) > 0 .and. idx1 .eq. 0) then
          idx1 = n
          idx2 = n+1
       end if
    end do

    if (trim(outvars(idx1)%var_pair) /= trim(outvars(idx2)%var_name)) then
       rc = 1
       write(errmsg,'(a)')'FATAL ERROR: vector pair for '//trim(outvars(idx1)%var_name)//' is not set correctly'
       return
    end if
    if (trim(outvars(idx2)%var_pair) /= trim(outvars(idx1)%var_name)) then
       rc = 1
       write(errmsg,'(a)')'FATAL ERROR: vector pair for '//trim(outvars(idx2)%var_name)//' is not set correctly'
       return
    end if

    ! check for u velocities on u-staggers and v-velocities on v-staggers
    if (outvars(idx1)%var_name(1:1) == 'u') then
       if ((outvars(idx1)%var_grid(1:2) /= 'Cu') .and. outvars(idx1)%var_grid(1:2) /= 'Bu') then
          rc = 1
          write(errmsg,'(a)')'FATAL ERROR: u-vector has wrong grid '
          return
       end if
    end if
    if (outvars(idx2)%var_name(1:1) == 'v') then
       if ((outvars(idx2)%var_grid(1:2) /= 'Cv') .and. outvars(idx2)%var_grid(1:2) /= 'Bu') then
          rc = 1
          write(errmsg,'(a)')'FATAL ERROR: v-vector has wrong grid '
          return
       end if
    end if

    ! all checks pass, continue
    write(errmsg,'(a)')'CSV successfully read, continue'
    rc = 0

  end subroutine readcsv
end module init_mod
