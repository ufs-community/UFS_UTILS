!=======================================================================
program global_equiv_resol
!=======================================================================

  use netcdf

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi_geom = 4.0*atan(1.0), &
                         radius_Earth = 6371000.0

  character(len=256) :: grid_fn
  integer :: ncid, nxSG_dimid, nySG_dimid, dASG_varid, num_args
  integer :: nxSG, nySG, nx, ny, RES_equiv
  real(dp) :: avg_cell_size, min_cell_size, max_cell_size
  real(dp), dimension(:,:), allocatable :: &
    quarter_dA_ll, quarter_dA_lr, quarter_dA_ur, quarter_dA_ul, &
    dASG, dA, sqrt_dA
!
!=======================================================================
!
! Read in the name of the file from the command line.  The command-line
! call to this program should have exactly one argument consisting of 
! the path to the NetCDF grid specification file to be read in.  If this
! is not the case, print out a usage message and exit.
!
!=======================================================================
!
  num_args = command_argument_count()
  if (num_args == 1) then
    call get_command_argument(1, grid_fn)
  else
    WRITE(*,500)
    WRITE(*,500) "Exactly one argument must be specified to program global_equiv_resol."
    WRITE(*,500) "Usage:"
    WRITE(*,500)
    WRITE(*,500) "  global_equiv_resol  path_to_grid_file"
    WRITE(*,500)
    WRITE(*,500) "where path_to_grid_file is the path to the NetCDF grid file.  Actual "
    WRITE(*,500) "number of specified command line arguments is:"
    WRITE(*,510) "  num_args = ", num_args
    WRITE(*,500) "Stopping."
500 FORMAT(A)
510 FORMAT(A, I3)
    STOP
  end if
!
!=======================================================================
!
! Open the grid file and read in the dimensions of the supergrid.  The
! supergrid is a grid that has twice the resolution of the actual/compu-
! tational grid.  In the file, the names of the supergrid dimensions are
! nx and ny.  Here, however, we reserve those names for the dimensions 
! of the actual grid (since in the FV3 code and in other data files, nx
! and ny are used to denote the dimensions of the actual grid) and in-
! stead use the variables nxSG and nySG to denote the dimensions of the
! supergrid.
!
!=======================================================================
!
  WRITE(*,500)
  WRITE(*,500) "Opening NetCDF grid file for reading/writing:"
  WRITE(*,500) "  grid_fn = " // trim(grid_fn)

  call check( nf90_open(trim(grid_fn), NF90_WRITE, ncid) )

  call check( nf90_inq_dimid(ncid, "nx", nxSG_dimid) )
  call check( nf90_inquire_dimension(ncid, nxSG_dimid, len=nxSG) )

  call check( nf90_inq_dimid(ncid, "ny", nySG_dimid) )
  call check( nf90_inquire_dimension(ncid, nySG_dimid, len=nySG) )

  WRITE(*,500)
  WRITE(*,500) "Dimensions of supergrid are:"
  WRITE(*,520) "  nxSG = ", nxSG
  WRITE(*,520) "  nySG = ", nySG
520 FORMAT(A, I7)
!
!=======================================================================
!
! Read in the cell areas on the supergrid.  Then add the areas of the 
! four supergrid cells that make up one grid cell to obtain the cell 
! areas on the actual grid.
!
!=======================================================================
!
  allocate(dASG(0:nxSG-1, 0:nySG-1))
  call check( nf90_inq_varid(ncid, "area", dASG_varid) )
  call check( nf90_get_var(ncid, dASG_varid, dASG) )

  nx = nxSG/2
  ny = nySG/2

  WRITE(*,500)
  WRITE(*,500) "Dimensions of (actual, i.e. computational) grid are:"
  WRITE(*,520) "  nx = ", nx
  WRITE(*,520) "  ny = ", ny

  allocate(quarter_dA_ll(0:nx-1, 0:ny-1))
  allocate(quarter_dA_lr(0:nx-1, 0:ny-1))
  allocate(quarter_dA_ul(0:nx-1, 0:ny-1))
  allocate(quarter_dA_ur(0:nx-1, 0:ny-1))

  quarter_dA_ll = dASG(0:nxSG-1:2, 0:nySG-1:2)
  quarter_dA_lr = dASG(0:nxSG-1:2, 1:nySG-1:2)
  quarter_dA_ur = dASG(1:nxSG-1:2, 1:nySG-1:2)
  quarter_dA_ul = dASG(1:nxSG-1:2, 0:nySG-1:2)

  allocate(dA(0:nx-1, 0:ny-1))
  allocate(sqrt_dA(0:nx-1, 0:ny-1))

  dA = quarter_dA_ll + quarter_dA_lr + quarter_dA_ur + quarter_dA_ul
!
!=======================================================================
!
! Calculate a typical/representative cell size for each cell by taking
! the square root of the area of the cell.  Then calculate the minimum,
! maximum, and average cell sizes over the whole grid.
!
!=======================================================================
!
  sqrt_dA = sqrt(dA)
  min_cell_size = minval(sqrt_dA)
  max_cell_size = maxval(sqrt_dA)
  avg_cell_size = sum(sqrt_dA)/(nx*ny)

  WRITE(*,500)
  WRITE(*,500) "Minimum, maximum, and average cell sizes are (based on square"
  WRITE(*,500) "root of cell area):"
  WRITE(*,530) "  min_cell_size = ", min_cell_size
  WRITE(*,530) "  max_cell_size = ", max_cell_size
  WRITE(*,530) "  avg_cell_size = ", avg_cell_size
530 FORMAT(A, G)
!
!=======================================================================
!
! Use the average cell size to calculate an equivalent global uniform 
! cubed-sphere resolution (in units of number of cells) for the regional
! grid.  This is the RES that a global uniform (i.e. stretch factor of
! 1) cubed-sphere grid would need to have in order to have the same no-
! minal cell size as the average cell size of the regional grid.
!
!=======================================================================
!
  RES_equiv = nint( (2.0*pi_geom*radius_Earth)/(4.0*avg_cell_size) )

  WRITE(*,500)
  WRITE(*,500) "Equivalent global uniform cubed-sphere resolution is:"
  WRITE(*,530) "  RES_equiv = ", RES_equiv
!
!=======================================================================
!
! Write the average cell size and equivalent global resolution to the 
! grid file as a global attributes.
!
!=======================================================================
!
  WRITE(*,500)
  WRITE(*,500) "Writing avg_cell_size and RES_equiv to the grid specification"
  WRITE(*,500) "file as global attributes..."

  call check( nf90_redef(ncid) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "avg_cell_size", avg_cell_size) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "RES_equiv", RES_equiv) )
  call check( nf90_enddef(ncid) )

  call check( nf90_close(ncid) )

  WRITE(*,500)
  WRITE(*,500) "Done."

end program global_equiv_resol


subroutine check(status)
  use netcdf
  integer,intent(in) :: status
!
  if(status /= nf90_noerr) then
    write(0,*) ' check netcdf status = ', status
    write(0,'("error ", a)') trim(nf90_strerror(status))
    stop "Stopped"
  endif
end subroutine check
