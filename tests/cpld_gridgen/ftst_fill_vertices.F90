!> Unit test for reshape_staggers routine
!!
!! This test checks the filling of vertex values for each stagger
!! location
!! @author Denise.Worthen@noaa.gov
program ftst_fill_vertices

  use assertion_mod, only: assert_equal
  use gengrid_kinds, only: dbl_kind, int_kind, CL
  use gengrid_utils, only: fill_vertices
  use grdvars      , only: nv, iVertCt, jVertCt

  implicit none

  integer, parameter :: nx = 5, ny = 4
  integer, parameter :: maxtests = 50, nresults = maxtests

  logical           :: ispassing(nresults)
  character(len=CL) :: testmsg(nresults) = ' '

  integer           :: iind(2), jind(2)
  ! test data
  real(dbl_kind)    :: lon(nx, ny), lat(nx, ny)
  integer, dimension(nv) :: iVertCx, jVertCx
  real(dbl_kind)    :: lonvert(nx, ny, nv), latvert(nx, ny, nv)
  real(dbl_kind)    :: xlonCx(nx), xlatCx(nx)
  ! result data
  !real(dbl_kind), allocatable    :: cnlons(:), cnlats(:)
  !integer(int_kind), allocatable :: cnmask(:)
  !real(dbl_kind), allocatable    :: crlons(:, :), crlats(:, :)

  integer :: ng, nt, ntests
  integer :: i, j, n
  integer :: ib, ie, jb, je, idim, jdim
  integer :: idx, jdx, idx1

  character(len=CL) :: msg, msg_out
  logical :: status

  ! Initialize global test data; coordinate encoded
  do j = 1, ny
     do i = 1, nx
        lon(i,j) =  10.0_dbl_kind * i + j
        lat(i,j) = -10.0_dbl_kind * i - j
     end do
  end do

  do j = 1,ny
     print '(i3,5f8.2)',j,(lon(i,j),i=1,nx)
  end do
  print *
  do j = 1,ny
     print '(i3,5f8.2)',j,(lat(i,j),i=1,nx)
  end do
  print *

   ntests = nt
  !do nt = 1,ntests
  !   print '(i5,a)',nt,' '//trim(testmsg(nt))
  !end do

  !if (all(ispassing(1:ntests))) then
  !   print '(a)', 'All unit tests passed '
  !else
  !   print '(a)', 'FAIL: At least one unit test failed '
  !   stop 1
  !end if

end program ftst_fill_vertices
