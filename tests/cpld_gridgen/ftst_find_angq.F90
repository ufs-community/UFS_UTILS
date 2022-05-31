! Unit test for cpld_gridgen routine "find_angq".
!
! Reads a sample MOM6 supergrid and calculates the
! rotation angle on corner points
!
! Author Denise Worthen 2/08/2022

  program ftst_find_angq
 
  use netcdf
  use grdvars,       only : ni,nj,nx,ny
  use grdvars,       only : x,y,xsgp1,ysgp1,sg_maxlat
  use grdvars,       only : angq
  use angles,        only : find_angq

  implicit none

  integer :: i,j,i1, i2
  integer :: rc, ncid, id

  logical :: mastertask = .false.
  logical :: debug = .false.

  ! pole locations on SG
  integer :: ipolesg(2)
  ! unit test values
  real(kind=8) :: puny = 1.0e-12
  real(kind=8) :: delta(15)
  real(kind=8) :: sumdelta

  print *,"Starting test of cpld_gridgen routine find_angq"

  ! 1deg MOM6 dimensions
  ni = 360
  nj = 320
  ! super grid dimensions
  nx = 2*ni
  ny = 2*nj

  ! supergrid x,y and angles on corners
  allocate (x(0:nx,0:ny), y(0:nx,0:ny), angq(0:nx,0:ny))
  ! supergrid "plus 1" arrays
  allocate (xsgp1(0:nx,0:ny+1), ysgp1(0:nx,0:ny+1))

  !open the supergrid file and read the x,y coords
  rc = nf90_open('./data/ocean_hgrid.nc', nf90_nowrite, ncid)
  rc = nf90_inq_varid(ncid, 'x', id)  !lon
  rc = nf90_get_var(ncid,    id,  x)

  rc = nf90_inq_varid(ncid, 'y', id)  !lat
  rc = nf90_get_var(ncid,    id,  y)
  rc = nf90_close(ncid)

  ! max lat on supergrid
  sg_maxlat = maxval(y)

  !pole index on supergrid
  ipolesg = -1
      j = ny
  do i = 1,nx/2
   if(y(i,j) .eq. sg_maxlat)ipolesg(1) = i
  enddo
  do i = nx/2+1,nx
   if(y(i,j) .eq. sg_maxlat)ipolesg(2) = i
  enddo

  ! test angleq calculation
  call find_angq

  ! required for checking longitudes across seam
  where(xsgp1 .lt. 0.0)xsgp1 = xsgp1 + 360.0

  j = ny+1
  i1 = ipolesg(1); i2 = ipolesg(2)-(ipolesg(1)-i1)
  delta = 0.0
  ! check lons match across seam
  delta( 1) = xsgp1(i1-2,j)-xsgp1(i2+2,j)
  delta( 2) = xsgp1(i1-1,j)-xsgp1(i2+1,j)
  delta( 3) = xsgp1(i1,  j)-xsgp1(i2,  j)
  delta( 4) = xsgp1(i1+1,j)-xsgp1(i2-1,j)
  delta( 5) = xsgp1(i1+2,j)-xsgp1(i2-2,j)
  ! check lats match across seam
  delta( 6) = ysgp1(i1-2,j)-ysgp1(i2+2,j)
  delta( 7) = ysgp1(i1-1,j)-ysgp1(i2+1,j)
  delta( 8) = ysgp1(i1,  j)-ysgp1(i2,  j)
  delta( 9) = ysgp1(i1+1,j)-ysgp1(i2-1,j)
  delta(10) = ysgp1(i1+2,j)-ysgp1(i2-2,j)
  ! check angq match across seam
  j = ny
  delta(11)=angq(i1-2,j)-angq(i2-2,j)
  delta(12)=angq(i1-1,j)-angq(i2-1,j)
  delta(13)=angq(i1,  j)-angq(i2,  j)
  delta(14)=angq(i1+1,j)-angq(i2+1,j)
  delta(15)=angq(i1+2,j)-angq(i2+2,j)

  sumdelta = 0.0
  sumdelta = sum(delta)

  if (sumdelta >= puny) then
    print *,'OK'
    print *,'SUCCESS!'
    deallocate(x,y,xsgp1,ysgp1,angq)
  else
    print *,'ftst_find_angq failed'
    stop 1
  endif

  end program ftst_find_angq
