!> @file
!! @brief Compute geo-referencing parameters for the Extended
!! Schmidt Gnomonic (ESG) regional grid.
!! @author R. J. Purser

!> Driver routine to compute geo-referencing parameters for
!! the Extended Schmidt Gnomonic (ESG) regional grid.
!! The parameters are:
!! - Geographic longitude (degrees)
!! - Geographic latitude (degrees)
!! - Grid edge 'x' distance (meters)
!! - Grid edge 'y' distance (meters)
!! - Area (meters squared)
!! - Grid vertex 'x' angle with respect to geographic east (degrees)
!! - Grid vertex 'y' angle with respect to geographic north (degrees)
!! @author R. J. Purser
!! @return 0 for success, error code otherwise.
program regional_grid

  use pkind, only: dp
  use pietc, only: dtor,rtod
  use pesg
  use netcdf

  implicit none

  ! namelist variables
  real(dp)                     :: plat,plon,pazi=0.0
  real(dp)                     :: delx,dely
  integer                      :: lx,ly
  namelist /regional_grid_nml/ plat,plon,pazi,delx,dely,lx,ly

  real(dp),parameter           :: re=6371000.0
  real(dp),parameter           :: lam=0.8

  integer                      :: nxh,nyh, nx,ny, nxm,nym
  logical                      :: ff

  real(dp),dimension(:,:),allocatable:: glat,glon
  real(dp),dimension(:,:),allocatable:: garea
  real(dp),dimension(:,:),allocatable:: dx,dy,angle_dx,angle_dy

  character(len=256)           :: nml_fname

  ! netcdf
  integer                      :: ncid
  integer                      :: string_dimid, nxp_dimid, nyp_dimid, nx_dimid, ny_dimid
  integer                      :: tile_varid, x_varid, y_varid, area_varid
  integer                      :: dx_varid, dy_varid, angle_dx_varid, angle_dy_varid
  integer, dimension(2)        :: dimids

  real(dp)                     :: a,k,m_arcx,m_arcy,q
  real(dp)                     :: m_delx,m_dely, delxre,delyre
  real(dp)                     :: arcx,arcy

!=============================================================================

  if (command_argument_count() == 1) then
    call get_command_argument(1, nml_fname)
  else
    nml_fname = "regional_grid.nml"
  end if

  open(10,file=trim(nml_fname),status="old",action="read")
  read(10,nml=regional_grid_nml)
  close(10)

  nxh=-lx
  nyh=-ly
  nx=2*nxh
  ny=2*nyh
  nxm=nx-1
  nym=ny-1

  allocate(glat(0:nx,0:ny))
  allocate(glon(0:nx,0:ny))
  allocate(garea(0:nxm,0:nym))

  allocate(dx(0:nxm,0:ny))
  allocate(dy(0:nx,0:nym))
  allocate(angle_dx(0:nx,0:ny))
  allocate(angle_dy(0:nx,0:ny))

  arcx=delx*nxh
  arcy=dely*nyh
  print'("arcx, arcy ",f8.4,f8.4)',arcx,arcy
  call bestesg_geo(lam,arcx,arcy, a,k,m_arcx,m_arcy,q,ff)
  if(ff)stop 'Failure flag returned from get_bestesg'
  print'("For lam=",f8.2," the best [smallest possible]")',lam
  print'("optimality criterion, Q, for this domain: ",e13.6 )',q
  print'("The corresponding optimal A and K:",2(1x,f8.4))',A,K
  print'("The corresponding m_arcx,y:",2(1x,e20.13))',m_arcx,m_arcy

  m_delx=m_arcx/nxh ! Map-space grid steps in x
  m_dely=m_arcy/nyh ! Map-space grid steps in y
  print'("x and y central grid resolution in map units:",2(1x,e12.5))',m_delx,m_dely

  print'("Get additional diagnostics from hgrid_ak_rr")'

  delxre=m_delx*re
  delyre=m_dely*re

  call hgrid_ak(lx,ly, nx,ny, A,K, plat*dtor,plon*dtor,pazi*dtor, re,delxre,delyre, &
                glat,glon,garea, dx,dy,angle_dx,angle_dy, ff)
  if(ff)stop 'Failure flag raised in hgrid routine'

  glon = glon*rtod
  glat = glat*rtod
  where (glon < 0.0) glon = glon + 360.0

  call check( nf90_create("regional_grid.nc", NF90_64BIT_OFFSET, ncid) )
  call check( nf90_def_dim(ncid, "string", 255, string_dimid) )
  call check( nf90_def_dim(ncid, "nx", nx, nx_dimid) )
  call check( nf90_def_dim(ncid, "ny", ny, ny_dimid) )
  call check( nf90_def_dim(ncid, "nxp", nx+1, nxp_dimid) )
  call check( nf90_def_dim(ncid, "nyp", ny+1, nyp_dimid) )

  call check( nf90_def_var(ncid, "tile", NF90_CHAR, [string_dimid], tile_varid) )
  call check( nf90_put_att(ncid, tile_varid, "standard_name", "grid_tile_spec") )

  dimids = (/ nxp_dimid, nyp_dimid /)
  call check( nf90_def_var(ncid, "x", NF90_DOUBLE, dimids, x_varid) )
  call check( nf90_put_att(ncid, x_varid, "standard_name", "geographic_longitude") )
  call check( nf90_put_att(ncid, x_varid, "units", "degree_east") )
  call check( nf90_put_att(ncid, x_varid, "hstagger", "C") )
  call check( nf90_def_var(ncid, "y", NF90_DOUBLE, dimids, y_varid) )
  call check( nf90_put_att(ncid, y_varid, "standard_name", "geographic_latitude") )
  call check( nf90_put_att(ncid, y_varid, "units", "degree_north") )
  call check( nf90_put_att(ncid, y_varid, "hstagger", "C") )
  dimids = (/ nx_dimid, ny_dimid /)
  call check( nf90_def_var(ncid, "area", NF90_DOUBLE, dimids, area_varid) )
  call check( nf90_put_att(ncid, area_varid, "standard_name", "grid_cell_area") )
  call check( nf90_put_att(ncid, area_varid, "units", "m2") )
  call check( nf90_put_att(ncid, area_varid, "hstagger", "H") )

  dimids = (/ nx_dimid, nyp_dimid /)
  call check( nf90_def_var(ncid, "dx", NF90_DOUBLE, dimids, dx_varid) )
  call check( nf90_put_att(ncid, dx_varid, "standard_name", "dx") )
  call check( nf90_put_att(ncid, dx_varid, "units", "m") )
  call check( nf90_put_att(ncid, dx_varid, "hstagger", "H") )

  dimids = (/ nxp_dimid, ny_dimid /)
  call check( nf90_def_var(ncid, "dy", NF90_DOUBLE, dimids, dy_varid) )
  call check( nf90_put_att(ncid, dy_varid, "standard_name", "dy") )
  call check( nf90_put_att(ncid, dy_varid, "units", "m") )
  call check( nf90_put_att(ncid, dy_varid, "hstagger", "H") )

  dimids = (/ nxp_dimid, nyp_dimid /)
  call check( nf90_def_var(ncid, "angle_dx", NF90_DOUBLE, dimids, angle_dx_varid) )
  call check( nf90_put_att(ncid, angle_dx_varid, "standard_name", "angle_dx") )
  call check( nf90_put_att(ncid, angle_dx_varid, "units", "deg") )
  call check( nf90_put_att(ncid, angle_dx_varid, "hstagger", "C") )
  call check( nf90_def_var(ncid, "angle_dy", NF90_DOUBLE, dimids, angle_dy_varid) )
  call check( nf90_put_att(ncid, angle_dy_varid, "standard_name", "angle_dy") )
  call check( nf90_put_att(ncid, angle_dy_varid, "units", "deg") )
  call check( nf90_put_att(ncid, angle_dy_varid, "hstagger", "C") )

  call check( nf90_put_att(ncid, NF90_GLOBAL, "history", "gnomonic_ed") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "source", "FV3GFS") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "grid", "akappa") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "plat", plat) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "plon", plon) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "pazi", pazi) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "delx", m_delx*rtod) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "dely", m_dely*rtod) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "lx", lx) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "ly", ly) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "a", a) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "k", k) )

  call check( nf90_enddef(ncid) )

  call check( nf90_put_var(ncid, tile_varid, "tile7") )
  call check( nf90_put_var(ncid, x_varid, glon) )
  call check( nf90_put_var(ncid, y_varid, glat) )
  call check( nf90_put_var(ncid, area_varid, garea) )
  call check( nf90_put_var(ncid, dx_varid, dx) )
  call check( nf90_put_var(ncid, dy_varid, dy) )
  call check( nf90_put_var(ncid, angle_dx_varid, angle_dx) )
  call check( nf90_put_var(ncid, angle_dy_varid, angle_dy) )

  call check( nf90_close(ncid) )

end program regional_grid

!> Check results of netCDF call.
!!
!! @param[in] status return code to check
!! @author R. J. Purser
subroutine check(status)
use netcdf
integer,intent(in) :: status
!
if(status /= nf90_noerr) then
  write(0,*)' check netcdf status=',status
  write(0,'("error ", a)')trim(nf90_strerror(status))
  stop "Stopped"
endif
end subroutine check
