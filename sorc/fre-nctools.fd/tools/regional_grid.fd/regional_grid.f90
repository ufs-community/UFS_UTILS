!=============================================================================
program regional_grid
!=============================================================================

  use pkind, only: dp
  use pietc, only: dtor,rtod
  use netcdf

  implicit none

  ! namelist variables
  real(dp)                     :: plat,plon,pazi=0.0
  real(dp)                     :: delx,dely
  integer                      :: lx,ly
  real(dp)                     :: a,k
  namelist /regional_grid_nml/ plat,plon,delx,dely,lx,ly,a,k

  real(dp),parameter           :: re=6371000.0

  real(dp)                     :: redelx,redely
  integer                      :: nx,nxm, ny,nym
  logical                      :: ff

  real(dp),dimension(:,:),allocatable:: glat,glon
  real(dp),dimension(:,:),allocatable:: garea

  character(len=256)           :: nml_fname

  ! netcdf
  integer                      :: ncid
  integer                      :: string_dimid, nxp_dimid, nyp_dimid, nx_dimid, ny_dimid
  integer                      :: tile_varid, x_varid, y_varid, area_varid
  integer, dimension(2)        :: dimids

!=============================================================================

  if (command_argument_count() == 1) then
    call get_command_argument(1, nml_fname)
  else
    nml_fname = "regional_grid.nml"
  end if

  open(10,file=trim(nml_fname),status="old",action="read")
  read(10,nml=regional_grid_nml)
  close(10)

  nx=-lx*2
  nxm=nx-1
  ny=-ly*2
  nym=ny-1

  redelx=re*(delx*dtor)
  redely=re*(dely*dtor)

  allocate(glat(0:nx,0:ny))
  allocate(glon(0:nx,0:ny))
  allocate(garea(0:nxm,0:nym))

  call hgrid_ak(lx,ly,nx,ny,a,k,plat*dtor,plon*dtor,pazi*dtor, &
                re,redelx,redely, glat,glon,garea, ff)
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

  call check( nf90_put_att(ncid, NF90_GLOBAL, "history", "gnomonic_ed") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "source", "FV3GFS") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "grid", "akappa") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "plat", plat) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "plon", plon) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "pazi", pazi) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "delx", delx) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "dely", dely) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "lx", lx) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "ly", ly) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "a", a) )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "k", k) )

  call check( nf90_enddef(ncid) )

  call check( nf90_put_var(ncid, tile_varid, "tile7") )
  call check( nf90_put_var(ncid, x_varid, glon) )
  call check( nf90_put_var(ncid, y_varid, glat) )
  call check( nf90_put_var(ncid, area_varid, garea) )

  call check( nf90_close(ncid) )

end program regional_grid

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
