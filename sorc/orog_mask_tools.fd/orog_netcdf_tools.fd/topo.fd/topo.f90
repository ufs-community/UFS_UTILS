 program topo_netcdf

! Convert the gmted20120 data to netcdf.

 use netcdf

 implicit none

 integer*4, parameter :: idim=43200
 integer*4, parameter :: jdim=21600
 integer*4, parameter :: idim_p1=43201
 integer*4, parameter :: jdim_p1=21601

 character(len=150) :: filenetcdf, fileraw

 integer :: i, istat, ncid, status, dim_i, dim_j
 integer :: dim_ip1, dim_jp1
 integer :: id_lon, id_lat, id_data
 integer :: id_lat_corner, id_lon_corner

 integer(kind=2), allocatable :: topo(:,:)

 real(kind=8), allocatable :: lats(:), lons(:)
 real(kind=8), allocatable :: lats_corner(:), lons_corner(:)
 real(kind=8)    :: lat11, lon11, dx, dy

 dx = 1.0_8/120.0_8
 dy = -(1.0_8/120.0_8)

 lat11 =  90.0_8 + dy*0.5_8
 lon11 = -180.0_8 + dx*0.5_8

 allocate(lons(idim),lats(jdim),topo(idim,jdim))
 allocate(lons_corner(idim_p1),lats_corner(jdim_p1))

 do i = 1, idim
   lons(i) = real((i-1),8) * dx + lon11
   print*,'lon ',i,lons(i)
 enddo

 do i = 1, jdim
   lats(i) = real((i-1),8) * dy + lat11
   print*,'lat ',i,lats(i)
 enddo

 lat11 =  90.0_8
 lon11 = -180.0_8

 do i = 1, idim_p1
   lons_corner(i) = real((i-1),8) * dx + lon11
   print*,'lon_corner ',i,lons_corner(i)
 enddo

 do i = 1, jdim_p1
   lats_corner(i) = real((i-1),8) * dy + lat11
   print*,'lat_corner ',i,lats_corner(i)
 enddo

 fileraw="/scratch1/NCEPDEV/global/glopara/fix/raw/orog/gmted2010.30sec.int"

 open(11, file=trim(fileraw), access='direct', recl=idim*jdim*2)
 read(11, rec=1, iostat=istat) topo
 if (istat /= 0) stop 99
 close(11)

 print*,'topo ', maxval(topo),minval(topo)

 filenetcdf="./gmted2010.30sec.nc"

 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid)
 if (status /= nf90_noerr) stop 1

 status=nf90_def_dim(ncid, 'idim', idim, dim_i)
 if (status /= nf90_noerr) stop 3

 status=nf90_def_dim(ncid, 'jdim', jdim, dim_j)
 if (status /= nf90_noerr) stop 2

 status=nf90_def_dim(ncid, 'idim_p1', (idim+1), dim_ip1)
 if (status /= nf90_noerr) stop 4

 status=nf90_def_dim(ncid, 'jdim_p1', (jdim+1), dim_jp1)
 if (status /= nf90_noerr) stop 5

 status=nf90_put_att(ncid, nf90_global, 'source', 'USGS GMTED2010 TOPOGRAPHY DATA')
 if (status /= nf90_noerr) stop 6

 status=nf90_put_att(ncid, nf90_global, 'projection', 'regular lat/lon')
 if (status /= nf90_noerr) stop 67

 status=nf90_def_var(ncid, 'lat', nf90_double, dim_j, id_lat)
 if (status /= nf90_noerr) stop 17

 status=nf90_put_att(ncid, id_lat, 'long_name', 'grid cell center latitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_put_att(ncid, id_lat, 'units', 'degrees')
 if (status /= nf90_noerr) stop 85

 status=nf90_def_var(ncid, 'lat_corner', nf90_double, dim_jp1, id_lat_corner)
 if (status /= nf90_noerr) stop 37

 status=nf90_put_att(ncid, id_lat_corner, 'long_name', 'grid cell corner latitude')
 if (status /= nf90_noerr) stop 38

 status=nf90_put_att(ncid, id_lat_corner, 'units', 'degrees')
 if (status /= nf90_noerr) stop 86

 status=nf90_def_var(ncid, 'lon', nf90_double, dim_i, id_lon)
 if (status /= nf90_noerr) stop 16

 status=nf90_put_att(ncid, id_lon, 'long_name', 'grid cell center longitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_put_att(ncid, id_lon, 'units', 'degrees')
 if (status /= nf90_noerr) stop 87

 status=nf90_def_var(ncid, 'lon_corner', nf90_double, dim_ip1, id_lon_corner)
 if (status /= nf90_noerr) stop 16

 status=nf90_put_att(ncid, id_lon_corner, 'long_name', 'grid cell corner longitude')
 if (status /= nf90_noerr) stop 40

 status=nf90_put_att(ncid, id_lon_corner, 'units', 'degrees')
 if (status /= nf90_noerr) stop 88

 status=nf90_def_var(ncid, 'topo', nf90_short, (/dim_i,dim_j/), id_data)
 if (status /= nf90_noerr) stop 20

 status=nf90_put_att(ncid, id_data, 'long_name', 'topography')
 if (status /= nf90_noerr) stop 65

 status=nf90_put_att(ncid, id_data, 'units', 'meters')
 if (status /= nf90_noerr) stop 55

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 22

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 19

 status=nf90_put_var(ncid, id_lon_corner, lons_corner)
 if (status /= nf90_noerr) stop 59

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 20

 status=nf90_put_var(ncid, id_lat_corner, lats_corner)
 if (status /= nf90_noerr) stop 57

 status=nf90_put_var(ncid, id_data, topo)
 if (status /= nf90_noerr) stop 24

 status=nf90_close(ncid)

 end program topo_netcdf
