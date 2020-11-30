 program vegt

! Read in nesdis data.  Interpolate it to a coarser
! resolution global grid.  Write result to a netcdf
! file readable by the UFS_UTILS sfc_climo_gen
! program.

 use netcdf

 implicit none

 character(len=150) :: filenetcdf, raw_file

 integer :: i, j, iraw, jraw
 integer :: lon, lat, nearest_i, nearest_j
 integer :: ncid_raw, ncid, status, ierr, iunit
 integer :: times, num_times, landice
 integer :: dim_lon, dim_lat, dim_lat_p1, dim_time
 integer :: id_lat_corner, id_lon_corner, id_var
 integer :: id_times, id_lat, id_lon, id_data
 integer(kind=8) :: offset
 integer, allocatable :: count(:,:,:)
 integer(kind=1), allocatable :: raw_data(:,:)
 integer :: num_vegt = 20
 integer :: vegtype

 real    :: missing = -999.9
 real(kind=8)    :: lat11, lon11, dx, dy, dx_raw, dy_raw
 real(kind=8)    :: lat11_raw, lon11_raw, x, y, lon_raw, lat_raw
 real(kind=8), allocatable :: lons(:), lons_corner(:)
 real(kind=8), allocatable :: lats(:), lats_corner(:)
 real, allocatable :: the_data(:,:)

 raw_file="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/fv3.vegt.new.tundra.netcdf/orig_data/NESDIS_VIIRS_SfcType.nc"

 print*,'open ',trim(raw_file)
 status=nf90_open(trim(raw_file), nf90_nowrite, ncid_raw)
 if (status /= nf90_noerr) stop 2

 iraw = 43200
 jraw = 21600

 dx_raw = 1.0_8/120.0_8
 dy_raw = 1.0_8/120.0_8
 
 lat11_raw =  -90.0_8 + dy_raw*0.5_8
 lon11_raw = -180.0_8 + dx_raw*0.5_8

 print*,'raw corner point ',lat11_raw,lon11_raw

 status=nf90_inq_varid(ncid_raw, 'surface_type', id_var)
 if (status /= nf90_noerr) stop 2

 allocate(raw_data(iraw,jraw))

 status=nf90_get_var(ncid_raw, id_var, raw_data)
 if (status /= nf90_noerr) stop 3

 print*,maxval(raw_data),minval(raw_data)

 status=nf90_close(ncid_raw)

 lon = 7200
 lat = 3600

 dx = 0.05_8
 dy = 0.05_8

 print*,lon, lat

 allocate(the_data(lon,lat))

 num_times = 1
 times = 0

 lat11 = -90.0_8 + dy*0.5_8
 lon11 = 0.0_8 + dx*0.5_8

 allocate(lons(lon))
 allocate(lats(lat))

 do i = 1, lon
   lons(i) = real((i-1),8) * dx + lon11
   print*,'lon ',i,lons(i)
 enddo

 do i = 1, lat
   lats(i) = real((i-1),8) * dy + lat11
   print*,'lat ',i,lats(i)
 enddo

! lat/lon of corner of grid box.

 allocate(lons_corner(lon))
 allocate(lats_corner(lat+1))  ! requires one extra corner

 do i = 1, lon
   lons_corner(i) = (real((i-1),8) * dx + lon11) - 0.5_8*dx
   print*,'corner lon ',i,lons_corner(i)
 enddo

 do i = 1, lat+1
   lats_corner(i) = (real((i-1),8) * dy + lat11) - 0.5_8*dy
   print*,'corner lat ',i,lats_corner(i)
 enddo

!---------------------------------------------------------------
! Interpolate data
!---------------------------------------------------------------

 allocate(count(lon,lat,num_vegt))
 count = 0

 do j = 1, jraw
   lat_raw = real((j-1),8) * dy_raw + lat11_raw
   y =  (lat_raw - lat11) / dy + 1.0
   nearest_j = nint(y)
   print*,'j is ',j,nearest_j
 do i = 1, iraw

   lon_raw = real((i-1),8) * dx_raw + lon11_raw
   x =  mod((lon_raw - lon11)+3600,360.0) / dx + 1.0
   nearest_i = nint(x)
   if (nearest_i > lon) nearest_i = nearest_i - lon

   vegtype = raw_data(i,j)

   count(nearest_i,nearest_j,vegtype) =  count(nearest_i,nearest_j,vegtype) + 1

 enddo
 enddo

! check for points with zero counts

 do j = 1, lat
 do i = 1, lon
   if (maxval(count(i,j,:)) == 0) then
     print*,'zero count at point ', i, j
     stop 33
   endif
 enddo
 enddo

 the_data = float(maxloc(count, dim=3))
 where(nint(the_data) == 17) the_data= missing

 filenetcdf = "./test.nc"
 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, ior(nf90_netcdf4, nf90_classic_model), ncid)
 if (status /= nf90_noerr) stop 1

 status=nf90_def_dim(ncid, 'idim', lon, dim_lon)
 if (status /= nf90_noerr) stop 3

 status=nf90_def_dim(ncid, 'jdim', lat, dim_lat)
 if (status /= nf90_noerr) stop 2

 status=nf90_def_dim(ncid, 'jdim_p1', (lat+1), dim_lat_p1)
 if (status /= nf90_noerr) stop 2

 status=nf90_def_dim(ncid, 'time', num_times, dim_time)
 if (status /= nf90_noerr) stop 5

 status=nf90_put_att(ncid, nf90_global, 'source', 'IGBP VEG TYPE')
 if (status /= nf90_noerr) stop 4

 status=nf90_put_att(ncid, nf90_global, 'projection', 'regular lat/lon')
 if (status /= nf90_noerr) stop 34

 status=nf90_def_var(ncid, 'time', nf90_float, dim_time, id_times)
 if (status /= nf90_noerr) stop 6

 status=nf90_def_var(ncid, 'lat', nf90_double, dim_lat, id_lat)
 if (status /= nf90_noerr) stop 17

 status=nf90_put_att(ncid, id_lat, 'long_name', 'grid cell center latitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'lon', nf90_double, dim_lon, id_lon)
 if (status /= nf90_noerr) stop 16

 status=nf90_put_att(ncid, id_lon, 'long_name', 'grid cell center longitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'lat_corner', nf90_double, dim_lat_p1, id_lat_corner)
 if (status /= nf90_noerr) stop 17

 status=nf90_put_att(ncid, id_lat_corner, 'long_name', 'grid cell corner latitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'lon_corner', nf90_double, dim_lon, id_lon_corner)
 if (status /= nf90_noerr) stop 16

 status=nf90_put_att(ncid, id_lon_corner, 'long_name', 'grid cell corner longitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'vegetation_type', nf90_float, (/dim_lon,dim_lat,dim_time/), id_data)
 if (status /= nf90_noerr) stop 6
 landice=15
 status=nf90_put_att(ncid, id_data, 'landice_category', landice)
 if (status /= nf90_noerr) stop 10

 status=nf90_put_att(ncid, id_data, 'missing_value', missing)
 if (status /= nf90_noerr) stop 10

 status=nf90_put_att(ncid, id_times, 'units', 'days since 2015-1-1')
 if (status /= nf90_noerr) stop 7

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 8


 status=nf90_put_var(ncid, id_times, times)
 if (status /= nf90_noerr) stop 9

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 19

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 20

 status=nf90_put_var(ncid, id_lon_corner, lons_corner)
 if (status /= nf90_noerr) stop 19

 status=nf90_put_var(ncid, id_lat_corner, lats_corner)
 if (status /= nf90_noerr) stop 20
 
 status=nf90_put_var(ncid, id_data, the_data, start=(/1,1,1/), count=(/lon,lat,1/))
 if (status /= nf90_noerr) stop 12

 status=nf90_close(ncid)

 stop

end program vegt
