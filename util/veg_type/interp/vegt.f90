 program vegt

!----------------------------------------------------
! Read in nesdis data.  Interpolate it to a coarser
! resolution global grid.  Write result to a netcdf
! file readable by the UFS_UTILS sfc_climo_gen
! program.
!----------------------------------------------------

 use netcdf

 implicit none

 character(len=150) :: filenetcdf, raw_file
 character(len=80), allocatable  :: att_name(:), gatt_name(:)
 character(len=150), allocatable :: att_value(:), gatt_value(:)

 integer :: i, j, iraw, jraw
 integer :: lon, lat, nearest_i, nearest_j
 integer :: ncid_raw, ncid, status
 integer :: times, num_times, landice
 integer :: dim_lon, dim_lat, dim_lat_p1, dim_time
 integer :: id_lat_corner, id_lon_corner, id_var
 integer :: id_times, id_lat, id_lon, id_data
 integer, allocatable :: count(:,:,:)
 integer :: vegtype, natts, ngatts
 integer, parameter           :: num_vegt = 20
 integer(kind=1), parameter   :: missing = -9
 integer(kind=1), allocatable :: raw_data(:,:)
 integer(kind=1), allocatable :: the_data(:,:)

 real(kind=8)    :: lat11, lon11, dx, dy, dx_raw, dy_raw
 real(kind=8)    :: lat11_raw, lon11_raw, x, y, lon_raw, lat_raw
 real(kind=8), allocatable :: lons(:), lons_corner(:)
 real(kind=8), allocatable :: lats(:), lats_corner(:)

!------------------------------------------------------------------------
! The raw file is the original 30-sec data from NESDIS/Mike B.
! Open file and read in data.
!------------------------------------------------------------------------

 raw_file="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/fv3.vegt.new.tundra.netcdf/orig_data/NESDIS_VIIRS_SfcType.nc"

 print*,'Open raw file: ',trim(raw_file)
 status=nf90_open(trim(raw_file), nf90_nowrite, ncid_raw)
 if (status /= nf90_noerr) stop 2

 iraw = 43200
 jraw = 21600

 dx_raw = 1.0_8/120.0_8
 dy_raw = 1.0_8/120.0_8
 
 lat11_raw =  -90.0_8 + dy_raw*0.5_8
 lon11_raw = -180.0_8 + dx_raw*0.5_8

 print*,'Raw file corner point ',lat11_raw,lon11_raw

 status=nf90_inq_varid(ncid_raw, 'surface_type', id_var)
 if (status /= nf90_noerr) stop 4

 allocate(raw_data(iraw,jraw))
 status=nf90_get_var(ncid_raw, id_var, raw_data)
 if (status /= nf90_noerr) stop 6

 print*,'Raw file max/min values: ', maxval(raw_data),minval(raw_data)

!------------------------------------------------------------------------
! Read global and vegetation class attributes and save.  These
! will be written later to the output file.
!------------------------------------------------------------------------

 status=nf90_inquire_variable(ncid_raw, id_var, natts=natts)
 if (status /= nf90_noerr) stop 8

 print*,'Number of variable attributes ', natts

 allocate(att_name(natts))
 att_name = " "
 allocate(att_value(natts))
 att_value = " "
 do i = 1, natts
   status=nf90_inq_attname(ncid_raw, id_var, i, att_name(i))
   if (status /= nf90_noerr) stop 10
   if (index(att_name(i), "class") /= 0) then  ! only save veg class attributes
     status=nf90_get_att(ncid_raw, id_var, att_name(i), att_value(i))
     if (status /= nf90_noerr) stop 12
     print*,'Class attribute ',i, trim(att_name(i)), ' ', trim(att_value(i))
   endif
 enddo

 status=nf90_inquire_variable(ncid_raw, nf90_global, natts=ngatts)
 if (status /= nf90_noerr) stop 14

 print*,'Number of global attributes ', ngatts

 allocate(gatt_name(ngatts))
 gatt_name = ' '
 allocate(gatt_value(ngatts))
 gatt_value = ' '
 do i = 1, ngatts
   status=nf90_inq_attname(ncid_raw, nf90_global, i, gatt_name(i))
   if (status /= nf90_noerr) stop 16
   status=nf90_get_att(ncid_raw, nf90_global, gatt_name(i), gatt_value(i))
   if (status /= nf90_noerr) stop 18
   print*,'Global attribute ',i, trim(gatt_name(i)), ' ', trim(gatt_value(i))
 enddo

 status=nf90_close(ncid_raw)

!------------------------------------------------------------------------
! Set grid specs for coarser output file.
!------------------------------------------------------------------------

! 0.03-degree global grid
 lon = 12000
 lat = 6000
 dx = 0.03_8
 dy = 0.03_8

! 0.05-degree global grid
!lon = 7200
!lat = 3600
!dx = 0.05_8
!dy = 0.05_8

! 0.10-degree global grid
!lon = 3600
!lat = 1800
!dx = 0.1_8
!dy = 0.1_8

 allocate(the_data(lon,lat))

 num_times = 1
 times = 0

 lat11 = -90.0_8 + dy*0.5_8
 lon11 = 0.0_8 + dx*0.5_8

 allocate(lons(lon))
 allocate(lats(lat))

 do i = 1, lon
   lons(i) = real((i-1),8) * dx + lon11
   print*,'Lon of output file ',i,lons(i)
 enddo

 do i = 1, lat
   lats(i) = real((i-1),8) * dy + lat11
   print*,'Lat of output file ',i,lats(i)
 enddo

! lat/lon of corner of grid box.

 allocate(lons_corner(lon))
 allocate(lats_corner(lat+1))  ! requires one extra corner

 do i = 1, lon
   lons_corner(i) = (real((i-1),8) * dx + lon11) - 0.5_8*dx
   print*,'Corner lon of output file ',i,lons_corner(i)
 enddo

 do i = 1, lat+1
   lats_corner(i) = (real((i-1),8) * dy + lat11) - 0.5_8*dy
   print*,'Corner lat of output file ',i,lats_corner(i)
 enddo

!---------------------------------------------------------------
! Interpolate data.  Accumulate the input grid points within
! each output grid box.  Then pick the predominate category.
!---------------------------------------------------------------

 allocate(count(lon,lat,num_vegt))
 count = 0

 do j = 1, jraw
   lat_raw = real((j-1),8) * dy_raw + lat11_raw
   y =  (lat_raw - lat11) / dy + 1.0
   nearest_j = nint(y)
   print*,'Process j ',j,nearest_j
 do i = 1, iraw

   lon_raw = real((i-1),8) * dx_raw + lon11_raw
   x =  mod((lon_raw - lon11)+3600,360.0) / dx + 1.0
   nearest_i = nint(x)
   if (nearest_i > lon) nearest_i = nearest_i - lon

   vegtype = raw_data(i,j)

   count(nearest_i,nearest_j,vegtype) =  count(nearest_i,nearest_j,vegtype) + 1

 enddo
 enddo

 deallocate(raw_data)

! Check for points with zero counts

 do j = 1, lat
 do i = 1, lon
   if (maxval(count(i,j,:)) == 0) then
     print*,'Error. Zero count at point ', i, j
     stop 33
   endif
 enddo
 enddo

! Mask out points that are water - cat 17.

 the_data = maxloc(count, dim=3)
 where(the_data == 17) the_data= missing

 deallocate(count)

!---------------------------------------------------------------
! Open output file and write data.
!---------------------------------------------------------------

 filenetcdf = "./vegt.nc"
 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, ior(nf90_netcdf4, nf90_classic_model), ncid)
 if (status /= nf90_noerr) stop 20

 status=nf90_def_dim(ncid, 'idim', lon, dim_lon)
 if (status /= nf90_noerr) stop 22

 status=nf90_def_dim(ncid, 'jdim', lat, dim_lat)
 if (status /= nf90_noerr) stop 24

 status=nf90_def_dim(ncid, 'jdim_p1', (lat+1), dim_lat_p1)
 if (status /= nf90_noerr) stop 26

 status=nf90_def_dim(ncid, 'time', num_times, dim_time)
 if (status /= nf90_noerr) stop 30

 status=nf90_put_att(ncid, nf90_global, 'source', 'VIIRS-BASED IGBP VEGETATION TYPE')
 if (status /= nf90_noerr) stop 32

 status=nf90_put_att(ncid, nf90_global, 'projection', 'regular lat/lon')
 if (status /= nf90_noerr) stop 34

 do i = 1, ngatts
   status=nf90_put_att(ncid, nf90_global, gatt_name(i), gatt_value(i))
   if (status /= nf90_noerr) stop 35
 enddo

 status=nf90_def_var(ncid, 'time', nf90_float, dim_time, id_times)
 if (status /= nf90_noerr) stop 37

 status=nf90_def_var(ncid, 'lat', nf90_double, dim_lat, id_lat)
 if (status /= nf90_noerr) stop 40

 status=nf90_put_att(ncid, id_lat, 'long_name', 'grid cell center latitude')
 if (status /= nf90_noerr) stop 42

 status=nf90_def_var(ncid, 'lon', nf90_double, dim_lon, id_lon)
 if (status /= nf90_noerr) stop 44

 status=nf90_put_att(ncid, id_lon, 'long_name', 'grid cell center longitude')
 if (status /= nf90_noerr) stop 46

 status=nf90_def_var(ncid, 'lat_corner', nf90_double, dim_lat_p1, id_lat_corner)
 if (status /= nf90_noerr) stop 48

 status=nf90_put_att(ncid, id_lat_corner, 'long_name', 'grid cell corner latitude')
 if (status /= nf90_noerr) stop 52

 status=nf90_def_var(ncid, 'lon_corner', nf90_double, dim_lon, id_lon_corner)
 if (status /= nf90_noerr) stop 54

 status=nf90_put_att(ncid, id_lon_corner, 'long_name', 'grid cell corner longitude')
 if (status /= nf90_noerr) stop 56

 status=nf90_def_var(ncid, 'vegetation_type', nf90_byte, (/dim_lon,dim_lat,dim_time/), id_data)
 if (status /= nf90_noerr) stop 60
 landice=15
 status=nf90_put_att(ncid, id_data, 'landice_category', landice)
 if (status /= nf90_noerr) stop 62

 status=nf90_put_att(ncid, id_data, 'missing_value', missing)
 if (status /= nf90_noerr) stop 64

 do i = 1, natts
   if (index(att_name(i), "class") /= 0) then
     status=nf90_put_att(ncid, id_data, att_name(i), att_value(i))
     if (status /= nf90_noerr) stop 68
  endif
 enddo

 status=nf90_put_att(ncid, id_times, 'units', 'days since 2015-1-1')
 if (status /= nf90_noerr) stop 70

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 80

 status=nf90_put_var(ncid, id_times, times)
 if (status /= nf90_noerr) stop 82

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 84

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 85

 status=nf90_put_var(ncid, id_lon_corner, lons_corner)
 if (status /= nf90_noerr) stop 87

 status=nf90_put_var(ncid, id_lat_corner, lats_corner)
 if (status /= nf90_noerr) stop 90
 
 status=nf90_put_var(ncid, id_data, the_data, start=(/1,1,1/), count=(/lon,lat,1/))
 if (status /= nf90_noerr) stop 92

 status=nf90_close(ncid)

 print*,'NORMAL COMPLETION.'

 stop

end program vegt
