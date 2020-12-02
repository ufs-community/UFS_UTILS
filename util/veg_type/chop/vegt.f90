 program chop

!---------------------------------------------------------
! Read in nesdis file.  Chop out a subset over a
! region such as N America.  Write the result to
! a netcdf file readable by the UFS_UTILS 
! sfc_climo_gen program.
!---------------------------------------------------------

 use netcdf

 implicit none

 character(len=150) :: filenetcdf, raw_file
 character(len=80), allocatable  :: att_name(:), gatt_name(:)
 character(len=150), allocatable :: att_value(:), gatt_value(:)

 integer, parameter :: isrc=43200
 integer, parameter :: jsrc=21600

 integer :: ncid, ncid_raw, status, dim_lat, dim_lon, dim_lat_p1, dim_lon_p1
 integer :: i, j, jj, lon, lat, dim_time
 integer :: id_lat, id_lon, landice
 integer :: id_lat_corner, id_lon_corner
 integer :: num_times, id_times, id_data, id_var
 integer :: times, natts, ngatts
 integer :: istart, iend, jstart, jend
 integer, parameter :: water_flag = 17
 integer(kind=1), allocatable :: global_dat(:,:), the_data(:,:)
 integer(kind=1), parameter    :: missing = -9

 real(kind=8)    :: lat11, lon11, dx, dy
 real(kind=8), allocatable :: lons(:), lons_corner(:)
 real(kind=8), allocatable :: lats(:), lats_corner(:)

!----------------------------------------------------------------------------
! Read in global raw file.
!----------------------------------------------------------------------------

 raw_file="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/fv3.vegt.new.tundra.netcdf/orig_data/NESDIS_VIIRS_SfcType.nc"

 print*,'Open: ',trim(raw_file)
 status=nf90_open(trim(raw_file), nf90_nowrite, ncid_raw)
 if (status /= nf90_noerr) stop 2

 allocate(global_dat(isrc,jsrc))
 status=nf90_inq_varid(ncid_raw, 'surface_type', id_var)
 if (status /= nf90_noerr) stop 4

 status=nf90_get_var(ncid_raw, id_var, global_dat)
 if (status /= nf90_noerr) stop 6

 print*,'Max/min values of raw data ',maxval(global_dat),minval(global_dat)

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

!---------------------------------------------------------------
! Subset the data to a regional grid. Note: logic may not
! work crossing the dateline.
!
! Set specs of the regional grid.
!---------------------------------------------------------------

! Conus
 istart = 5251
 iend = 14500
 jstart = 11501
 jend   = 18000

 lon = iend - istart + 1
 lat = jend - jstart + 1

 allocate(the_data(istart:iend,jstart:jend))
 do j = jstart, jend
   jj = j
 do i = istart, iend
   if (global_dat(i,jj) == water_flag) then
     the_data(i,j) = missing
   else
     the_data(i,j) = global_dat(i,jj)
   endif
 enddo
 enddo

 deallocate(global_dat)

 num_times = 1
 times = 0
 
 dx = 1.0_8/120.0_8
 dy = 1.0_8/120.0_8

 lat11 = -90.0_8 + dy*0.5_8
 lon11 = -180.0_8 + dx*0.5_8

 allocate(lons(istart:iend))
 allocate(lats(jstart:jend))

 do i = istart, iend
   lons(i) = real((i-1),8) * dx + lon11
   print*,'Lon of output grid ',i,lons(i)
 enddo

 do i = jstart, jend
   lats(i) = real((i-1),8) * dy + lat11
   print*,'Lat of output grid ',i,lats(i)
 enddo

! lat/lon of corner of grid box.

 allocate(lons_corner(istart:iend+1))  ! required one extra corner
                                       ! for non-periodic grids.
 allocate(lats_corner(jstart:jend+1))  ! requires one extra corner

 do i = istart, iend+1
   lons_corner(i) = (real((i-1),8) * dx + lon11) - 0.5_8*dx
   print*,'Corner lon of output grid ',i,lons_corner(i)

 enddo

 do i = jstart, jend+1
   lats_corner(i) = (real((i-1),8) * dy + lat11) - 0.5_8*dy
   print*,'Corner lat of output grid ',i,lats_corner(i)
 enddo

 filenetcdf="./vegt.nc"
 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, ior(nf90_netcdf4, nf90_classic_model), ncid)
 if (status /= nf90_noerr) stop 20

 status=nf90_def_dim(ncid, 'idim', lon, dim_lon)
 if (status /= nf90_noerr) stop 22

 status=nf90_def_dim(ncid, 'jdim', lat, dim_lat)
 if (status /= nf90_noerr) stop 24

 status=nf90_def_dim(ncid, 'idim_p1', (lon+1), dim_lon_p1)
 if (status /= nf90_noerr) stop 26

 status=nf90_def_dim(ncid, 'jdim_p1', (lat+1), dim_lat_p1)
 if (status /= nf90_noerr) stop 30

 status=nf90_def_dim(ncid, 'time', num_times, dim_time)
 if (status /= nf90_noerr) stop 32

 status=nf90_put_att(ncid, nf90_global, 'source', 'VIIRS-BASED IGBP VEGETATION TYPE')
 if (status /= nf90_noerr) stop 34

 status=nf90_put_att(ncid, nf90_global, 'projection', 'regular lat/lon')
 if (status /= nf90_noerr) stop 36

 do i = 1, ngatts
   status=nf90_put_att(ncid, nf90_global, gatt_name(i), gatt_value(i))
   if (status /= nf90_noerr) stop 40
 enddo

 status=nf90_def_var(ncid, 'time', nf90_float, dim_time, id_times)
 if (status /= nf90_noerr) stop 42

 status=nf90_def_var(ncid, 'lat', nf90_double, dim_lat, id_lat)
 if (status /= nf90_noerr) stop 44

 status=nf90_put_att(ncid, id_lat, 'long_name', 'grid cell center latitude')
 if (status /= nf90_noerr) stop 46

 status=nf90_def_var(ncid, 'lon', nf90_double, dim_lon, id_lon)
 if (status /= nf90_noerr) stop 50

 status=nf90_put_att(ncid, id_lon, 'long_name', 'grid cell center longitude')
 if (status /= nf90_noerr) stop 60

 status=nf90_def_var(ncid, 'lat_corner', nf90_double, dim_lat_p1, id_lat_corner)
 if (status /= nf90_noerr) stop 62

 status=nf90_put_att(ncid, id_lat_corner, 'long_name', 'grid cell corner latitude')
 if (status /= nf90_noerr) stop 64

 status=nf90_def_var(ncid, 'lon_corner', nf90_double, dim_lon_p1, id_lon_corner)
 if (status /= nf90_noerr) stop 66

 status=nf90_put_att(ncid, id_lon_corner, 'long_name', 'grid cell corner longitude')
 if (status /= nf90_noerr) stop 70

 status=nf90_def_var(ncid, 'vegetation_type', nf90_byte, (/dim_lon,dim_lat,dim_time/), id_data)
 if (status /= nf90_noerr) stop 76
 landice=15
 status=nf90_put_att(ncid, id_data, 'landice_category', landice)
 if (status /= nf90_noerr) stop 80

 status=nf90_put_att(ncid, id_data, 'missing_value', missing)
 if (status /= nf90_noerr) stop 90

 do i = 1, natts
   if (index(att_name(i), "class") /= 0) then
     status=nf90_put_att(ncid, id_data, att_name(i), att_value(i))
     if (status /= nf90_noerr) stop 92
  endif
 enddo

 status=nf90_put_att(ncid, id_times, 'units', 'days since 2015-1-1')
 if (status /= nf90_noerr) stop 94

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 96

 status=nf90_put_var(ncid, id_times, times)
 if (status /= nf90_noerr) stop 98

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 100

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 104

 status=nf90_put_var(ncid, id_lon_corner, lons_corner)
 if (status /= nf90_noerr) stop 109

 status=nf90_put_var(ncid, id_lat_corner, lats_corner)
 if (status /= nf90_noerr) stop 114

 status=nf90_put_var(ncid, id_data, the_data, start=(/1,1,1/), count=(/lon,lat,1/))
 if (status /= nf90_noerr) stop 117

 deallocate(the_data)

 status=nf90_close(ncid)

 end program chop
