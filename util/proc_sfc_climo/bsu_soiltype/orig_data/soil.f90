 program bnu_soil

! Read each 1200x1200 tile file. Then piece the data together and
! output to a netcdf file readable by UFS_UTILS.

 use netcdf

 implicit none

 integer, parameter :: idim=43200
 integer, parameter :: jdim=21600
 integer, parameter :: itile=1200
 integer, parameter :: jtile=1200

 character(len=5) :: istart_ch, iend_ch, jstart_ch, jend_ch
 character(len=23) :: file_name
 character(len=100) :: file_path
 character(len=140) :: raw_file
 character(len=9) :: filenetcdf

 integer         :: i, j, ierr, ncid, status
 integer         :: dim_lon, dim_lat, dim_time
 integer         :: dim_lon_p1, dim_lat_p1
 integer         :: id_data, id_lat, id_lon
 integer         :: id_lat_corner, id_lon_corner, id_times
 integer         :: istart, iend, jstart, jend
 integer         :: num_times
 integer(kind=1) :: soil_tile(itile,jtile)
 integer(kind=1) :: soil(idim,jdim)
 integer(kind=1), parameter :: missing = -9
 integer, parameter :: water_flag = 14
 integer, parameter :: landice=16

 real(kind=8)    :: lat11, lon11, dx, dy
 real(kind=8), allocatable :: lons(:), lons_corner(:)
 real(kind=8), allocatable :: lats(:), lats_corner(:)

 file_path="/lfs/h2/emc/global/noscrub/George.Gayno/ufs_utils.git/bnu.soil/orig_data/bnu_soiltype_top/"

! Loop to read the tile files.

 do i=1, idim, itile
  istart=i
  iend=istart+itile-1
  do j=1, jdim, jtile

    jstart=j
    jend=jstart+jtile-1
    write(istart_ch,'(i5.5)') istart
    write(iend_ch,'(i5.5)') iend
    write(jstart_ch,'(i5.5)') jstart
    write(jend_ch,'(i5.5)') jend
    file_name= istart_ch // "-" // iend_ch // "." // jstart_ch // "-" // jend_ch
    raw_file= trim(file_path) // file_name

    print*,'open and read file ',trim(raw_file)
    open(14, file=raw_file, &
         access='direct', recl=itile*jtile, iostat=ierr, status='old', form='unformatted')
    if (ierr /= 0) stop 2

    read(14,rec=1,iostat=ierr) soil_tile
    if (ierr /= 0) stop 3

    close(14)

    print*,'max/min values ',maxval(soil_tile),minval(soil_tile)

    where(soil_tile == water_flag) soil_tile = missing

    soil(istart:iend,jstart:jend) = soil_tile

   enddo
 enddo

! Compute lat/lons of the grid.

 dx = 1.0_8/120.0_8
 dy = 1.0_8/120.0_8

 lat11 = -90.0_8 + dy*0.5_8
 lon11 = -180.0_8 + dx*0.5_8

 allocate(lons(idim))
 allocate(lats(jdim))

 do i = 1, idim
   lons(i) = real((i-1),8) * dx + lon11
   print*,'Lon of output grid ',i,lons(i)
 enddo

 do i = 1, jdim
   lats(i) = real((i-1),8) * dy + lat11
   print*,'Lat of output grid ',i,lats(i)
 enddo

 allocate(lons_corner(idim+1))  ! required one extra corner
                                       ! for non-periodic grids.
 allocate(lats_corner(jdim+1))  ! requires one extra corner

 do i = 1, idim+1
   lons_corner(i) = (real((i-1),8) * dx + lon11) - 0.5_8*dx
   print*,'Corner lon of output grid ',i,lons_corner(i)

 enddo

 do i = 1, jdim+1
   lats_corner(i) = (real((i-1),8) * dy + lat11) - 0.5_8*dy
   print*,'Corner lat of output grid ',i,lats_corner(i)
 enddo

 filenetcdf="./soil.nc"
 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, nf90_netcdf4, ncid)
 if (status /= nf90_noerr) stop 20

 status=nf90_def_dim(ncid, 'idim', idim, dim_lon)
 if (status /= nf90_noerr) stop 22

 status=nf90_def_dim(ncid, 'jdim', jdim, dim_lat)
 if (status /= nf90_noerr) stop 24
 
 status=nf90_def_dim(ncid, 'idim_p1', (idim+1), dim_lon_p1)
 if (status /= nf90_noerr) stop 26

 status=nf90_def_dim(ncid, 'jdim_p1', (jdim+1), dim_lat_p1)
 if (status /= nf90_noerr) stop 30

 num_times=1
 status=nf90_def_dim(ncid, 'time', num_times, dim_time)
 if (status /= nf90_noerr) stop 32

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

 status=nf90_def_var(ncid, 'time', nf90_float, dim_time, id_times)
 if (status /= nf90_noerr) stop 42

 status=nf90_put_att(ncid, id_times, 'units', 'days since 2015-1-1')
 if (status /= nf90_noerr) stop 94

 status=nf90_def_var(ncid, 'soil_type', nf90_byte, (/dim_lon,dim_lat,dim_time/), id_data)
 if (status /= nf90_noerr) stop 76

 status=nf90_put_att(ncid, id_data, 'landice_category', landice)
 if (status /= nf90_noerr) stop 80

 status=nf90_put_att(ncid, id_data, 'missing_value', missing)
 if (status /= nf90_noerr) stop 90

 status=nf90_put_att(ncid, id_data, 'class_01', 'Sand')
 if (status /= nf90_noerr) stop 91

 status=nf90_put_att(ncid, id_data, 'class_02', 'Loamy Sand')
 if (status /= nf90_noerr) stop 121

 status=nf90_put_att(ncid, id_data, 'class_03', 'Sandy Loam')
 if (status /= nf90_noerr) stop 122

 status=nf90_put_att(ncid, id_data, 'class_04', 'Silt Loam')
 if (status /= nf90_noerr) stop 123

 status=nf90_put_att(ncid, id_data, 'class_05', 'Silt')
 if (status /= nf90_noerr) stop 124

 status=nf90_put_att(ncid, id_data, 'class_06', 'Loam')
 if (status /= nf90_noerr) stop 125

 status=nf90_put_att(ncid, id_data, 'class_07', 'Sandy Clay Loam')
 if (status /= nf90_noerr) stop 126

 status=nf90_put_att(ncid, id_data, 'class_08', 'Silty Clay Loam')
 if (status /= nf90_noerr) stop 127

 status=nf90_put_att(ncid, id_data, 'class_09', 'Clay Loam')
 if (status /= nf90_noerr) stop 128

 status=nf90_put_att(ncid, id_data, 'class_10', 'Sandy Clay')
 if (status /= nf90_noerr) stop 129

 status=nf90_put_att(ncid, id_data, 'class_11', 'Silty Clay')
 if (status /= nf90_noerr) stop 130

 status=nf90_put_att(ncid, id_data, 'class_12', 'Clay')
 if (status /= nf90_noerr) stop 131

 status=nf90_put_att(ncid, id_data, 'class_13', 'Organic Material')
 if (status /= nf90_noerr) stop 132

 status=nf90_put_att(ncid, id_data, 'class_14', 'Water')
 if (status /= nf90_noerr) stop 133

 status=nf90_put_att(ncid, id_data, 'class_15', 'Bedrock')
 if (status /= nf90_noerr) stop 134

 status=nf90_put_att(ncid, id_data, 'class_16', 'Permanent Ice')
 if (status /= nf90_noerr) stop 135

 status=nf90_put_att(ncid, nf90_global, 'source', 'Beijing Normal University (BNU) SOIL TYPE')
 if (status /= nf90_noerr) stop 34

 status=nf90_put_att(ncid, nf90_global, 'projection', 'regular lat/lon')
 if (status /= nf90_noerr) stop 36

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 96

 status=nf90_put_var(ncid, id_times, 0)
 if (status /= nf90_noerr) stop 98

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 100

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 104

 status=nf90_put_var(ncid, id_lon_corner, lons_corner)
 if (status /= nf90_noerr) stop 109

 status=nf90_put_var(ncid, id_lat_corner, lats_corner)
 if (status /= nf90_noerr) stop 114

 status=nf90_put_var(ncid, id_data, soil, start=(/1,1,1/), count=(/idim,jdim,num_times/))
 if (status /= nf90_noerr) stop 117
 
 status=nf90_close(ncid)

 print*,'DONE'

 stop
 end program bnu_soil
