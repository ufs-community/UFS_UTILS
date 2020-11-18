 program soilveg

! Read in nesdis file.  Chop out a subset over a
! region such as N America.  Write the result to
! a netcdf file readable by the UFS_UTILS 
! sfc_climo_gen program.

 use netcdf

 implicit none

 include 'mpif.h'

 character(len=150) :: filenetcdf, raw_file

 integer*4, parameter :: isrc=43200
 integer*4, parameter :: jsrc=21600

 integer :: ncid, status, dim_lat, dim_lon, dim_lat_p1, dim_lon_p1
 integer :: i, j, jj, lon, lat, dim_time
 integer :: id_lat, id_lon, landice
 integer :: id_lat_corner, id_lon_corner
 integer :: num_times, id_times, id_data
 integer :: times, ierr, iunit
 integer :: istart, iend, jstart, jend
 integer(kind=8) :: offset
 integer(kind=1), allocatable :: global_dat(:,:)
 integer :: water_flag

 real    :: missing = -999.9
 real(kind=8)    :: lat11, lon11, dx, dy
 real(kind=8), allocatable :: lons(:), lons_corner(:)
 real(kind=8), allocatable :: lats(:), lats_corner(:)
 real, allocatable :: the_data(:,:)


 raw_file="/work/noaa/da/ggayno/save/ufs_utils.git/fv3.vegt.new.tundra.netcdf/interp/vegtype_igbp1a.30s.bin"
 filenetcdf="./vegetation_type.igbp.0.01.nc"

 call mpi_init(ierr)
 print*,'open: ',trim(raw_file)
 iunit = 9
 call mpi_file_open(mpi_comm_world, raw_file, mpi_mode_rdonly, &
                    mpi_info_null, iunit, ierr)

 if (ierr /= 0) then
   print*,'bad open', ierr
   stop
 endif

 allocate(global_dat(isrc,jsrc))
 offset = 48_8
 call mpi_file_read_at(iunit, offset, global_dat, (isrc*jsrc), &
                        mpi_integer1, mpi_status_ignore, ierr)
 if (ierr /= 0) then
   print*,'bad read ', ierr
   stop
 endif

 print*,'global_data ',maxval(global_dat),minval(global_dat)

 water_flag = 17

! conus
 istart = 5251
 iend = 14500
 jstart = 11501
 jend   = 18000

 lon = iend - istart + 1
 lat = jend - jstart + 1

 allocate(the_data(istart:iend,jstart:jend))
 do j = jstart, jend
   jj = jsrc - j + 1
 do i = istart, iend
   if (global_dat(i,jj) == water_flag) then
     the_data(i,j) = missing
   else
     the_data(i,j) = float(global_dat(i,jj))
   endif
 enddo
 enddo

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
   print*,'lon ',i,lons(i)
 enddo

 do i = jstart, jend
   lats(i) = real((i-1),8) * dy + lat11
   print*,'lat ',i,lats(i)
 enddo

! lat/lon of corner of grid box.

 allocate(lons_corner(istart:iend+1))  ! required one extra corner
                                       ! for non-periodic grids.
 allocate(lats_corner(jstart:jend+1))  ! requires one extra corner

 do i = istart, iend+1
   lons_corner(i) = (real((i-1),8) * dx + lon11) - 0.5_8*dx
   print*,'corner lon ',i,lons_corner(i)
 enddo

 do i = jstart, jend+1
   lats_corner(i) = (real((i-1),8) * dy + lat11) - 0.5_8*dy
   print*,'corner lat ',i,lats_corner(i)
 enddo

 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, ior(nf90_netcdf4, nf90_classic_model), ncid)
 if (status /= nf90_noerr) stop 1

 status=nf90_def_dim(ncid, 'idim', lon, dim_lon)
 if (status /= nf90_noerr) stop 3

 status=nf90_def_dim(ncid, 'jdim', lat, dim_lat)
 if (status /= nf90_noerr) stop 2

 status=nf90_def_dim(ncid, 'idim_p1', (lon+1), dim_lon_p1)
 if (status /= nf90_noerr) stop 3

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

 status=nf90_def_var(ncid, 'lon_corner', nf90_double, dim_lon_p1, id_lon_corner)
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

 end program soilveg
