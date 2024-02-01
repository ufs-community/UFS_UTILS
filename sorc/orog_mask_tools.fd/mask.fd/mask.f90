 program mask_netcdf

! Convert the UMD land use data to netcdf.

 use netcdf

 implicit none

 integer*4, parameter :: idim=43200
 integer*4, parameter :: jdim=21600

 character(len=150) :: filenetcdf, fileraw

 integer :: i, istat, ncid, status, dim_idim, dim_jdim
 integer :: dim_idim_p1, dim_jdim_p1
 integer :: id_lon, id_lat, id_data

 integer(kind=1), allocatable :: mask(:,:)

 real(kind=8), allocatable :: lats(:), lons(:)
 real(kind=8)    :: lat11, lon11, dx, dy

 dx = 1.0_8/120.0_8
 dy = -(1.0_8/120.0_8)

 lat11 =  90.0_8 + dy*0.5_8
 lon11 = -180.0_8 + dx*0.5_8

 allocate(lons(idim),lats(jdim),mask(idim,jdim))

 do i = 1, idim
   lons(i) = real((i-1),8) * dx + lon11
   print*,'lon ',i,lons(i)
 enddo

 do i = 1, jdim
   lats(i) = real((i-1),8) * dy + lat11
   print*,'lat ',i,lats(i)
 enddo

 fileraw="/scratch1/NCEPDEV/global/glopara/fix/raw/orog/landcover30.fixed"

 open(11, file=trim(fileraw), access='direct', recl=idim*jdim)
 read(11, rec=1, iostat=istat) mask
 if (istat /= 0) stop 99
 close(11)

 print*,'mask ', maxval(mask),minval(mask)
 where(mask > 0) mask = 1

 filenetcdf="./mask.nc"

 print*,"- CREATE FILE: ", trim(filenetcdf)
 status=nf90_create(filenetcdf, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid)
 if (status /= nf90_noerr) stop 1

 status=nf90_def_dim(ncid, 'idim', idim, dim_idim)
 if (status /= nf90_noerr) stop 3

 status=nf90_def_dim(ncid, 'jdim', jdim, dim_jdim)
 if (status /= nf90_noerr) stop 2

 status=nf90_def_dim(ncid, 'idim_p1', (idim+1), dim_idim_p1)
 if (status /= nf90_noerr) stop 4

 status=nf90_def_dim(ncid, 'jdim_p1', (jdim+1), dim_jdim_p1)
 if (status /= nf90_noerr) stop 5

 status=nf90_put_att(ncid, nf90_global, 'source', 'UMD land use')
 if (status /= nf90_noerr) stop 6

 status=nf90_def_var(ncid, 'lat', nf90_double, dim_jdim, id_lat)
 if (status /= nf90_noerr) stop 17

 status=nf90_put_att(ncid, id_lat, 'long_name', 'grid cell center latitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'lon', nf90_double, dim_idim, id_lon)
 if (status /= nf90_noerr) stop 16

 status=nf90_put_att(ncid, id_lon, 'long_name', 'grid cell center longitude')
 if (status /= nf90_noerr) stop 10

 status=nf90_def_var(ncid, 'land_mask', nf90_byte, (/dim_idim,dim_jdim/), id_data)
 if (status /= nf90_noerr) stop 20

 status=nf90_enddef(ncid)
 if (status /= nf90_noerr) stop 22

 status=nf90_put_var(ncid, id_lon, lons)
 if (status /= nf90_noerr) stop 19

 status=nf90_put_var(ncid, id_lat, lats)
 if (status /= nf90_noerr) stop 20

 status=nf90_put_var(ncid, id_data, mask)
 if (status /= nf90_noerr) stop 24

 status=nf90_close(ncid)

 end program mask_netcdf
