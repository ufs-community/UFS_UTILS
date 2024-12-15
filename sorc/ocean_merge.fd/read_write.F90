!> Read the grid dimensions from the MOM6 ocean mask NetCDF file.
!!
!! @param[in] pth1 Directory path to file.
!! @param[in] atmres Atmospheric resolution.
!! @param[in] ocnres Ocean resolution in decimal percent.
!! @param[in] tile Tile number.
!! @param[out] lon E/W dimension of tile.
!! @param[out] lat N/S dimension of tile.
!!
!! @author Shan Sun
!! @author Rahul Mahajan
 subroutine read_grid_dims(pth1, atmres, ocnres, tile, lon, lat)

 use netcdf

 implicit none

 character(len=*), intent(in)   :: pth1
 character(len=*), intent(in)   :: atmres
 character(len=*), intent(in)   :: ocnres

 integer, intent(in)            :: tile
 integer, intent(out)           :: lon, lat

 character(len=250)             :: flnm

 integer                        :: ncid, ndims, nvars, natts
 integer                        :: latid, lonid

 write(flnm,'(5a,i1,a)') trim(pth1),trim(atmres),'.',trim(ocnres),'.tile',tile,'.nc'

 print*,'- READ GRID DIMESIONS FROM: ',trim(flnm)

 call handle_err (nf90_open (flnm, NF90_NOWRITE, ncid))
 call handle_err (nf90_inquire (ncid, ndimensions=ndims, nvariables=nvars, nattributes=natts))
 call handle_err (nf90_inquire (ncid, ndimensions=ndims, nvariables=nvars, nattributes=natts))
 call handle_err (nf90_inq_dimid (ncid, 'grid_xt', latid))  ! RM: lon is no longer in this file; try grid_xt
 call handle_err (nf90_inq_dimid (ncid, 'grid_yt', lonid))  ! RM: lat is no longer in this file; try grid_yt
 call handle_err (nf90_inquire_dimension (ncid, latid, len=lat))
 call handle_err (nf90_inquire_dimension (ncid, lonid, len=lon))
 call handle_err (nf90_close (ncid))

 print*,'- DIMENSIONS ARE: ',lon, lat

 end subroutine read_grid_dims

!> Read the ocean fraction from the MOM6 ocean NetCDF file.
!!
!! @param[in] pth1 Directory path to file.
!! @param[in] atmres Atmospheric resolution.
!! @param[in] ocnres Ocean resolution.
!! @param[in] tile Tile number.
!! @param[in] lon E/W dimension of tile.
!! @param[in] lat N/S dimension of tile.
!! @param[out] ocn_frac ocean fraction in decimal percent.
!!
!! @author Shan Sun
!! @author Rahul Mahajan
 subroutine read_ocean_frac(pth1,atmres,ocnres,tile,lon,lat,ocn_frac)

 use netcdf

 implicit none

 character(len=*), intent(in)   :: pth1
 character(len=*), intent(in)   :: atmres
 character(len=*), intent(in)   :: ocnres
 
 integer,          intent(in)   :: lat
 integer,          intent(in)   :: lon
 integer,          intent(in)   :: tile

 real, intent(out)              :: ocn_frac(lon,lat)

 character(len=300)             :: flnm

 integer                        :: ncid, v1id, start(2), count(2)

 write(flnm,'(5a,i1,a)') trim(pth1),trim(atmres),'.',trim(ocnres),'.tile',tile,'.nc'

 print*,'-READ OCEAN FRACTION FROM: ',trim(flnm)

 start(1:2) = (/1,1/)
 count(1:2) = (/lon,lat/)

 call handle_err (nf90_open (flnm, NF90_NOWRITE, ncid))

! The file record is named 'land_frac', but the data is ocean fraction.

 call handle_err (nf90_inq_varid(ncid, 'land_frac', v1id))
 call handle_err (nf90_get_var (ncid, v1id, ocn_frac, start=start, count=count))
 call handle_err (nf90_close (ncid))

 end subroutine read_ocean_frac

!> Read lake fraction, lake depth and latitude from a NetCDF file.
!!
!! @param[in] pth2 Directory path to the file.
!! @param[in] atmres Atmospheric resolution.
!! @param[in] tile Tile number.
!! @param[in] lon E/W dimension of tile.
!! @param[in] lat N/S dimension of tile.
!! @param[out] lake_frac Lake fraction in decimal percent.
!! @param[out] lake_depth Lake depth in meters.
!! @param[out] lat2d Latitude in degrees.
!!
!! @author Shan Sun
!! @author Rahul Mahajan
 subroutine read_lake_mask(pth2,atmres,tile,lon,lat,lake_frac, &
                           lake_depth,lat2d)

 use netcdf

 implicit none

 character(len=*), intent(in) :: pth2, atmres

 integer, intent(in)          :: tile, lon, lat

 real, intent(out)            :: lake_frac(lon,lat)
 real, intent(out)            :: lake_depth(lon,lat)
 real, intent(out)            :: lat2d(lon,lat)

 character(len=250)           :: flnm

 integer                      :: ncid, ndims, nvars, natts
 integer                      :: v2id, v3id, vlat
 integer                      :: start(2), count(2)

 write(flnm,'(4a,i1,a)') trim(pth2),'oro.',trim(atmres),'.tile',tile,'.nc'
 print *,'- READ LAKE DEPTH, FRACTION AND LATITUDE FROM: ',trim(flnm)
 call handle_err (nf90_open (flnm, NF90_NOWRITE, ncid))
 call handle_err (nf90_inquire (ncid, ndimensions=ndims, nvariables=nvars, nattributes=natts))
 call handle_err (nf90_inq_varid(ncid, 'lake_frac', v2id))
 call handle_err (nf90_inq_varid(ncid, 'lake_depth',v3id))
 call handle_err (nf90_inq_varid(ncid, 'geolat'    ,vlat))
 start(1:2) = (/1,1/)
 count(1:2) = (/lon,lat/)
 call handle_err (nf90_get_var (ncid, v2id, lake_frac, start=start, count=count))
 call handle_err (nf90_get_var (ncid, v3id, lake_depth,start=start, count=count))
 call handle_err (nf90_get_var (ncid, vlat, lat2d,     start=start, count=count))
 call handle_err (nf90_close (ncid))

 end subroutine read_lake_mask
!> Write the merged data to a NetCDF file. Each tile is in
!! its own file.
!!
!! @param[in] atmres Atmospheric resolution.
!! @param[in] ocnres Ocean resolution.
!! @param[in] pth3 Directory path to output file.
!! @param[in] tile Tile number.
!! @param[in] lon E/W dimension of tile.
!! @param[in] lat N/S dimension of tile.
!! @param[in] land_frac Land fraction in decimal percent.
!! @param[in] lake_frac Lake fraction.
!! @param[in] lake_depth Lake depth in meters.
!! @param[in] slmsk Land/sea mask - 0-non-land; 1-land.
!!
!! @author Shan Sun
!! @author Rahul Mahajan
!! @author Sanath Kumar
 subroutine write_data(atmres,ocnres,pth3,tile,lon,lat,land_frac, &
                       lake_frac,lake_depth,slmsk)

 use netcdf

 implicit none

 character(len=*), intent(in)       :: atmres, ocnres, pth3
 
 integer, intent(in)                :: tile, lon, lat

 real, intent(in)                   :: land_frac(lon,lat), lake_frac(lon,lat)
 real, intent(in)                   :: lake_depth(lon,lat), slmsk(lon,lat)

 character(len=250) :: flnm

 integer            :: ncid4, dims(2), v1id, v2id, v3id, v4id
 
 write(flnm,'(4a,i1,a)') trim(atmres),'.',trim(ocnres),'.tile',tile,'.nc'
 print *,'- OUTPUT DATA TO FILE: ',trim(flnm)
 call handle_err (nf90_create (path=trim(pth3)//trim(flnm), &
    cmode=or(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid=ncid4))   ! netcdf3

 call handle_err (nf90_def_dim (ncid4,'lon', lon, dims(1)))
 call handle_err (nf90_def_dim (ncid4,'lat', lat, dims(2)))
 call handle_err (nf90_def_var (ncid4,'land_frac', nf90_float, dims(1:2), v1id))
 call handle_err (nf90_def_var (ncid4,'lake_frac', nf90_float, dims(1:2), v2id))
 call handle_err (nf90_def_var (ncid4,'lake_depth',nf90_float, dims(1:2), v3id))
 call handle_err (nf90_def_var (ncid4,'slmsk',     nf90_float, dims(1:2), v4id))

 call handle_err (nf90_enddef  (ncid4))

 call handle_err (nf90_put_var (ncid4, v1id,land_frac))
 call handle_err (nf90_put_var (ncid4, v2id,lake_frac))
 call handle_err (nf90_put_var (ncid4, v3id,lake_depth))
 call handle_err (nf90_put_var (ncid4, v4id,slmsk))
 call handle_err (nf90_close(ncid4))

 end subroutine write_data
