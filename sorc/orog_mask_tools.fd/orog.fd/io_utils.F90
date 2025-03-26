!> @file
!! @brief i/o utilities
!! @author George Gayno NOAA/EMC

!> Module containing utilities that read and write data.
!!
!! @author George Gayno NOAA/EMC

 module io_utils

 implicit none

 private

 public :: qc_orog_by_ramp
 public :: read_global_mask
 public :: read_global_orog
 public :: read_mask
 public :: read_mdl_dims
 public :: read_mdl_grid_file
 public :: write_mask_netcdf
 public :: write_netcdf

 contains

!> Write out orography file in netcdf format.
!!
!! @param[in] im 'i' dimension of a model grid tile.
!! @param[in] jm 'j' dimension of a model grid tile.
!! @param[in] slm Land-sea mask.
!! @param[in] land_frac Land fraction.
!! @param[in] oro Orography
!! @param[in] hprime The gravity wave drag fields on the model grid tile.
!! @param[in] ntiles Number of tiles to output.
!! @param[in] tile Tile number to output.
!! @param[in] geolon Longitude on the model grid tile.
!! @param[in] geolat Latitude on the model grid tile.
!! @param[in] lon Longitude of the first row of the model grid tile.
!! @param[in] lat Latitude of the first column of the model grid tile.
!! @author Jordan Alpert NOAA/EMC GFDL Programmer
  subroutine write_netcdf(im, jm, slm, land_frac, oro, hprime, ntiles, tile, geolon, geolat, lon, lat)
    use netcdf
    implicit none
    integer, intent(in):: im, jm, ntiles, tile
    real, intent(in) :: lon(im), lat(jm)
    real, intent(in), dimension(im,jm)  :: slm, oro, geolon, geolat, land_frac
    real, intent(in), dimension(im,jm,14):: hprime
    character(len=128) :: outfile
    integer            :: error, ncid
    integer            :: header_buffer_val = 16384      
    integer            :: fsize=65536, inital = 0  
    integer            :: dim1, dim2
    integer            :: dim_lon, dim_lat
    integer            :: id_geolon,id_geolat
    integer            :: id_slmsk,id_orog_raw,id_orog_filt,id_land_frac
    integer            :: id_stddev,id_convex
    integer            :: id_oa1,id_oa2,id_oa3,id_oa4
    integer            :: id_ol1,id_ol2,id_ol3,id_ol4
    integer            :: id_theta,id_gamma,id_sigma,id_elvmax

    if(ntiles > 1) then
      write(outfile, '(a,i4.4,a)') 'out.oro.tile', tile, '.nc'
    else
      outfile = "out.oro.nc"
    endif

    dim1=size(lon,1)
    dim2=size(lat,1)
      
    !--- open the file
    error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, &
            initialsize=inital, chunksize=fsize)
    call netcdf_err(error, 'Creating file '//trim(outfile) )
    !--- define dimension
    error = nf90_def_dim(ncid, 'lon', dim1, dim_lon)
    call netcdf_err(error, 'define dimension lon for file='//trim(outfile) )
    error = nf90_def_dim(ncid, 'lat', dim2, dim_lat)
    call netcdf_err(error, 'define dimension lat for file='//trim(outfile) )  

    !--- define field
!---geolon
    error = nf90_def_var(ncid, 'geolon', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolon)
    call netcdf_err(error, 'define var geolon for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolon, "long_name", "Longitude")
    call netcdf_err(error, 'define geolon name for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolon, "units", "degrees_east")
    call netcdf_err(error, 'define geolon units for file='//trim(outfile) )
!---geolat
    error = nf90_def_var(ncid, 'geolat', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolat)
    call netcdf_err(error, 'define var geolat for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolat, "long_name", "Latitude")
    call netcdf_err(error, 'define geolat name for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolat, "units", "degrees_north")
    call netcdf_err(error, 'define geolat units for file='//trim(outfile) )
!---slmsk
    error = nf90_def_var(ncid, 'slmsk', NF90_FLOAT,(/dim_lon,dim_lat/), id_slmsk)
    call netcdf_err(error, 'define var slmsk for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_slmsk, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define slmsk coordinates for file='//trim(outfile) )
!--- land_frac
    error = nf90_def_var(ncid, 'land_frac', NF90_FLOAT, (/dim_lon,dim_lat/), id_land_frac)
    call netcdf_err(error, 'define var land_frac for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_land_frac, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define land_frac coordinates for file='//trim(outfile) )
!---orography - raw
    error = nf90_def_var(ncid, 'orog_raw', NF90_FLOAT, (/dim_lon,dim_lat/), id_orog_raw)
    call netcdf_err(error, 'define var orog_raw for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_orog_raw, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define orog_raw coordinates for file='//trim(outfile) )
!---orography - filtered
    error = nf90_def_var(ncid, 'orog_filt', NF90_FLOAT, (/dim_lon,dim_lat/), id_orog_filt)
    call netcdf_err(error, 'define var orog_filt for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_orog_filt, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define orog_filt coordinates for file='//trim(outfile) )
!---stddev
    error = nf90_def_var(ncid, 'stddev', NF90_FLOAT, (/dim_lon,dim_lat/), id_stddev)
    call netcdf_err(error, 'define var stddev for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_stddev, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define stddev coordinates for file='//trim(outfile) )
!---convexity
    error = nf90_def_var(ncid, 'convexity', NF90_FLOAT, (/dim_lon,dim_lat/), id_convex)
    call netcdf_err(error, 'define var convexity for file='//trim(outfile) )      
    error = nf90_put_att(ncid, id_convex, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define convexity coordinates for file='//trim(outfile) )
!---oa1 -> oa4
    error = nf90_def_var(ncid, 'oa1', NF90_FLOAT, (/dim_lon,dim_lat/), id_oa1)
    call netcdf_err(error, 'define var oa1 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_oa1, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define oa1 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'oa2', NF90_FLOAT, (/dim_lon,dim_lat/), id_oa2)
    call netcdf_err(error, 'define var oa2 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_oa2, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define oa2 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'oa3', NF90_FLOAT, (/dim_lon,dim_lat/), id_oa3)
    call netcdf_err(error, 'define var oa3 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_oa3, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define oa3 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'oa4', NF90_FLOAT, (/dim_lon,dim_lat/), id_oa4)
    call netcdf_err(error, 'define var oa4 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_oa4, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define oa4 coordinates for file='//trim(outfile) )
!---ol1 -> ol4
    error = nf90_def_var(ncid, 'ol1', NF90_FLOAT, (/dim_lon,dim_lat/), id_ol1)
    call netcdf_err(error, 'define var ol1 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_ol1, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define ol1 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'ol2', NF90_FLOAT, (/dim_lon,dim_lat/), id_ol2)
    call netcdf_err(error, 'define var ol2 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_ol2, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define ol2 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'ol3', NF90_FLOAT, (/dim_lon,dim_lat/), id_ol3)
    call netcdf_err(error, 'define var ol3 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_ol3, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define ol3 coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'ol4', NF90_FLOAT, (/dim_lon,dim_lat/), id_ol4)
    call netcdf_err(error, 'define var ol4 for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_ol4, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define ol4 coordinates for file='//trim(outfile) )
!---theta gamma sigma elvmax
    error = nf90_def_var(ncid, 'theta', NF90_FLOAT, (/dim_lon,dim_lat/), id_theta)
    call netcdf_err(error, 'define var theta for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_theta, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define theta coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'gamma', NF90_FLOAT, (/dim_lon,dim_lat/), id_gamma)
    call netcdf_err(error, 'define var gamma for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_gamma, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define gamma coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'sigma', NF90_FLOAT, (/dim_lon,dim_lat/), id_sigma)
    call netcdf_err(error, 'define var sigma for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_sigma, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define sigma coordinates for file='//trim(outfile) )
    error = nf90_def_var(ncid, 'elvmax', NF90_FLOAT, (/dim_lon,dim_lat/), id_elvmax)
    call netcdf_err(error, 'define var elvmax for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_elvmax, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define elvmax coordinates for file='//trim(outfile) )

    error = nf90_enddef(ncid, header_buffer_val,4,0,4)
    call netcdf_err(error, 'end meta define for file='//trim(outfile) )
      
    !--- write out data
    error = nf90_put_var( ncid, id_geolon, geolon(:dim1,:dim2))
    call netcdf_err(error, 'write var geolon for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_geolat, geolat(:dim1,:dim2))
    call netcdf_err(error, 'write var geolat for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_slmsk, slm(:dim1,:dim2))
    call netcdf_err(error, 'write var slmsk for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_land_frac, land_frac(:dim1,:dim2))
    call netcdf_err(error, 'write var land_frac for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_orog_raw, oro(:dim1,:dim2))
    call netcdf_err(error, 'write var orog_raw for file='//trim(outfile) )
! We no longer filter the orog, so the raw and filtered records are the same.
    error = nf90_put_var( ncid, id_orog_filt, oro(:dim1,:dim2))
    call netcdf_err(error, 'write var orog_filt for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_stddev, hprime(:dim1,:dim2,1))
    call netcdf_err(error, 'write var stddev for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_convex, hprime(:dim1,:dim2,2))
    call netcdf_err(error, 'write var convex for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_oa1, hprime(:dim1,:dim2,3))
    call netcdf_err(error, 'write var oa1 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_oa2, hprime(:dim1,:dim2,4))
    call netcdf_err(error, 'write var oa2 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_oa3, hprime(:dim1,:dim2,5))
    call netcdf_err(error, 'write var oa3 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_oa4, hprime(:dim1,:dim2,6))
    call netcdf_err(error, 'write var oa4 for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_ol1, hprime(:dim1,:dim2,7))
    call netcdf_err(error, 'write var ol1 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_ol2, hprime(:dim1,:dim2,8))
    call netcdf_err(error, 'write var ol2 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_ol3, hprime(:dim1,:dim2,9))
    call netcdf_err(error, 'write var ol3 for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_ol4, hprime(:dim1,:dim2,10))
    call netcdf_err(error, 'write var ol4 for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_theta, hprime(:dim1,:dim2,11))
    call netcdf_err(error, 'write var theta for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_gamma, hprime(:dim1,:dim2,12))
    call netcdf_err(error, 'write var gamma for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_sigma, hprime(:dim1,:dim2,13))
    call netcdf_err(error, 'write var sigma for file='//trim(outfile) )
    error = nf90_put_var( ncid, id_elvmax, hprime(:dim1,:dim2,14))
    call netcdf_err(error, 'write var elvmax for file='//trim(outfile) )

    error = nf90_close(ncid) 
    call netcdf_err(error, 'close file='//trim(outfile) )  
      
  end subroutine write_netcdf

!> Check NetCDF error code and output the error message.
!!
!! @param[in] err NetCDF error code
!! @param[in] string The NetCDF error message
!! @author Jordan Alpert NOAA/EMC
  subroutine netcdf_err( err, string )
      use netcdf
      implicit none
      integer, intent(in) :: err
      character(len=*), intent(in) :: string
      character(len=256) :: errmsg

      if( err.EQ.NF90_NOERR )return
      errmsg = NF90_STRERROR(err)
      print*, 'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
      call abort

      return
  end subroutine netcdf_err

!> Write the land mask file
!!
!! @param[in] im 'i' dimension of a model grid tile.
!! @param[in] jm 'j' dimension of a model grid tile.
!! @param[in] slm Land-sea mask.
!! @param[in] land_frac Land fraction.
!! @param[in] ntiles Number of tiles to output.
!! @param[in] tile Tile number to output.
!! @param[in] geolon Longitude on the model grid tile.
!! @param[in] geolat Latitude on the model grid tile.
!! @author George Gayno NOAA/EMC

  subroutine write_mask_netcdf(im, jm, slm, land_frac, ntiles, tile, geolon, geolat)
    use netcdf
    implicit none
    integer, intent(in):: im, jm, ntiles, tile
    real, intent(in), dimension(im,jm)  :: slm, geolon, geolat, land_frac
    character(len=128) :: outfile
    integer            :: error, ncid
    integer            :: header_buffer_val = 16384      
    integer            :: fsize=65536, inital = 0  
    integer            :: dim1, dim2
    integer            :: dim_lon, dim_lat
    integer            :: id_geolon,id_geolat
    integer            :: id_slmsk,id_land_frac

    if(ntiles > 1) then
      write(outfile, '(a,i4.4,a)') 'out.oro.tile', tile, '.nc'
    else
      outfile = "out.oro.nc"
    endif

    dim1=im
    dim2=jm
      
    !--- open the file
    error = nf90_create(outfile, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, &
            initialsize=inital, chunksize=fsize)
    call netcdf_err(error, 'Creating file '//trim(outfile) )
    !--- define dimension
    error = nf90_def_dim(ncid, 'lon', dim1, dim_lon)
    call netcdf_err(error, 'define dimension lon for file='//trim(outfile) )
    error = nf90_def_dim(ncid, 'lat', dim2, dim_lat)
    call netcdf_err(error, 'define dimension lat for file='//trim(outfile) )  

    !--- define field
!---geolon
    error = nf90_def_var(ncid, 'geolon', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolon)
    call netcdf_err(error, 'define var geolon for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolon, "long_name", "Longitude")
    call netcdf_err(error, 'define geolon name for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolon, "units", "degrees_east")
    call netcdf_err(error, 'define geolon units for file='//trim(outfile) )
!---geolat
    error = nf90_def_var(ncid, 'geolat', NF90_FLOAT, (/dim_lon,dim_lat/), id_geolat)
    call netcdf_err(error, 'define var geolat for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolat, "long_name", "Latitude")
    call netcdf_err(error, 'define geolat name for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_geolat, "units", "degrees_north")
    call netcdf_err(error, 'define geolat units for file='//trim(outfile) )
!---slmsk
    error = nf90_def_var(ncid, 'slmsk', NF90_FLOAT, (/dim_lon,dim_lat/), id_slmsk)
    call netcdf_err(error, 'define var slmsk for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_slmsk, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define slmsk coordinates for file='//trim(outfile) )
!--- land_frac
    error = nf90_def_var(ncid, 'land_frac', NF90_FLOAT, (/dim_lon,dim_lat/), id_land_frac)
    call netcdf_err(error, 'define var land_frac for file='//trim(outfile) )
    error = nf90_put_att(ncid, id_land_frac, "coordinates", "geolon geolat")
    call netcdf_err(error, 'define land_frac coordinates for file='//trim(outfile) )

    error = nf90_enddef(ncid, header_buffer_val,4,0,4)
    call netcdf_err(error, 'end meta define for file='//trim(outfile) )
      
    !--- write out data
    error = nf90_put_var( ncid, id_geolon, geolon(:dim1,:dim2))
    call netcdf_err(error, 'write var geolon for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_geolat, geolat(:dim1,:dim2))
    call netcdf_err(error, 'write var geolat for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_slmsk, slm(:dim1,:dim2))
    call netcdf_err(error, 'write var slmsk for file='//trim(outfile) )

    error = nf90_put_var( ncid, id_land_frac, land_frac(:dim1,:dim2))
    call netcdf_err(error, 'write var land_frac for file='//trim(outfile) )

    error = nf90_close(ncid) 
    call netcdf_err(error, 'close file='//trim(outfile) )  
      
  end subroutine write_mask_netcdf
 
!> Read the land mask file
!!
!! @param[in] merge_file path 
!! @param[out] slm Land-sea mask.
!! @param[out] land_frac Land fraction.
!! @param[out] lake_frac Lake fraction
!! @param[in] im 'i' dimension of a model grid tile.
!! @param[in] jm 'j' dimension of a model grid tile.
!! @author George Gayno NOAA/EMC

  subroutine read_mask(merge_file,slm,land_frac,lake_frac,im,jm)

  use netcdf

  implicit none

  character(len=*), intent(in) :: merge_file

  integer, intent(in) :: im, jm

  real, intent(out) :: land_frac(im,jm)
  real, intent(out) :: lake_frac(im,jm)
  real, intent(out) :: slm(im,jm)

  integer ncid, error, id_var

  print*,'- READ IN EXTERNAL LANDMASK FILE: ',trim(merge_file)
  error=nf90_open(merge_file,nf90_nowrite,ncid)
  call netcdf_err(error, 'Open file '//trim(merge_file) )

  error=nf90_inq_varid(ncid, 'land_frac', id_var)
  call netcdf_err(error, 'inquire varid of land_frac')
  error=nf90_get_var(ncid, id_var, land_frac)
  call netcdf_err(error, 'inquire data of land_frac')

  error=nf90_inq_varid(ncid, 'slmsk', id_var)
  call netcdf_err(error, 'inquire varid of slmsk')
  error=nf90_get_var(ncid, id_var, slm)
  call netcdf_err(error, 'inquire data of slmsk')

  error=nf90_inq_varid(ncid, 'lake_frac', id_var)
  call netcdf_err(error, 'inquire varid of lake_frac')
  error=nf90_get_var(ncid, id_var, lake_frac)
  call netcdf_err(error, 'inquire data of lake_frac')

  error = nf90_close(ncid) 

  end subroutine read_mask

!> Read the grid dimensions from the model 'grid' file
!!
!! @param[in] mdl_grid_file path/name of model 'grid' file.
!! @param[out] im 'i' dimension of a model grid tile.
!! @param[out] jm 'j' dimension of a model grid tile.
!! @author George Gayno NOAA/EMC
  subroutine read_mdl_dims(mdl_grid_file, im, jm)

  use netcdf

  implicit none

  character(len=*), intent(in) :: mdl_grid_file

  integer, intent(out)         :: im, jm

  integer ncid, error, id_dim, nx, ny

  print*, "- READ MDL GRID DIMENSIONS FROM= ", trim(mdl_grid_file)

  error=nf90_open(mdl_grid_file, nf90_nowrite, ncid)
  call netcdf_err(error, 'Opening file '//trim(mdl_grid_file) )

  error=nf90_inq_dimid(ncid, 'nx', id_dim)
  call netcdf_err(error, 'inquire dimension nx from file '// trim(mdl_grid_file) )
  error=nf90_inquire_dimension(ncid, id_dim, len=nx)
  call netcdf_err(error, 'inquire nx from file '//trim(mdl_grid_file) )

  error=nf90_inq_dimid(ncid, 'ny', id_dim)
  call netcdf_err(error, 'inquire dimension ny from file '// trim(mdl_grid_file) )
  error=nf90_inquire_dimension(ncid, id_dim, len=ny)
  call netcdf_err(error, 'inquire ny from file '//trim(mdl_grid_file) )

  error=nf90_close(ncid)

  IM = nx/2
  JM = ny/2

  print*,"- MDL GRID DIMENSIONS ", im, jm

  end subroutine read_mdl_dims

!> Read the grid dimensions from the model 'grid' file
!!
!! @param[in] mdl_grid_file Path/name of model 'grid' file.
!! @param[in] im 'i' Dimension of a model grid tile.
!! @param[in] jm 'j' Dimension of a model grid tile.
!! @param[out] geolon Longitude at the grid point centers.
!! @param[out] geolon_c Longitude at the grid point corners.
!! @param[out] geolat Latitude at the grid point centers.
!! @param[out] geolat_c Latitude at the grid point corners.
!! @param[out] dx Length of model grid points in the 'x' direction.
!! @param[out] dy Length of model grid points in the 'y' direction.
!! @param[out] is_north_pole 'true' for points surrounding the north pole.
!! @param[out] is_south_pole 'true' for points surrounding the south pole.
!! @author George Gayno NOAA/EMC
  subroutine read_mdl_grid_file(mdl_grid_file, im, jm, &
             geolon, geolon_c, geolat, geolat_c, dx, dy, &
             is_north_pole, is_south_pole)

  use netcdf

  use orog_utils, only : find_poles, find_nearest_pole_points

  implicit none

  character(len=*), intent(in) :: mdl_grid_file

  integer, intent(in)          :: im, jm

  logical, intent(out)         :: is_north_pole(im,jm)
  logical, intent(out)         :: is_south_pole(im,jm)

  real, intent(out)            :: geolat(im,jm)
  real, intent(out)            :: geolat_c(im+1,jm+1)
  real, intent(out)            :: geolon(im,jm)
  real, intent(out)            :: geolon_c(im+1,jm+1)
  real, intent(out)            :: dx(im,jm), dy(im,jm)

  integer                      :: i, j
  integer                      :: ncid, error, id_var, nx, ny
  integer                      :: i_south_pole,j_south_pole
  integer                      :: i_north_pole,j_north_pole

  real, allocatable     :: tmpvar(:,:)

  nx = 2*im
  ny = 2*jm

  allocate(tmpvar(nx+1,ny+1))

  print*, "- OPEN AND READ= ", trim(mdl_grid_file)

  error=nf90_open(mdl_grid_file, nf90_nowrite, ncid)
  call netcdf_err(error, 'Opening file '//trim(mdl_grid_file) )

  error=nf90_inq_varid(ncid, 'x', id_var)
  call netcdf_err(error, 'inquire varid of x from file ' // trim(mdl_grid_file))
  error=nf90_get_var(ncid, id_var, tmpvar)
  call netcdf_err(error, 'inquire data of x from file ' // trim(mdl_grid_file))

! Adjust lontitude to be between 0 and 360.
  do j = 1,ny+1
  do i = 1,nx+1
    if(tmpvar(i,j) .GT. 360) tmpvar(i,j) = tmpvar(i,j) - 360
    if(tmpvar(i,j) .LT. 0) tmpvar(i,j) = tmpvar(i,j) + 360
  enddo
  enddo

  geolon(1:IM,1:JM) = tmpvar(2:nx:2,2:ny:2)
  geolon_c(1:IM+1,1:JM+1) = tmpvar(1:nx+1:2,1:ny+1:2)

  error=nf90_inq_varid(ncid, 'y', id_var)
  call netcdf_err(error, 'inquire varid of y from file ' // trim(mdl_grid_file))
  error=nf90_get_var(ncid, id_var, tmpvar)
  call netcdf_err(error, 'inquire data of y from file ' // trim(mdl_grid_file))

  geolat(1:IM,1:JM) = tmpvar(2:nx:2,2:ny:2)
  geolat_c(1:IM+1,1:JM+1) = tmpvar(1:nx+1:2,1:ny+1:2)

  call find_poles(tmpvar, nx, ny, i_north_pole, j_north_pole, &
                  i_south_pole, j_south_pole)

  deallocate(tmpvar)

  call find_nearest_pole_points(i_north_pole, j_north_pole, &
       i_south_pole, j_south_pole, im, jm, is_north_pole, &
       is_south_pole)

  allocate(tmpvar(nx,ny))

  error=nf90_inq_varid(ncid, 'area', id_var)
  call netcdf_err(error, 'inquire varid of area from file ' // trim(mdl_grid_file))
  error=nf90_get_var(ncid, id_var, tmpvar)
  call netcdf_err(error, 'inquire data of area from file ' // trim(mdl_grid_file))

  error = nf90_close(ncid)

  do j = 1, jm
    do i = 1, im
      dx(i,j) = sqrt(tmpvar(2*i-1,2*j-1)+tmpvar(2*i,2*j-1)   &
                + tmpvar(2*i-1,2*j  )+tmpvar(2*i,2*j  ))
      dy(i,j) = dx(i,j)
    enddo
  enddo

  deallocate(tmpvar)

  end subroutine read_mdl_grid_file

!> Read input global 30-arc second orography data.
!!
!! @param[in] imn i-dimension of orography data.
!! @param[in] jmn j-dimension of orography data.
!! @param[out] glob The orography data.
!! @author Jordan Alpert NOAA/EMC
 subroutine read_global_orog(imn,jmn,glob)

 use orog_utils, only : transpose_orog
 use netcdf

 implicit none

 integer, intent(in)    :: imn, jmn
 integer*2, intent(out) :: glob(imn,jmn)

 integer :: ncid, error, id_dim, id_var, idim, jdim

 print*,"- OPEN AND READ ./topography.gmted2010.30s.nc"

 error=nf90_open("./topography.gmted2010.30s.nc", &
                nf90_nowrite, ncid)
 call netcdf_err(error, 'Opening file topography.gmted2010.30s.nc' )

 error=nf90_inq_dimid(ncid, 'idim', id_dim)
 call netcdf_err(error, 'Inquire dimid of idim' )

 error=nf90_inquire_dimension(ncid,id_dim,len=idim)
 call netcdf_err(error, 'Reading idim' )

 if (imn /= idim) then
   print*,"FATAL ERROR: i-dimensions do not match."
 endif

 error=nf90_inq_dimid(ncid, 'jdim', id_dim)
 call netcdf_err(error, 'Inquire dimid of jdim' )

 error=nf90_inquire_dimension(ncid,id_dim,len=jdim)
 call netcdf_err(error, 'Reading jdim' )

 if (jmn /= jdim) then
   print*,"FATAL ERROR: j-dimensions do not match."
 endif

 error=nf90_inq_varid(ncid, 'topo', id_var)
 call netcdf_err(error, 'Inquire varid of topo')

 error=nf90_get_var(ncid, id_var, glob)
 call netcdf_err(error, 'Reading topo')

 error = nf90_close(ncid)

 print*,"- MAX/MIN OF OROGRAPHY DATA ",maxval(glob),minval(glob)

 call transpose_orog(imn,jmn,glob)

 return
 end subroutine read_global_orog

!> Read input global 30-arc second land mask data.
!!
!! @param[in] imn i-dimension of orography data.
!! @param[in] jmn j-dimension of orography data.
!! @param[out] mask The land mask data.
!! @author G. Gayno NOAA/EMC
 subroutine read_global_mask(imn, jmn, mask)

 use orog_utils, only : transpose_mask
 use netcdf

 implicit none

 integer, intent(in)        :: imn, jmn

 integer(1), intent(out)    :: mask(imn,jmn)

 integer   :: ncid, id_var, id_dim, error, idim, jdim

 print*,"- OPEN AND READ ./landcover.umd.30s.nc"

 error=nf90_open("./landcover.umd.30s.nc",nf90_nowrite,ncid)
 call netcdf_err(error, 'Opening file landcover.umd.30s.nc' )

 error=nf90_inq_dimid(ncid, 'idim', id_dim)
 call netcdf_err(error, 'Inquire dimid of idim' )

 error=nf90_inquire_dimension(ncid,id_dim,len=idim)
 call netcdf_err(error, 'Reading idim' )

 if (imn /= idim) then
   print*,"FATAL ERROR: i-dimensions do not match."
 endif

 error=nf90_inq_dimid(ncid, 'jdim', id_dim)
 call netcdf_err(error, 'Inquire dimid of jdim' )

 error=nf90_inquire_dimension(ncid,id_dim,len=jdim)
 call netcdf_err(error, 'Reading jdim' )

 if (jmn /= jdim) then
   print*,"FATAL ERROR: j-dimensions do not match."
 endif

 error=nf90_inq_varid(ncid, 'land_mask', id_var)
 call netcdf_err(error, 'Inquire varid of land_mask')

 error=nf90_get_var(ncid, id_var, mask)
 call netcdf_err(error, 'Inquire data of land_mask')

 error = nf90_close(ncid)

 call transpose_mask(imn,jmn,mask)

 end subroutine read_global_mask

!> Quality control the global orography and landmask
!! data over Antarctica using RAMP data.
!!
!! @param[in] imn i-dimension of the global data.
!! @param[in] jmn j-dimension of the global data.
!! @param[inout] zavg The global orography data.
!! @param[inout] zslm The global landmask data.
!! @author G. Gayno
 subroutine qc_orog_by_ramp(imn, jmn, zavg, zslm)

 use netcdf

 implicit none

 integer, intent(in)      :: imn, jmn
 integer, intent(inout)   :: zavg(imn,jmn)
 integer, intent(inout)   :: zslm(imn,jmn)

 integer                  :: i, j, error, ncid, id_var, id_dim, jramp

 real(4), allocatable     :: gice(:,:)

! Read 30-sec Antarctica RAMP data. Points scan from South
! to North, and from Greenwich to Greenwich.

 print*,"- OPEN/READ RAMP DATA ./topography.antarctica.ramp.30s.nc"

 error=nf90_open("./topography.antarctica.ramp.30s.nc", &
                 nf90_nowrite, ncid)
 call netcdf_err(error, 'Opening RAMP topo file' )

 error=nf90_inq_dimid(ncid, 'jdim', id_dim)
 call netcdf_err(error, 'Inquire dimid of jdim' )

 error=nf90_inquire_dimension(ncid, id_dim, len=jramp)
 call netcdf_err(error, 'Reading jdim' )

 allocate (GICE(IMN+1,jramp))

 error=nf90_inq_varid(ncid, 'topo', id_var)
 call netcdf_err(error, 'Inquire varid of RAMP topo')

 error=nf90_get_var(ncid, id_var, GICE)
 call netcdf_err(error, 'Inquire data of RAMP topo')

 error = nf90_close(ncid)

 print*,"- QC GLOBAL OROGRAPHY DATA WITH RAMP."

! If RAMP values are valid, replace the global value with the RAMP
! value. Invalid values are less than or equal to 0 (0, -1, or -99).

 do j = 1, jramp
 do i = 1, IMN
   if ( GICE(i,j) .gt. 0.) then
     ZAVG(i,j) = int( GICE(i,j) + 0.5 )
     ZSLM(i,j) = 1
   endif
 enddo
 enddo

 deallocate (GICE)

 end subroutine qc_orog_by_ramp

 end module io_utils
