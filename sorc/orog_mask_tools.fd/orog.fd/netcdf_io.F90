!> @file
!! @brief Write out data in netcdf format
!! @author Jordan Alpert NOAA/EMC

!> Write out orography file in netcdf format.
!!
!! @param[in] im 'i' dimension of a model grid tile.
!! @param[in] jm 'j' dimension of a model grid tile.
!! @param[in] slm Land-sea mask.
!! @param[in] land_frac Land fraction.
!! @param[in] oro Orography
!! @param[in] orf Filtered orography. Currently the same as 'oro'.
!! @param[in] hprime The gravity wave drag fields on the model grid tile.
!! @param[in] ntiles Number of tiles to output.
!! @param[in] tile Tile number to output.
!! @param[in] geolon Longitude on the model grid tile.
!! @param[in] geolat Latitude on the model grid tile.
!! @param[in] lon Longitude of the first row of the model grid tile.
!! @param[in] lat Latitude of the first column of the model grid tile.
!! @author Jordan Alpert NOAA/EMC GFDL Programmer
  subroutine write_netcdf(im, jm, slm, land_frac, oro, orf, hprime, ntiles, tile, geolon, geolat, lon, lat)
    implicit none
    integer, intent(in):: im, jm, ntiles, tile
    real, intent(in) :: lon(im), lat(jm)
    real, intent(in), dimension(im,jm)  :: slm, oro, orf, geolon, geolat, land_frac
    real, intent(in), dimension(im,jm,14):: hprime
    character(len=128) :: outfile
    integer            :: error, ncid, i
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
    include "netcdf.inc"

    if(ntiles > 1) then
      write(outfile, '(a,i4.4,a)') 'out.oro.tile', tile, '.nc'
    else
      outfile = "out.oro.nc"
    endif

    dim1=size(lon,1)
    dim2=size(lat,1)
    write(6,*) ' netcdf dims are: ',dim1, dim2
      
    !--- open the file
    error = NF__CREATE(outfile, IOR(NF_NETCDF4,NF_CLASSIC_MODEL), inital, fsize, ncid)
    call netcdf_err(error, 'Creating file '//trim(outfile) )
    !--- define dimension
    error = nf_def_dim(ncid, 'lon', dim1, dim_lon)
    call netcdf_err(error, 'define dimension lon for file='//trim(outfile) )
    error = nf_def_dim(ncid, 'lat', dim2, dim_lat)
    call netcdf_err(error, 'define dimension lat for file='//trim(outfile) )  

    !--- define field
!---geolon
    error = nf_def_var(ncid, 'geolon', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_geolon)
    call netcdf_err(error, 'define var geolon for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_geolon, "long_name", 9, "Longitude")
    call netcdf_err(error, 'define geolon name for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_geolon, "units", 12, "degrees_east")
    call netcdf_err(error, 'define geolon units for file='//trim(outfile) )
!---geolat
    error = nf_def_var(ncid, 'geolat', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_geolat)
    call netcdf_err(error, 'define var geolat for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_geolat, "long_name", 8, "Latitude")
    call netcdf_err(error, 'define geolat name for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_geolat, "units", 13, "degrees_north")
    call netcdf_err(error, 'define geolat units for file='//trim(outfile) )
!---slmsk
    error = nf_def_var(ncid, 'slmsk', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_slmsk)
    call netcdf_err(error, 'define var slmsk for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_slmsk, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define slmsk coordinates for file='//trim(outfile) )
!--- land_frac
    error = nf_def_var(ncid, 'land_frac', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_land_frac)
    call netcdf_err(error, 'define var land_frac for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_land_frac, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define land_frac coordinates for file='//trim(outfile) )
!---orography - raw
    error = nf_def_var(ncid, 'orog_raw', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_orog_raw)
    call netcdf_err(error, 'define var orog_raw for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_orog_raw, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define orog_raw coordinates for file='//trim(outfile) )
!---orography - filtered
    error = nf_def_var(ncid, 'orog_filt', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_orog_filt)
    call netcdf_err(error, 'define var orog_filt for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_orog_filt, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define orog_filt coordinates for file='//trim(outfile) )
!---stddev
    error = nf_def_var(ncid, 'stddev', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_stddev)
    call netcdf_err(error, 'define var stddev for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_stddev, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define stddev coordinates for file='//trim(outfile) )
!---convexity
    error = nf_def_var(ncid, 'convexity', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_convex)
    call netcdf_err(error, 'define var convexity for file='//trim(outfile) )      
    error = nf_put_att_text(ncid, id_convex, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define convexity coordinates for file='//trim(outfile) )
!---oa1 -> oa4
    error = nf_def_var(ncid, 'oa1', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_oa1)
    call netcdf_err(error, 'define var oa1 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_oa1, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define oa1 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'oa2', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_oa2)
    call netcdf_err(error, 'define var oa2 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_oa2, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define oa2 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'oa3', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_oa3)
    call netcdf_err(error, 'define var oa3 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_oa3, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define oa3 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'oa4', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_oa4)
    call netcdf_err(error, 'define var oa4 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_oa4, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define oa4 coordinates for file='//trim(outfile) )
!---ol1 -> ol4
    error = nf_def_var(ncid, 'ol1', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_ol1)
    call netcdf_err(error, 'define var ol1 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_ol1, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define ol1 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'ol2', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_ol2)
    call netcdf_err(error, 'define var ol2 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_ol2, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define ol2 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'ol3', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_ol3)
    call netcdf_err(error, 'define var ol3 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_ol3, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define ol3 coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'ol4', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_ol4)
    call netcdf_err(error, 'define var ol4 for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_ol4, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define ol4 coordinates for file='//trim(outfile) )
!---theta gamma sigma elvmax
    error = nf_def_var(ncid, 'theta', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_theta)
    call netcdf_err(error, 'define var theta for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_theta, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define theta coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'gamma', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_gamma)
    call netcdf_err(error, 'define var gamma for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_gamma, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define gamma coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'sigma', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_sigma)
    call netcdf_err(error, 'define var sigma for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_sigma, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define sigma coordinates for file='//trim(outfile) )
    error = nf_def_var(ncid, 'elvmax', NF_FLOAT, 2, (/dim_lon,dim_lat/), id_elvmax)
    call netcdf_err(error, 'define var elvmax for file='//trim(outfile) )
    error = nf_put_att_text(ncid, id_elvmax, "coordinates", 13, "geolon geolat")
    call netcdf_err(error, 'define elvmax coordinates for file='//trim(outfile) )

    error = nf__enddef(ncid, header_buffer_val,4,0,4)
    call netcdf_err(error, 'end meta define for file='//trim(outfile) )
      
    !--- write out data
    error = nf_put_var_double( ncid, id_geolon, geolon(:dim1,:dim2))
    call netcdf_err(error, 'write var geolon for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_geolat, geolat(:dim1,:dim2))
    call netcdf_err(error, 'write var geolat for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_slmsk, slm(:dim1,:dim2))
    call netcdf_err(error, 'write var slmsk for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_land_frac, land_frac(:dim1,:dim2))
    call netcdf_err(error, 'write var land_frac for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_orog_raw, oro(:dim1,:dim2))
    call netcdf_err(error, 'write var orog_raw for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_orog_filt, orf(:dim1,:dim2))
    call netcdf_err(error, 'write var orog_filt for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_stddev, hprime(:dim1,:dim2,1))
    call netcdf_err(error, 'write var stddev for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_convex, hprime(:dim1,:dim2,2))
    call netcdf_err(error, 'write var convex for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_oa1, hprime(:dim1,:dim2,3))
    call netcdf_err(error, 'write var oa1 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_oa2, hprime(:dim1,:dim2,4))
    call netcdf_err(error, 'write var oa2 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_oa3, hprime(:dim1,:dim2,5))
    call netcdf_err(error, 'write var oa3 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_oa4, hprime(:dim1,:dim2,6))
    call netcdf_err(error, 'write var oa4 for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_ol1, hprime(:dim1,:dim2,7))
    call netcdf_err(error, 'write var ol1 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_ol2, hprime(:dim1,:dim2,8))
    call netcdf_err(error, 'write var ol2 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_ol3, hprime(:dim1,:dim2,9))
    call netcdf_err(error, 'write var ol3 for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_ol4, hprime(:dim1,:dim2,10))
    call netcdf_err(error, 'write var ol4 for file='//trim(outfile) )

    error = nf_put_var_double( ncid, id_theta, hprime(:dim1,:dim2,11))
    call netcdf_err(error, 'write var theta for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_gamma, hprime(:dim1,:dim2,12))
    call netcdf_err(error, 'write var gamma for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_sigma, hprime(:dim1,:dim2,13))
    call netcdf_err(error, 'write var sigma for file='//trim(outfile) )
    error = nf_put_var_double( ncid, id_elvmax, hprime(:dim1,:dim2,14))
    call netcdf_err(error, 'write var elvmax for file='//trim(outfile) )

    error = nf_close(ncid) 
    call netcdf_err(error, 'close file='//trim(outfile) )  
      
  end subroutine

!> Check NetCDF error code and output the error message.
!!
!! @param[in] err NetCDF error code
!! @param[in] string The NetCDF error message
!! @author Jordan Alpert NOAA/EMC
  subroutine netcdf_err( err, string )
      integer, intent(in) :: err
      character(len=*), intent(in) :: string
      character(len=256) :: errmsg
      include "netcdf.inc"

      if( err.EQ.NF_NOERR )return
      errmsg = NF_STRERROR(err)
      print*, trim(string), ': ', trim(errmsg)

      return
    end subroutine netcdf_err
