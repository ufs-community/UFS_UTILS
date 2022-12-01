!> @file
!! @brief Write model categorical data for a single tile.
!! @author George Gayno NCEP/EMC @date 2022

!> Output categorical data such as vegetation type. Include
!! percentage of each category within a model grid box and
!! the dominate category.
!!
!! @author George Gayno NCEP/EMC @date 2022
 module output_frac_cats

 implicit none

 private

 public :: output_driver

 contains

!> Driver routine to output model categorical data.
!!
!! @param[in] data_one_tile The percentage of each category within a model grid cell.
!! @param[in] dom_cat_one_tile The dominate category within a model grid cell.
!! @param[in] lat_one_tile Latitude of each model grid cell.
!! @param[in] lon_one_tile Longitude of each model grid cell.
!! @param[in] i_mdl i dimension of model grid.
!! @param[in] j_mdl j dimension of model grid.
!! @param[in] num_categories Number of categories.
!! @param[in] tile Tile number.
!! @author George Gayno @date 2022
 subroutine output_driver(data_one_tile, dom_cat_one_tile, lat_one_tile, lon_one_tile, &
                          i_mdl, j_mdl, num_categories, tile)

 use mpi
 use esmf
 use source_grid, only             : field_names
 use model_grid, only              : grid_tiles
 use program_setup, only           : halo

 implicit none

 integer, intent(in)              :: i_mdl, j_mdl, tile, num_categories

 real(esmf_kind_r4), intent(in)   :: data_one_tile(i_mdl,j_mdl,num_categories)
 real(esmf_kind_r4), intent(in)   :: dom_cat_one_tile(i_mdl,j_mdl)
 real(esmf_kind_r4), intent(in)   :: lat_one_tile(i_mdl,j_mdl)
 real(esmf_kind_r4), intent(in)   :: lon_one_tile(i_mdl,j_mdl)
 
 character(len=200)               :: out_file
 character(len=200)               :: out_file_with_halo

 integer                          :: field_idx
 integer                          :: ierr
 integer                          :: i_out, j_out
 integer                          :: i_start, i_end, j_start, j_end

 field_idx = 1

 select case (field_names(field_idx))
   case ('soil_type')
     out_file = "./soil_type." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./soil_type." // grid_tiles(tile) // ".halo.nc"
   case ('vegetation_type')
     out_file = "./vegetation_type." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./vegetation_type." // grid_tiles(tile) // ".halo.nc"
   case default
     print*,'- FATAL ERROR IN ROUTINE OUTPUT.  UNIDENTIFIED FIELD : ', field_names(field_idx)
     call mpi_abort(mpi_comm_world, 67, ierr)
 end select

!----------------------------------------------------------------------
! If user specified a halo (for running stand-alone regional grid),
! remove it.
!----------------------------------------------------------------------

 if (halo > 0) then
   print*,"- WILL WRITE WITHOUT HALO REGION OF ", halo, " ROWS/COLS."
   i_start = 1 + halo
   i_end   = i_mdl - halo
   j_start = 1 + halo
   j_end   = j_mdl - halo
   i_out   = i_end - i_start + 1
   j_out   = j_end - j_start + 1
   call writeit(out_file, i_out, j_out, num_categories, &
              lat_one_tile(i_start:i_end,j_start:j_end), &
              lon_one_tile(i_start:i_end,j_start:j_end), &
              data_one_tile(i_start:i_end,j_start:j_end,:), &
              dom_cat_one_tile(i_start:i_end,j_start:j_end) )
   print*,"- WILL WRITE FULL DOMAIN INCLUDING HALO."
   call writeit(out_file_with_halo, i_mdl, j_mdl, num_categories, &
                lat_one_tile, lon_one_tile, data_one_tile, dom_cat_one_tile)
 else
   print*,"- WILL WRITE DATA."
   call writeit(out_file, i_mdl, j_mdl, num_categories, &
                lat_one_tile, lon_one_tile, data_one_tile, dom_cat_one_tile)
 endif

 return

 end subroutine output_driver

!> Write data to a netcdf file.
!!
!! @param[in] out_file Output file name.
!! @param[in] iout i-dimension of data.
!! @param[in] jout j-dimension of data.
!! @param[in] num_categories Number of categories.
!! @param[in] latitude Latitude of data.
!! @param[in] longitude Longitude of data.
!! @param[in] data_pct Percentage of each category in each model grid cell.
!! @param[in] dominate_cat Dominate category in each model grid cell.
 subroutine writeit(out_file, iout, jout, num_categories, &
                    latitude, longitude, data_pct, dominate_cat)

 use esmf
 use netcdf
 use utils
 use source_grid, only  : day_of_rec, source, field_names, num_time_recs
 use model_grid, only   : missing

 implicit none

 character(len=*), intent(in) :: out_file

 integer, intent(in) :: iout, jout, num_categories

 real(esmf_kind_r4), intent(in)  :: latitude(iout,jout)
 real(esmf_kind_r4), intent(in)  :: longitude(iout,jout)
 real(esmf_kind_r4), intent(in)  :: data_pct(iout,jout,num_categories)
 real(esmf_kind_r4), intent(in)  :: dominate_cat(iout,jout)

 character(len=200)  :: field_names_pct
 integer             :: header_buffer_val = 16384
 integer             :: ncid, dim_x, dim_y, dim_z, dim_time
 integer             :: id_times, id_lat, id_lon, id_data_pct
 integer             :: id_data_dom_cat, id_sum
 integer             :: error

 real :: sum_all(iout,jout)

 print*,"- OPEN AND WRITE: ",trim(out_file)
 error = nf90_create(out_file, NF90_NETCDF4, ncid)
 call netcdf_err(error, 'ERROR IN NF90_CREATE' )
 error = nf90_def_dim(ncid, 'nx', iout, dim_x)
 call netcdf_err(error, 'DEFINING NX DIMENSION' )
 error = nf90_def_dim(ncid, 'ny', jout, dim_y)
 call netcdf_err(error, 'DEFINING NY DIMENSION' )
 error = nf90_def_dim(ncid, 'num_categories', num_categories, dim_z)
 call netcdf_err(error, 'DEFINING NZ DIMENSION' )
 error = nf90_def_dim(ncid, 'time', num_time_recs, dim_time)
 call netcdf_err(error, 'DEFINING TIME DIMENSION' )

 error = nf90_def_var(ncid, 'time', NF90_FLOAT, dim_time, id_times)
 call netcdf_err(error, 'DEFINING TIME VARIABLE' )
 error = nf90_put_att(ncid, id_times, "units", "days since 2015-1-1")
 call netcdf_err(error, 'DEFINING TIME ATTRIBUTE' )
 if (len_trim(source) > 0) then
   error = nf90_put_att(ncid, nf90_global, 'source', source)
   call netcdf_err(error, 'DEFINING GLOBAL SOURCE ATTRIBUTE' )
 endif

 error = nf90_def_var(ncid, 'geolat', NF90_FLOAT, (/dim_x,dim_y/), id_lat)
 call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
 error = nf90_put_att(ncid, id_lat, "long_name", "Latitude")
 call netcdf_err(error, 'DEFINING GEOLAT NAME ATTRIBUTE' )
 error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
 call netcdf_err(error, 'DEFINING GEOLAT UNIT ATTRIBUTE' )
 error = nf90_def_var(ncid, 'geolon', NF90_FLOAT, (/dim_x,dim_y/), id_lon)
 call netcdf_err(error, 'DEFINING GEOLON FIELD' )
 error = nf90_put_att(ncid, id_lon, "long_name", "Longitude")
 call netcdf_err(error, 'DEFINING GEOLON NAME ATTRIBUTE' )
 error = nf90_put_att(ncid, id_lon, "units", "degrees_east")
 call netcdf_err(error, 'DEFINING GEOLON UNIT ATTRIBUTE' )

 field_names_pct = trim(field_names(1)) // "_pct"
 error = nf90_def_var(ncid, trim(field_names_pct), NF90_FLOAT, (/dim_x,dim_y,dim_z,dim_time/), id_data_pct)
 call netcdf_err(error, 'DEFINING FIELD' )
 error = nf90_put_att(ncid, id_data_pct, "units", "percent coverage each category")
 call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
 error = nf90_put_att(ncid, id_data_pct, "missing_value", missing)
 call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
 error = nf90_put_att(ncid, id_data_pct, "coordinates", "geolon geolat")
 call netcdf_err(error, 'DEFINING COORD ATTRIBUTE' )

 error = nf90_def_var(ncid, trim(field_names(1)), NF90_FLOAT, (/dim_x,dim_y,dim_time/), id_data_dom_cat)
 call netcdf_err(error, 'DEFINING FIELD' )
 error = nf90_put_att(ncid, id_data_dom_cat, "units", "dominate category")
 call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
 error = nf90_put_att(ncid, id_data_dom_cat, "missing_value", missing)
 call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
 error = nf90_put_att(ncid, id_data_dom_cat, "coordinates", "geolon geolat")
 call netcdf_err(error, 'DEFINING COORD ATTRIBUTE' )

 error = nf90_def_var(ncid, 'sum', NF90_FLOAT, (/dim_x,dim_y,dim_time/), id_sum)
 call netcdf_err(error, 'DEFINING FIELD' )

 error = nf90_enddef(ncid, header_buffer_val,4,0,4)

 error = nf90_put_var( ncid, id_times, day_of_rec) 
 call netcdf_err(error, 'WRITING TIME FIELD' )

 error = nf90_put_var( ncid, id_lat, latitude)
 call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLAT' )

 error = nf90_put_var( ncid, id_lon, longitude)
 call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLON' )

 error = nf90_put_var( ncid, id_data_pct, data_pct)
 call netcdf_err(error, 'IN NF90_PUT_VAR' )
  
 error = nf90_put_var( ncid, id_data_dom_cat, dominate_cat)
 call netcdf_err(error, 'IN NF90_PUT_VAR' )

! Temporary output of sum of %.
 sum_all = sum(data_pct, dim=3)
 error = nf90_put_var( ncid, id_sum, sum_all)

 error = nf90_close(ncid)

 end subroutine writeit

 end module output_frac_cats
