!> @file
!! @brief Output model data for a single tile and a single record.
!! @author George Gayno @date 2018

!> Output model data for a single tile and a single
!! record in netcdf format.
!!
!! @param[in] data_one_tile Data to be output (single tile).
!! @param[in] lat_one_tile Latitude of tile.
!! @param[in] lon_one_tile Longitude of tile.
!! @param[in] field_idx Index of field within field name array.
!! @param[in] i_mdl i dimensions of tile.
!! @param[in] j_mdl j dimensions of tile.
!! @param[in] record Record number to be output.
!! @param[in] tile Tile number.
!! @param[in] time Time period to be output.
!! @author George Gayno @date 2018
 subroutine output(data_one_tile, lat_one_tile, lon_one_tile, i_mdl, j_mdl, &
                   tile, record, time, field_idx)

 use mpi
 use esmf
 use netcdf
 use utils
 use source_grid, only             : field_names, source, num_fields, &
                                     num_time_recs, num_records, day_of_rec
 use model_grid, only              : missing, grid_tiles
 use program_setup, only           : halo

 implicit none

 integer, intent(in)              :: i_mdl, j_mdl, tile
 integer, intent(in)              :: record, time, field_idx

 real(esmf_kind_r4), intent(in)   :: data_one_tile(i_mdl,j_mdl)
 real(esmf_kind_r4), intent(in)   :: lat_one_tile(i_mdl,j_mdl)
 real(esmf_kind_r4), intent(in)   :: lon_one_tile(i_mdl,j_mdl)

 character(len=200)               :: out_file
 character(len=200)               :: out_file_with_halo

 integer                          :: initialsiz, fsize, error, j
 integer                          :: dim_x, dim_y, id_data
 integer                          :: dim_time, id_times, ierr
 integer                          :: header_buffer_val = 16384
 integer                          :: i_out, j_out, id_lat, id_lon
 integer                          :: i_start, i_end, j_start, j_end
 integer, save                    :: ncid(6), ncid_with_halo

 select case (field_names(field_idx))
   case ('substrate_temperature')
     out_file = "./substrate_temperature." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./substrate_temperature." // grid_tiles(tile) // ".halo.nc"
   case ('vegetation_greenness')
     out_file = "./vegetation_greenness." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./vegetation_greenness." // grid_tiles(tile) // ".halo.nc"
   case ('maximum_snow_albedo')
     out_file = "./maximum_snow_albedo." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./maximum_snow_albedo." // grid_tiles(tile) // ".halo.nc"
   case ('leaf_area_index')
     out_file = "./leaf_area_index." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./leaf_area_index." // grid_tiles(tile) // ".halo.nc"
   case ('visible_black_sky_albedo', 'visible_white_sky_albedo', 'near_IR_black_sky_albedo', 'near_IR_white_sky_albedo')
     out_file = "./snowfree_albedo." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./snowfree_albedo." // grid_tiles(tile) // ".halo.nc"
   case ('facsf')
     out_file = "./facsf." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./facsf." // grid_tiles(tile) // ".halo.nc"
   case ('slope_type')
     out_file = "./slope_type." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./slope_type." // grid_tiles(tile) // ".halo.nc"
   case ('soil_type')
     out_file = "./soil_type." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./soil_type." // grid_tiles(tile) // ".halo.nc"
   case ('soil_color')
     out_file = "./soil_color." // grid_tiles(tile) // ".nc"
     out_file_with_halo = "./soil_color." // grid_tiles(tile) // ".halo.nc"
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
   print*,"- WILL REMOVE HALO REGION OF ", halo, " ROWS/COLS."
   i_start = 1 + halo
   i_end   = i_mdl - halo
   j_start = 1 + halo
   j_end   = j_mdl - halo
 else
   i_start = 1
   i_end   = i_mdl
   j_start = 1
   j_end   = j_mdl
 endif

 i_out = i_end - i_start + 1
 j_out = j_end - j_start + 1

 if (record == 1) then

   initialsiz = 0
   fsize      = 65536
   error = nf90_create(out_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                       ncid(tile), initialsize=initialsiz, chunksize=fsize)
   call netcdf_err(error, 'ERROR IN NF90_CREATE' )
   error = nf90_def_dim(ncid(tile), 'nx', i_out, dim_x)
   call netcdf_err(error, 'DEFINING NX DIMENSION' )
   error = nf90_def_dim(ncid(tile), 'ny', j_out, dim_y)
   call netcdf_err(error, 'DEFINING NY DIMENSION' )
   error = nf90_def_dim(ncid(tile), 'time', num_time_recs, dim_time)
   call netcdf_err(error, 'DEFINING TIME DIMENSION' )
   error = nf90_def_var(ncid(tile), 'time', NF90_FLOAT, dim_time, id_times)
   call netcdf_err(error, 'DEFINING TIME VARIABLE' )
   error = nf90_put_att(ncid(tile), id_times, "units", "days since 2015-1-1")
   call netcdf_err(error, 'DEFINING TIME ATTRIBUTE' )
   if (len_trim(source) > 0) then
     error = nf90_put_att(ncid(tile), nf90_global, 'source', source)
     call netcdf_err(error, 'DEFINING GLOBAL SOURCE ATTRIBUTE' )
   endif

   error = nf90_def_var(ncid(tile), 'geolat', NF90_FLOAT, (/dim_x,dim_y/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid(tile), id_lat, "long_name", "Latitude")
   call netcdf_err(error, 'DEFINING GEOLAT NAME ATTRIBUTE' )
   error = nf90_put_att(ncid(tile), id_lat, "units", "degrees_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNIT ATTRIBUTE' )
   error = nf90_def_var(ncid(tile), 'geolon', NF90_FLOAT, (/dim_x,dim_y/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid(tile), id_lon, "long_name", "Longitude")
   call netcdf_err(error, 'DEFINING GEOLON NAME ATTRIBUTE' )
   error = nf90_put_att(ncid(tile), id_lon, "units", "degrees_east")
   call netcdf_err(error, 'DEFINING GEOLON UNIT ATTRIBUTE' )

   do j = 1, num_fields
     error = nf90_def_var(ncid(tile), trim(field_names(j)), NF90_FLOAT, (/dim_x,dim_y,dim_time/), id_data)
     call netcdf_err(error, 'DEFINING FIELD' )
     error = nf90_put_att(ncid(tile), id_data, "missing_value", missing)
     call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
     error = nf90_put_att(ncid(tile), id_data, "coordinates", "geolon geolat")
     call netcdf_err(error, 'DEFINING COORD ATTRIBUTE' )
   enddo

   error = nf90_enddef(ncid(tile), header_buffer_val,4,0,4)
   call netcdf_err(error, 'IN NF90_ENDDEF' )

   error = nf90_put_var( ncid(tile), id_times, day_of_rec) 
   call netcdf_err(error, 'WRITING TIME FIELD' )

   error = nf90_put_var( ncid(tile), id_lat, lat_one_tile(i_start:i_end,j_start:j_end),  &
                         start=(/1,1/), count=(/i_out,j_out/))
   call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLAT' )

   error = nf90_put_var( ncid(tile), id_lon, lon_one_tile(i_start:i_end,j_start:j_end),  &
                         start=(/1,1/), count=(/i_out,j_out/))
   call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLON' )

 endif

 print*,'- WRITE DATA FOR RECORD: ',record
 error = nf90_inq_varid( ncid(tile), field_names(field_idx), id_data)
 call netcdf_err(error, 'IN NF90_INQ_VARID' )
 error = nf90_put_var( ncid(tile), id_data, data_one_tile(i_start:i_end,j_start:j_end),  &
                            start=(/1,1,time/), count=(/i_out,j_out,1/))
 call netcdf_err(error, 'IN NF90_PUT_VAR' )
  
 if (record == num_records) then
   error = nf90_close(ncid(tile))
 endif

!----------------------------------------------------------------------
! For regional nests, also output files including the halo
!----------------------------------------------------------------------

 if (halo == 0) return

 print*,"- WRITE OUT FILES THAT INCLUDE HALO REGION."

 if (record == 1) then

   initialsiz = 0
   fsize      = 65536
   error = nf90_create(out_file_with_halo, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                       ncid_with_halo, initialsize=initialsiz, chunksize=fsize)
   call netcdf_err(error, 'IN NF90_CREATE' )
   error = nf90_def_dim(ncid_with_halo, 'nx', i_mdl, dim_x)
   call netcdf_err(error, 'DEFINING NX DIMENSION' )
   error = nf90_def_dim(ncid_with_halo, 'ny', j_mdl, dim_y)
   call netcdf_err(error, 'DEFINING NY DIMENSION' )
   error = nf90_def_dim(ncid_with_halo, 'time', num_time_recs, dim_time)
   call netcdf_err(error, 'DEFINING TIME DIMENSION' )
   error = nf90_def_var(ncid_with_halo, 'time', NF90_FLOAT, dim_time, id_times)
   call netcdf_err(error, 'DEFINING TIME VARIABLE' )
   error = nf90_put_att(ncid_with_halo, id_times, "units", "days since 2015-1-1")
   call netcdf_err(error, 'DEFINING TIME ATTRIBUTE' )
   if (len_trim(source) > 0) then
     error = nf90_put_att(ncid_with_halo, nf90_global, 'source', source)
     call netcdf_err(error, 'DEFINING GLOBAL SOURCE ATTRIBUTE' )
   endif

   error = nf90_def_var(ncid_with_halo, 'geolat', NF90_FLOAT, (/dim_x,dim_y/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid_with_halo, id_lat, "long_name", "Latitude")
   call netcdf_err(error, 'DEFINING GEOLAT NAME ATTRIBUTE' )
   error = nf90_put_att(ncid_with_halo, id_lat, "units", "degrees_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNIT ATTRIBUTE' )
   error = nf90_def_var(ncid_with_halo, 'geolon', NF90_FLOAT, (/dim_x,dim_y/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid_with_halo, id_lon, "long_name", "Longitude")
   call netcdf_err(error, 'DEFINING GEOLON NAME ATTRIBUTE' )
   error = nf90_put_att(ncid_with_halo, id_lon, "units", "degrees_east")
   call netcdf_err(error, 'DEFINING GEOLON UNIT ATTRIBUTE' )

   do j = 1, num_fields
     error = nf90_def_var(ncid_with_halo, field_names(j), NF90_FLOAT, (/dim_x,dim_y,dim_time/), id_data)
     call netcdf_err(error, 'DEFINING FIELD VARIABLE' )
     error = nf90_put_att(ncid_with_halo, id_data, "missing_value", missing)
     call netcdf_err(error, 'DEFINING FIELD ATTRIBUTE' )
     error = nf90_put_att(ncid_with_halo, id_data, "coordinates", "geolon geolat")
     call netcdf_err(error, 'DEFINING COORD ATTRIBUTE' )
   enddo

   error = nf90_enddef(ncid_with_halo, header_buffer_val,4,0,4)
   call netcdf_err(error, 'WRITING HEADER ENDDEF' )

   error = nf90_put_var(ncid_with_halo, id_times, day_of_rec) 
   call netcdf_err(error, 'WRITING TIME VARIABLE' )

   error = nf90_put_var( ncid_with_halo, id_lat, lat_one_tile,  &
                         start=(/1,1/), count=(/i_mdl,j_mdl/))
   call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLAT' )

   error = nf90_put_var( ncid_with_halo, id_lon, lon_one_tile,  &
                         start=(/1,1/), count=(/i_mdl,j_mdl/))
   call netcdf_err(error, 'IN NF90_PUT_VAR FOR GEOLON' )

 endif

 print*,'- WRITE DATA FOR RECORD: ',record
 error = nf90_inq_varid(ncid_with_halo, field_names(field_idx), id_data)
 call netcdf_err(error, 'IN NF90_INQ_VARID' )

 error = nf90_put_var(ncid_with_halo, id_data, data_one_tile,  &
                      start=(/1,1,time/), count=(/i_mdl,j_mdl,1/))
 call netcdf_err(error, 'IN NF90_PUT_VAR' )
  
 if (record == num_records) then
   error = nf90_close(ncid_with_halo)
 endif

 return

 end subroutine output
