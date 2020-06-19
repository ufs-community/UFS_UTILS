module grib2_util

!--------------------------------------------------------------------------
! Module: grib2_util
!
! Abstract: Utilities for use when reading grib2 data.
!
!--------------------------------------------------------------------------

use esmf

use model_grid, only      : i_input, j_input

implicit none

contains 

 subroutine rh2spfh(rh_sphum,p,t)
    
  implicit none
  real,parameter      :: alpha=-9.477E-4 , & !K^-1,
                         Tnot=273.15, &  !K
                         Lnot=2.5008E6, & !JKg^-1
                         Rv=461.51, & !JKg^-1K^-1
                         esnot=611.21 !Pa
  
  real(esmf_kind_r4), intent(inout), dimension(i_input,j_input) ::rh_sphum
  real(esmf_kind_r8), intent(in)                  :: p, t(i_input,j_input)

  real, dimension(i_input,j_input)  :: es, e, rh

  print*,"- CONVERT RH TO SPFH AT LEVEL ", p

  rh = rh_sphum
  !print *, 'T = ', T, ' RH = ', RH, ' P = ', P
  es = esnot * exp( Lnot/Rv * ((t-Tnot)/(t*tnot) + alpha * LOG(t/Tnot) - alpha * (t-Tnot)/ t))
  !print *, 'es = ', es
  e = rh * es / 100.0
  !print *, 'e = ', e
  rh_sphum = 0.622 * e / p
  !print *, 'q = ', sphum
  
  !if (P .eq. 100000.0) THEN
  ! print *, 'T = ', T, ' RH = ', RH, ' P = ', P, ' es = ', es, ' e = ', e, ' q = ', sphum
  !end if

end subroutine RH2SPFH

subroutine convert_omega(omega,p,t,q,clb,cub)

  implicit none
  real(esmf_kind_r8), pointer     :: omega(:,:,:), p(:,:,:), t(:,:,:), q(:,:,:),omtmp,ptmp
  
  integer                         :: clb(3), cub(3), i ,j, k
  
  real, parameter                 :: Rd = 287.15_esmf_kind_r8, &  !JKg^-1K^-1
                                     Rv=461.51_esmf_kind_r8, & !JKg^-1K^-1
                                     g = 9.81_esmf_kind_r8 ! ms^-2
                                     
  real(esmf_kind_r8)              :: tv, w
  
  do k = clb(3),cub(3)
    do j = clb(2),cub(2)
      do i = clb(1),cub(1)
        tv = t(i,j,k)*(1+Rd/Rv*q(i,j,k))
        omtmp=>omega(i,j,k)
        ptmp=>p(i,j,k)

        w = -1 * omtmp * Rd * tv / (ptmp * g)
        omega(i,j,k)=w
      enddo
    enddo
  enddo

end subroutine convert_omega

function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

subroutine write_temp(filename,fieldname,localpet,field_input,field_target)

 use esmf
 use netcdf


 use model_grid, only              : num_tiles_target_grid, &
                                     i_target, j_target, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     i_input, j_input, &
                                     longitude_input_grid, &
                                     latitude_input_grid

 implicit none

 integer, intent(in)              :: localpet
 
 character(len=128), intent(in) :: filename, fieldname
 
 real(esmf_kind_r8), intent(in), optional   :: field_input(i_input,j_input)
 real(esmf_kind_r8), intent(in), optional   :: field_target(i_target,j_target)


 integer                          :: error, ncid, tile
 integer                          :: fsize=65536, initial = 0
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_lon, dim_lat
 integer                          :: id_lon, id_lat, id_ps
 integer                          :: i_target_out, j_target_out
 
 logical                          :: istarget = .true.


 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)


 if (localpet < num_tiles_target_grid) then
   if (present(field_target)) then
     allocate(data_one_tile(i_target,j_target))
     i_target_out = i_target
     j_target_out = j_target
   elseif (present(field_input)) then
     allocate(data_one_tile(i_input,j_input))
     i_target_out = i_input
     j_target_out = j_input
     istarget=.false.
   endif
 else
   allocate(data_one_tile(0,0))
 endif


 HEADER : if (localpet < num_tiles_target_grid) then

   tile = localpet + 1

!--- open the file
   error = nf90_create(filename, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                       ncid, initialsize=initial, chunksize=fsize)
   call netcdf_err(error, 'CREATING FILE='//trim(filename) )

!--- define dimension
   error = nf90_def_dim(ncid, 'lon', i_target_out, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'lat', j_target_out, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )



!--- define field
   error = nf90_def_var(ncid, 'lon', NF90_FLOAT, (/dim_lon/), id_lon)
   call netcdf_err(error, 'DEFINING LON FIELD' )
   error = nf90_put_att(ncid, id_lon, "cartesian_axis", "X")
   call netcdf_err(error, 'WRITING LON FIELD' )
   error = nf90_def_var(ncid, 'lat', NF90_FLOAT, (/dim_lat/), id_lat)
   call netcdf_err(error, 'DEFINING LAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "cartesian_axis", "Y")
   call netcdf_err(error, 'WRITING LAT FIELD' )
   error = nf90_def_var(ncid, trim(fieldname), NF90_FLOAT, (/dim_lon,dim_lat/), id_ps)
   call netcdf_err(error, 'DEFINING FIELD' )


   error = nf90_enddef(ncid, header_buffer_val,4,0,4)
   call netcdf_err(error, 'DEFINING HEADER' )

 endif HEADER

!  longitude

 tile = 1

   if (istarget) then
     print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE FOR TILE: ", tile
     call ESMF_FieldGather(longitude_target_grid, data_one_tile, rootPet=tile-1, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)
   else
     print*,"- CALL FieldGather FOR INPUT GRID LONGITUDE FOR TILE: ", tile
     call ESMF_FieldGather(longitude_input_grid, data_one_tile, rootPet=tile-1, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)
   endif

 if (localpet < num_tiles_target_grid) then
   error = nf90_put_var( ncid, id_lon, data_one_tile(:,1))
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

!  latitude

   if (istarget) then
     print*,"- CALL FieldGather FOR TARGET GRID LATITUDE FOR TILE: ", tile
     call ESMF_FieldGather(latitude_target_grid, data_one_tile, rootPet=tile-1, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)
   else
     print*,"- CALL FieldGather FOR INPUT GRID LATITUDE FOR TILE: ", tile
     call ESMF_FieldGather(latitude_input_grid, data_one_tile, rootPet=tile-1, tile=tile, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGather", error)
   endif

 if (localpet < num_tiles_target_grid) then
   error = nf90_put_var( ncid, id_lat, data_one_tile(1,:))
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif

!  surface pressure

 print*, "- WRITING FIELD TO FILE "
 if (localpet < num_tiles_target_grid) then
   if (istarget) then
     error = nf90_put_var( ncid, id_ps, field_target)
   else
     error = nf90_put_var( ncid, id_ps, field_input)
   endif
   call netcdf_err(error, 'WRITING FIELD RECORD' )
 endif

 deallocate(data_one_tile)
 
  if (localpet < num_tiles_target_grid) error = nf90_close(ncid)

 end subroutine write_temp

 end module grib2_util
