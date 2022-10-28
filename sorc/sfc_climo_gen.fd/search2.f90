!> @file
!! @brief Replace undefined values on the model grid with a valid
!! value at a nearby neighbor.
!! @author George Gayno @date 2018

!> Replace undefined values on the model grid with a valid
!! value at a nearby neighbor. Undefined values are typically
!! associated with isolated islands where there is no source data.
!! Routine searches a neighborhood with a radius of 100 grid points.
!! If no valid value is found, a default value is used. This
!! routine works for one tile of a cubed sphere grid. It does
!! not consider valid values at adjacent faces. That is a future
!! upgrade.
!!
!! @param[inout] field - input: field before missing values are replaced
!!                     - output: field after missing values are replaced
!! @param[in] mask field bitmap. Field defined where mask=1
!! @param[inout] idim i dimension of tile
!! @param[inout] jdim j dimension of tile
!! @param[in] tile tile number
!! @param[in] field_name field name
!! @author George Gayno @date 2018
 subroutine search2 (field, mask, idim, jdim, tile, field_name)

 use mpi
 use esmf

 implicit none

 character(len=*)                  :: field_name

 integer, intent(in)               :: idim, jdim, tile
 integer(esmf_kind_i4), intent(in) :: mask(idim,jdim)

 real(esmf_kind_r4), intent(inout) :: field(idim,jdim,20)

 integer                           :: i, j, krad, ii, jj
 integer                           :: istart, iend
 integer                           :: jstart, jend
 integer                           :: ierr

 real                              :: default_value
 real(esmf_kind_r4), allocatable   :: field_save(:,:,:)

!-----------------------------------------------------------------------
! Set default value.
!-----------------------------------------------------------------------


 select case (field_name)
   case ('substrate_temperature') ! soil substrate_temperature
     default_value = 280.0
   case ('vegetation_greenness') ! vegetation greenness
     default_value = 0.5
   case ('maximum_snow_albedo') ! maximum snow albedo
     default_value = 0.5
   case ('leaf_area_index') ! leaf area index
     default_value = 1.0
   case ('visible_black_sky_albedo') ! visible black sky albedo
     default_value = 0.1
   case ('visible_white_sky_albedo') ! visible white sky albedo
     default_value = 0.1
   case ('near_IR_black_sky_albedo') ! near IR black sky albedo
     default_value = 0.2
   case ('near_IR_white_sky_albedo') ! near IR white sky albedo
     default_value = 0.2
   case ('facsf') ! facsf
     default_value = 0.5
   case ('facwf') ! facwf
     default_value = 0.5
   case ('slope_type') ! slope type
     default_value = float(1)
   case ('soil_type') ! soil type
     default_value = float(2)
   case ('vegetation_type') ! vegetation type
     default_value = float(3)
   case default
     print*,'- FATAL ERROR IN ROUTINE SEARCH.  UNIDENTIFIED FIELD : ', field
     call mpi_abort(mpi_comm_world, 77, ierr)
 end select

!-----------------------------------------------------------------------
! Perform search and replace.
!-----------------------------------------------------------------------

 allocate (field_save(idim,jdim,20))
 field_save = field

 J_LOOP : do j = 1, jdim
   I_LOOP : do i = 1, idim

     if (mask(i,j) == 1 .and. field_save(i,j,1) < -9999.0) then

       KRAD_LOOP : do krad = 1, 100

         istart = i - krad
         iend   = i + krad
         jstart = j - krad
         jend   = j + krad

         JJ_LOOP : do jj = jstart, jend
         II_LOOP : do ii = istart, iend

!-----------------------------------------------------------------------
!          Search only along outer square.
!-----------------------------------------------------------------------

           if ((jj == jstart) .or. (jj == jend) .or.   &
               (ii == istart) .or. (ii == iend)) then

             if (jj < 1 .or. jj > jdim) cycle JJ_LOOP
             if (ii < 1 .or. ii > idim) cycle II_LOOP

               print*,'in search ',ii,jj,mask(ii,jj),maxval(field_save(ii,jj,:))
               if (mask(ii,jj) == 1  .and. maxval(field_save(ii,jj,:)) > 0.0) then
                 field(i,j,:) = field_save(ii,jj,:)
                 write(6,100) tile,i,j,ii,jj,field(i,j,1)
                 cycle I_LOOP
               endif

           endif

         enddo II_LOOP
         enddo JJ_LOOP

       enddo KRAD_LOOP

       field(i,j,:) = 0.0
       field(i,j,nint(default_value)) = 1.0  ! Search failed.  Use default value.

       write(6,101) tile,i,j,default_value

     endif
   enddo I_LOOP
 enddo J_LOOP

 print*,'after search 59/166 ',field(59,166,:)
 print*,'after search 60/167 ',field(60,167,:)
 print*,'after search 55/168 ',field(55,168,:)
 print*,'after search 56/169 ',field(55,168,:)
 deallocate(field_save)

 100 format(1x,"- MISSING2 POINT TILE: ",i2," I/J: ",i5,i5," SET TO VALUE AT: ",i5,i5,". NEW VALUE IS: ",f8.3)
 101 format(1x,"- MISSING2 POINT TILE: ",i2," I/J: ",i5,i5," SET TO DEFAULT VALUE OF: ",f8.3)

 end subroutine search2
