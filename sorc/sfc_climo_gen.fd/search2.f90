!> @file
!! @brief Replace undefined values on the model grid with a valid
!! value at a nearby neighbor. This routine works for fractional
!! categorical fields.
!! @author George Gayno @date 2022

!> Replace undefined values on the model grid with a valid
!! value at a nearby neighbor. Undefined values are typically
!! associated with isolated islands where there is no source data.
!! Routine searches a neighborhood with a radius of 100 grid points.
!! If no valid value is found, a default value is used. This
!! routine works for one tile of a cubed sphere grid. It does
!! not consider valid values at adjacent faces. That is a future
!! upgrade.
!!
!! @note This routine works for fractional categorical fields.
!!
!! @param[inout] field - input: field before missing values are replaced
!!                     - output: field after missing values are replaced
!! @param[in] mask field bitmap. Field defined where mask=1
!! @param[in] idim i dimension of tile
!! @param[in] jdim j dimension of tile
!! @param[in] num_categories number of veg/soil categories
!! @param[in] tile tile number
!! @param[in] field_name field name
!! @author George Gayno @date 2022
 subroutine search2 (field, mask, idim, jdim, num_categories, tile, field_name)

 use mpi
 use esmf

 implicit none

 character(len=*)                  :: field_name

 integer, intent(in)               :: idim, jdim, tile, num_categories
 integer(esmf_kind_i4), intent(in) :: mask(idim,jdim)

 real(esmf_kind_r4), intent(inout) :: field(idim,jdim,num_categories)

 integer                           :: i, j, krad, ii, jj
 integer                           :: istart, iend
 integer                           :: jstart, jend
 integer                           :: ierr
 integer                              :: default_category

 real(esmf_kind_r4), allocatable   :: field_save(:,:,:)

!-----------------------------------------------------------------------
! Set default category.
!-----------------------------------------------------------------------

 select case (field_name)
   case ('soil_type') ! soil type
     default_category = 3
   case ('vegetation_type') ! vegetation type
     default_category = 3
   case default
     print*,'- FATAL ERROR IN ROUTINE SEARCH.  UNIDENTIFIED FIELD : ', field
     call mpi_abort(mpi_comm_world, 77, ierr)
 end select

!-----------------------------------------------------------------------
! Perform search and replace.
!-----------------------------------------------------------------------

 allocate (field_save(idim,jdim,num_categories))
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
       field(i,j,default_category) = 1.0  ! Search failed.  Use 100% of default category.

       write(6,101) tile,i,j,default_category

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
