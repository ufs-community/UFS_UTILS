 program check

! How to use.
!
! 1) Run chgres in fractional grid mode.
!
! 2) Run chgres again using the updated orography files where 'slmsk' is
! updated to be 'floor(land fraction)'. Here chgres is run
! in non-fractional mode.
!
! 3) Compare the fields at non-land from (2) to the fields with some
! non-land from (1). They should match.

 use netcdf

 implicit none

 integer, parameter :: num_var=39
 integer, parameter :: num_var3d=1

 character(len=150) :: file_floor, file_frac, file_orog_frac

 character(len=15) :: varname_frac(num_var), varname3d_frac(num_var3d)
 character(len=15) :: varname_floor(num_var), varname3d_floor(num_var3d)
 character(len=16) :: the_files(6)
 character(len=21) :: oro_files(6)

 integer :: error, id_dim, idim, jdim, tiles
 integer :: varid, varid_floor, varid_frac, var
 integer :: i, j, n, ncid_floor, ncid_orog_frac, ncid_frac

 real*8, allocatable :: slmsk(:,:), land_frac(:,:), dummy_frac(:,:)
 real*8, allocatable :: dummy_floor(:,:)
 real*8, allocatable :: dummy_frac3d(:,:,:), dummy_floor3d(:,:,:)

 data varname_frac /'alvsf_nl', 'alvwf_nl', 'alnsf_nl', 'alnwf_nl', &
               'f10m',  'ffhh', 'ffmm', 'fice', 'hice', &
               'q2m', 'sheleg_ice', 'snwdph_ice', 'srflag', &
               't2m', 'tg3_ice', 'tisfc', 'tprcp',         & 
               'tsea', 'uustar', 'zorl_ice', 'zorl', &
               'c_0',  'c_d', 'd_conv', 'dt_cool', 'ifd', &
               'qrain', 'tref', 'w_0', 'w_d', 'xs', &
               'xt',  'xtts', 'xu', 'xv', 'xz', &
               'xzts', 'z_c', 'zm'/ 

 data varname_floor /'alvsf', 'alvwf', 'alnsf', 'alnwf', &
               'f10m',  'ffhh', 'ffmm', 'fice', 'hice', &
               'q2m', 'sheleg', 'snwdph', 'srflag', &
               't2m', 'tg3', 'tisfc', 'tprcp',         &
               'tsea', 'uustar', 'zorl', 'zorl', &
               'c_0',  'c_d', 'd_conv', 'dt_cool', 'ifd', &
               'qrain', 'tref', 'w_0', 'w_d', 'xs', &
               'xt',  'xtts', 'xu', 'xv', 'xz', &
               'xzts', 'z_c', 'zm'/ 

 data varname3d_frac /'stc_ice'/

 data varname3d_floor /'stc'/

 data the_files /'out.sfc.tile1.nc', 'out.sfc.tile2.nc', 'out.sfc.tile3.nc', &
                 'out.sfc.tile4.nc', 'out.sfc.tile5.nc', 'out.sfc.tile6.nc'/
 
 data oro_files /'C96_oro_data.tile1.nc', 'C96_oro_data.tile2.nc', 'C96_oro_data.tile3.nc', &
                 'C96_oro_data.tile4.nc', 'C96_oro_data.tile5.nc', 'C96_oro_data.tile6.nc' /

 do tiles = 1, 6

 file_frac="/scratch2/NCEPDEV/stmp1/George.Gayno/chgres_fractional/" // the_files(tiles)
 file_orog_frac="/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/chgres_cube.fractional/my_grids_fract/C96/" &
                 // oro_files(tiles)
 file_floor="/scratch2/NCEPDEV/stmp1/George.Gayno/chgres_floor/" // the_files(tiles)

!file_frac="/gpfs/dell1/stmp/George.Gayno/chgres_fractional/" // the_files(tiles)
!file_orog_frac="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/chgres_cube.fractional/my_grids_fract/C96/" &
!                // oro_files(tiles)
!file_floor="/gpfs/dell1/stmp/George.Gayno/chgres_floor/" // the_files(tiles)

! Open the file created using the non-fractional logic, but with
! slmsk modified to be floor of the land fraction.

 print*,'- OPEN FILE ',trim(file_floor)
 error=nf90_open(trim(file_floor),nf90_nowrite,ncid_floor)
 call netcdf_err(error, 'opening file_floor' )

 error=nf90_inq_dimid(ncid_floor, 'xaxis_1', id_dim)
 call netcdf_err(error, 'getting xaxis' )
 error=nf90_inquire_dimension(ncid_floor,id_dim,len=idim)
 call netcdf_err(error, 'reading xaxis' )

 error=nf90_inq_dimid(ncid_floor, 'yaxis_1', id_dim)
 call netcdf_err(error, 'getting yaxis' )
 error=nf90_inquire_dimension(ncid_floor,id_dim,len=jdim)
 call netcdf_err(error, 'reading yaxis' )

 if (.not. allocated(slmsk)) allocate(slmsk(idim,jdim))

 error=nf90_inq_varid(ncid_floor,"slmsk",varid)
 call netcdf_err(error, 'reading slmsk id' )
 error=nf90_get_var(ncid_floor,varid,slmsk)
 call netcdf_err(error, 'reading slmsk' )

 print*,'slmsk ',maxval(slmsk),minval(slmsk)

! Open the fractional grid orography file. Read the
! land fraction.

 print*,'- OPEN FILE ',trim(file_orog_frac)
 error=nf90_open(trim(file_orog_frac),nf90_nowrite,ncid_orog_frac)
 call netcdf_err(error, 'opening file_orog_frac' )

 if (.not. allocated(land_frac)) allocate(land_frac(idim,jdim))

 error=nf90_inq_varid(ncid_orog_frac,"land_frac",varid_frac)
 call netcdf_err(error, 'reading land_frac id' )
 error=nf90_get_var(ncid_orog_frac,varid_frac,land_frac)
 call netcdf_err(error, 'reading land_frac' )

 print*,'land_frac ',maxval(land_frac),minval(land_frac)

! Are the land fraction from the fractional orography file
! and the 'floor' slmsk created from the non-fractional
! logic consistent with each other?

 do j = 1, jdim
 do i = 1, idim
   if (nint(slmsk(i,j)) == 0 .or. nint(slmsk(i,j)) == 2) then
     if(floor(land_frac(i,j)) /= 0.0_8) then
       print*,'bad mask point 1', i,j,slmsk(i,j),floor(land_frac(i,j))
       stop
     endif
   endif
   if (floor(land_frac(i,j)) == 0.0_8) then
     if (nint(slmsk(i,j)) == 1) then
       print*,'bad mask point 2', i,j,slmsk(i,j),floor(land_frac(i,j))
       stop
     endif
   endif
 enddo
 enddo

! Check data.

 if(.not.allocated(dummy_frac)) allocate(dummy_frac(idim,jdim))
 if(.not.allocated(dummy_floor)) allocate(dummy_floor(idim,jdim))

 print*,'- OPEN FILE ',trim(file_frac)
 error=nf90_open(trim(file_frac),nf90_nowrite,ncid_frac)
 call netcdf_err(error, 'opening file_frac' )

 print*,'process tile number ',tiles

 do var = 1, num_var
   
   print*,'CHECK FIELD ', trim(varname_frac(var))

   error=nf90_inq_varid(ncid_frac, varname_frac(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac)
   call netcdf_err(error, 'reading frac field' )

   print*,'frac field ',maxval(dummy_frac),minval(dummy_frac)

   error=nf90_inq_varid(ncid_floor, varname_floor(var), varid_floor)
   call netcdf_err(error, 'reading floor id' )
   error=nf90_get_var(ncid_floor,varid_floor,dummy_floor)
   call netcdf_err(error, 'reading floor field' )

   print*,'floor field ',maxval(dummy_floor),minval(dummy_floor)

   if (trim(varname_frac(var)) == 'tg3_ice' .or. &
       trim(varname_frac(var)) == 'tisfc' .or. &
       trim(varname_frac(var)) == 'snwdph_ice' .or. &
       trim(varname_frac(var)) == 'sheleg_ice' .or. &
       trim(varname_frac(var)) == 'zorl_ice') then ! check at ice only.
     do j = 1, jdim
     do i = 1, idim
       if (nint(slmsk(i,j)) == 2) then
         if (dummy_floor(i,j) /= dummy_frac(i,j)) then
           print*,'bad ice pt ',i,j,dummy_floor(i,j),dummy_frac(i,j)
           stop
         endif
       endif
     enddo
     enddo
   elseif (trim(varname_frac(var)) == 'zorl' .or. &
           trim(varname_frac(var)) == 'tsea'  .or. &
           trim(varname_frac(var)) == 'c_0'  .or. &
           trim(varname_frac(var)) == 'c_d'  .or. &
           trim(varname_frac(var)) == 'd_conv'  .or. &
           trim(varname_frac(var)) == 'dt_cool'  .or. &
           trim(varname_frac(var)) == 'qrain'  .or. &
           trim(varname_frac(var)) == 'w_0'  .or. &
           trim(varname_frac(var)) == 'w_d'  .or. &
           trim(varname_frac(var)) == 'xs'  .or. &
           trim(varname_frac(var)) == 'xt'  .or. &
           trim(varname_frac(var)) == 'xtts'  .or. &
           trim(varname_frac(var)) == 'xu'  .or. &
           trim(varname_frac(var)) == 'xv'  .or. &
           trim(varname_frac(var)) == 'xz'  .or. &
           trim(varname_frac(var)) == 'xzts'  .or. &
           trim(varname_frac(var)) == 'z_c'  .or. &
           trim(varname_frac(var)) == 'zm'  .or. &
           trim(varname_frac(var)) == 'tref') then ! check at open water only.
     do j = 1, jdim
     do i = 1, idim
       if (nint(slmsk(i,j)) == 0) then
         if (dummy_floor(i,j) /= dummy_frac(i,j)) then
           print*,'bad water pt ',i,j,dummy_floor(i,j),dummy_frac(i,j)
           stop
         endif
       endif
     enddo
     enddo
   else      ! Check at ice and open water.
     do j = 1, jdim
     do i = 1, idim
       if (nint(slmsk(i,j)) == 0 .or. nint(slmsk(i,j)) == 2) then
         if (dummy_floor(i,j) /= dummy_frac(i,j)) then
           print*,'bad pt ',i,j,dummy_floor(i,j),dummy_frac(i,j)
           stop
         endif
       endif
     enddo
     enddo
   endif

 enddo

 if(.not.allocated(dummy_frac3d)) allocate(dummy_frac3d(idim,jdim,2))
 if(.not.allocated(dummy_floor3d)) allocate(dummy_floor3d(idim,jdim,4))

 do var = 1, num_var3d
   
   print*,'CHECK FIELD ', trim(varname3d_frac(var))

   error=nf90_inq_varid(ncid_frac, varname3d_frac(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac3d)
   call netcdf_err(error, 'reading frac field' )

   error=nf90_inq_varid(ncid_floor, varname3d_floor(var), varid_floor)
   call netcdf_err(error, 'reading floor id' )
   error=nf90_get_var(ncid_floor,varid_floor,dummy_floor3d)
   call netcdf_err(error, 'reading floor field' )

   do n = 1, 2

     print*,'frac field level ',n,maxval(dummy_frac3d(:,:,n)),minval(dummy_frac3d(:,:,n))
     print*,'floor field level ',n,maxval(dummy_floor3d(:,:,n)),minval(dummy_floor3d(:,:,n))

     do j = 1, jdim
     do i = 1, idim

       if (nint(slmsk(i,j)) == 2) then
         if (dummy_floor3d(i,j,n) /= dummy_frac3d(i,j,n)) then
           print*,'bad 3d pt ',i,j,n,dummy_floor3d(i,j,n),dummy_frac3d(i,j,n)
           stop
         endif
       endif

     enddo
     enddo

   enddo

 enddo

 error=nf90_close(ncid_floor)
 error=nf90_close(ncid_orog_frac)
 error=nf90_close(ncid_frac)

 enddo

 print*,'DONE'

 end program check

 subroutine netcdf_err( err, string )

 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg
 integer :: iret

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 stop 6

 return
 end subroutine netcdf_err
