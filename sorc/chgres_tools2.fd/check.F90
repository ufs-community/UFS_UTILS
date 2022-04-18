 program check

! How to use.
!
! 1) Run chgres in fractional grid mode.
!
! 2) Run chgres again using the updated orography files where 'slmsk' is
! updated to be 'floor(land fraction)'. Here chgres is run
! in non-fractional mode.
!
! 3) Compare the fields at land from (2) to the fields with at least
! some non-land from (1). They should match.

 use netcdf

 implicit none

 integer, parameter :: num_var=7
 integer, parameter :: num_var3d=3

 character(len=150) :: file_floor, file_frac, file_orog_frac

 character(len=15) :: varname(num_var), varname3d(num_var3d)

 integer :: error, id_dim, idim, jdim
 integer :: varid, varid_floor, varid_frac, var
 integer :: i, j, n, ncid_floor, ncid_orog_frac, ncid_frac

 real*8, allocatable :: slmsk(:,:), land_frac(:,:), dummy_frac(:,:)
 real*8, allocatable :: dummy_floor(:,:)
 real*8, allocatable :: dummy_frac3d(:,:,:), dummy_floor3d(:,:,:)

 data varname /'alvsf_nl', 'alvwf_nl', 'alnsf_nl', 'alnwf_nl', &
               'f10m',  'ffhh', 'ffmm'  /

 data varname3d /'stc', 'slc', 'smc'/

 file_frac="/gpfs/dell1/stmp/George.Gayno/chgres_fractional/out.sfc.tile1.nc"
 file_orog_frac="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/chgres_cube.fractional/my_grids_fract/C96/C96_oro_data.tile1.nc"
 file_floor="/gpfs/dell1/stmp/George.Gayno/chgres_floor/out.sfc.tile1.nc"

! Open the file created using the non-fractional logic, but with
! slmsk modified to be floor of the land fraction.

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

 allocate(slmsk(idim,jdim))

 error=nf90_inq_varid(ncid_floor,"slmsk",varid)
 call netcdf_err(error, 'reading slmsk id' )
 error=nf90_get_var(ncid_floor,varid,slmsk)
 call netcdf_err(error, 'reading slmsk' )

 print*,'slmsk ',maxval(slmsk),minval(slmsk)

! Open the fractional grid orography file. Read the
! land fraction.

 error=nf90_open(trim(file_orog_frac),nf90_nowrite,ncid_orog_frac)
 call netcdf_err(error, 'opening file_orog_frac' )

 allocate(land_frac(idim,jdim))

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
   if (nint(slmsk(i,j)) == 0) then
     if(floor(land_frac(i,j)) /= 0.0_8) then
       print*,'bad mask point 1', i,j,slmsk(i,j),floor(land_frac(i,j))
       stop
     endif
   endif
   if (floor(land_frac(i,j)) == 0.0_8) then
     if (nint(slmsk(i,j)) /= 0) then
       print*,'bad mask point 2', i,j,slmsk(i,j),floor(land_frac(i,j))
       stop
     endif
   endif
 enddo
 enddo

! Check data.

 allocate(dummy_frac(idim,jdim))
 allocate(dummy_floor(idim,jdim))

 error=nf90_open(trim(file_frac),nf90_nowrite,ncid_frac)
 call netcdf_err(error, 'opening file_frac' )

 do var = 1, num_var
   
   print*,'CHECK FIELD ', trim(varname(var))

   error=nf90_inq_varid(ncid_frac, varname(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac)
   call netcdf_err(error, 'reading frac field' )

   print*,'frac field ',maxval(dummy_frac),minval(dummy_frac)

   if (trim(varname(var)) == 'alvsf_nl') varname(var) = 'alvsf'
   if (trim(varname(var)) == 'alvwf_nl') varname(var) = 'alvwf'
   if (trim(varname(var)) == 'alnsf_nl') varname(var) = 'alnsf'
   if (trim(varname(var)) == 'alnwf_nl') varname(var) = 'alnwf'
   error=nf90_inq_varid(ncid_floor, varname(var), varid_floor)
   call netcdf_err(error, 'reading floor id' )
   error=nf90_get_var(ncid_floor,varid_floor,dummy_floor)
   call netcdf_err(error, 'reading floor field' )

   print*,'floor field ',maxval(dummy_floor),minval(dummy_floor)

   do j = 1, jdim
   do i = 1, idim

     if (nint(slmsk(i,j)) == 0) then
       if (dummy_floor(i,j) /= dummy_frac(i,j)) then
         print*,'bad pt ',i,j,dummy_floor(i,j),dummy_frac(i,j)
         stop
       endif
     endif

   enddo
   enddo

 enddo

 stop

 allocate(dummy_frac3d(idim,jdim,4))
 allocate(dummy_floor3d(idim,jdim,4))

 do var = 1, num_var3d
   
   print*,'CHECK FIELD ', trim(varname3d(var))

   error=nf90_inq_varid(ncid_frac, varname3d(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac3d)
   call netcdf_err(error, 'reading frac field' )

   error=nf90_inq_varid(ncid_floor, varname3d(var), varid_floor)
   call netcdf_err(error, 'reading floor id' )
   error=nf90_get_var(ncid_floor,varid_floor,dummy_floor3d)
   call netcdf_err(error, 'reading floor field' )

   do n = 1, 4

     print*,'frac field level ',n,maxval(dummy_frac3d(:,:,n)),minval(dummy_frac3d(:,:,n))
     print*,'floor field level ',n,maxval(dummy_floor3d(:,:,n)),minval(dummy_floor3d(:,:,n))

     do j = 1, jdim
     do i = 1, idim

       if (nint(slmsk(i,j)) == 1) then
         if (dummy_floor3d(i,j,n) /= dummy_frac3d(i,j,n)) then
           print*,'bad pt ',i,j,n,dummy_floor3d(i,j,n),dummy_frac3d(i,j,n)
         endif
       endif

     enddo
     enddo

   enddo

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
