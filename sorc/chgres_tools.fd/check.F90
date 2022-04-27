 program check

! How to use.
!
! 1) Run chgres in fractional grid mode.
!
! 2) Run chgres again using the updated orography files where 'slmsk' is
! updated to be 'ceiling(land fraction)'. Here chgres is run
! in non-fractional mode.
!
! 3) Compare the fields at land from (2) to the fields with at least
! some land from (1). They should match.

 use netcdf

 implicit none

 integer, parameter :: num_var=27
 integer, parameter :: num_var3d=3

 character(len=150) :: file_ceiling, file_frac, file_orog_frac

 character(len=15) :: varname_frac(num_var), varname3d(num_var3d)
 character(len=15) :: varname_ceil(num_var)
 character(len=21) :: oro_files(6)
 character(len=16) :: the_files(6)

 integer :: error, id_dim, idim, jdim, tiles
 integer :: varid, varid_ceil, varid_frac, var
 integer :: i, j, n, ncid_ceil, ncid_orog_frac, ncid_frac

 real*8, allocatable :: slmsk(:,:), land_frac(:,:), dummy_frac(:,:)
 real*8, allocatable :: dummy_ceil(:,:)
 real*8, allocatable :: dummy_frac3d(:,:,:), dummy_ceil3d(:,:,:)

 data varname_frac /'alvsf', 'alvwf', 'alnsf', 'alnwf', 'canopy', &
               'f10m' , 'facsf', 'facwf', 'ffhh' , 'ffmm',   &
               'q2m'  , 'shdmax','shdmin', 'sheleg', 'slope', &
               'snoalb', 'snwdph', 'srflag','stype', 't2m', &
               'tg3',  'tprcp', 'tsfcl',   'uustar', 'vfrac', &
               'vtype', 'zorll'  /

 data varname_ceil /'alvsf', 'alvwf', 'alnsf', 'alnwf', 'canopy', &
               'f10m' , 'facsf', 'facwf', 'ffhh' , 'ffmm',   &
               'q2m'  , 'shdmax','shdmin', 'sheleg', 'slope', &
               'snoalb', 'snwdph', 'srflag','stype', 't2m', &
               'tg3',  'tprcp', 'tsea',   'uustar', 'vfrac', &
               'vtype', 'zorl'  /

 data the_files /'out.sfc.tile1.nc', 'out.sfc.tile2.nc', 'out.sfc.tile3.nc', &
                 'out.sfc.tile4.nc', 'out.sfc.tile5.nc', 'out.sfc.tile6.nc'/

 data oro_files /'C96_oro_data.tile1.nc', 'C96_oro_data.tile2.nc', 'C96_oro_data.tile3.nc', &
                 'C96_oro_data.tile4.nc', 'C96_oro_data.tile5.nc', 'C96_oro_data.tile6.nc' /

 data varname3d /'stc', 'slc', 'smc'/

 TILE : do tiles = 1, 6

 file_frac="/gpfs/dell1/stmp/George.Gayno/chgres_fractional/" // the_files(tiles)
 file_orog_frac="/gpfs/dell2/emc/modeling/noscrub/George.Gayno/ufs_utils.git/chgres_cube.fractional/my_grids_fract/C96/" &
                 // oro_files(tiles)
 file_ceiling="/gpfs/dell1/stmp/George.Gayno/chgres_ceiling/" // the_files(tiles)

! Open the file created using the non-fractional logic, but with
! slmsk modified to be ceiling of the land fraction.

 print*,"- OPEN FILE ", trim(file_ceiling)
 error=nf90_open(trim(file_ceiling),nf90_nowrite,ncid_ceil)
 call netcdf_err(error, 'opening file_ceiling' )

 error=nf90_inq_dimid(ncid_ceil, 'xaxis_1', id_dim)
 call netcdf_err(error, 'getting xaxis' )
 error=nf90_inquire_dimension(ncid_ceil,id_dim,len=idim)
 call netcdf_err(error, 'reading xaxis' )

 error=nf90_inq_dimid(ncid_ceil, 'yaxis_1', id_dim)
 call netcdf_err(error, 'getting yaxis' )
 error=nf90_inquire_dimension(ncid_ceil,id_dim,len=jdim)
 call netcdf_err(error, 'reading yaxis' )

 if(.not.allocated(slmsk)) allocate(slmsk(idim,jdim))

 error=nf90_inq_varid(ncid_ceil,"slmsk",varid)
 call netcdf_err(error, 'reading slmsk id' )
 error=nf90_get_var(ncid_ceil,varid,slmsk)
 call netcdf_err(error, 'reading slmsk' )

 print*,'slmsk ',maxval(slmsk),minval(slmsk)

! Open the fractional grid orography file. Read the
! land fraction.

 print*,"- OPEN FILE ", trim(file_orog_frac)
 error=nf90_open(trim(file_orog_frac),nf90_nowrite,ncid_orog_frac)
 call netcdf_err(error, 'opening file_orog_frac' )

 if(.not.allocated(land_frac)) allocate(land_frac(idim,jdim))

 error=nf90_inq_varid(ncid_orog_frac,"land_frac",varid_frac)
 call netcdf_err(error, 'reading land_frac id' )
 error=nf90_get_var(ncid_orog_frac,varid_frac,land_frac)
 call netcdf_err(error, 'reading land_frac' )

 print*,'land_frac ',maxval(land_frac),minval(land_frac)

! Are the land fraction from the fractional orography file
! and the 'ceiling' slmsk created from the non-fractional
! logic consistent with each other?

 do j = 1, jdim
 do i = 1, idim
   if (nint(slmsk(i,j)) == 1) then
     if(ceiling(land_frac(i,j)) /= 1.0_8) then
       print*,'bad mask point 1', i,j,slmsk(i,j),ceiling(land_frac(i,j))
       stop
     endif
   endif
   if (ceiling(land_frac(i,j)) == 1.0_8) then
     if (nint(slmsk(i,j)) /= 1) then
       print*,'bad mask point 2', i,j,slmsk(i,j),ceiling(land_frac(i,j))
       stop
     endif
   endif
 enddo
 enddo

! Check data.

 if(.not.allocated(dummy_frac)) allocate(dummy_frac(idim,jdim))
 if(.not.allocated(dummy_ceil)) allocate(dummy_ceil(idim,jdim))

 print*,"- OPEN FILE ", trim(file_frac)
 error=nf90_open(trim(file_frac),nf90_nowrite,ncid_frac)
 call netcdf_err(error, 'opening file_frac' )

 do var = 1, num_var
   
   print*,'CHECK FIELD ', trim(varname_frac(var))

   error=nf90_inq_varid(ncid_frac, varname_frac(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac)
   call netcdf_err(error, 'reading frac field' )

   print*,'frac field ',maxval(dummy_frac),minval(dummy_frac)

   print*,'read variable ',varname_ceil(var)
   error=nf90_inq_varid(ncid_ceil, varname_ceil(var), varid_ceil)
   call netcdf_err(error, 'reading ceil id' )
   error=nf90_get_var(ncid_ceil,varid_ceil,dummy_ceil)
   call netcdf_err(error, 'reading ceil field' )

   print*,'ceil field ',maxval(dummy_ceil),minval(dummy_ceil)

   do j = 1, jdim
   do i = 1, idim

     if (nint(slmsk(i,j)) == 1) then
       if (dummy_ceil(i,j) /= dummy_frac(i,j)) then
         print*,'bad pt ',i,j,dummy_ceil(i,j),dummy_frac(i,j)
         stop
       endif
     endif

   enddo
   enddo

 enddo

 if(.not.allocated(dummy_frac3d)) allocate(dummy_frac3d(idim,jdim,4))
 if(.not.allocated(dummy_ceil3d)) allocate(dummy_ceil3d(idim,jdim,4))

 do var = 1, num_var3d
   
   print*,'CHECK FIELD ', trim(varname3d(var))

   error=nf90_inq_varid(ncid_frac, varname3d(var), varid_frac)
   call netcdf_err(error, 'reading frac id' )
   error=nf90_get_var(ncid_frac,varid_frac,dummy_frac3d)
   call netcdf_err(error, 'reading frac field' )

   error=nf90_inq_varid(ncid_ceil, varname3d(var), varid_ceil)
   call netcdf_err(error, 'reading ceil id' )
   error=nf90_get_var(ncid_ceil,varid_ceil,dummy_ceil3d)
   call netcdf_err(error, 'reading ceil field' )

   do n = 1, 4

     print*,'frac field level ',n,maxval(dummy_frac3d(:,:,n)),minval(dummy_frac3d(:,:,n))
     print*,'ceil field level ',n,maxval(dummy_ceil3d(:,:,n)),minval(dummy_ceil3d(:,:,n))

     do j = 1, jdim
     do i = 1, idim

       if (nint(slmsk(i,j)) == 1) then
         if (dummy_ceil3d(i,j,n) /= dummy_frac3d(i,j,n)) then
           print*,'bad 3d pt ',i,j,n,dummy_ceil3d(i,j,n),dummy_frac3d(i,j,n)
           stop
         endif
       endif

     enddo
     enddo

   enddo

 enddo

 enddo TILE

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
