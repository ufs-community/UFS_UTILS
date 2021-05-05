! Unit test for chres_cube utility, write_data.F90
!
! Ed Hartnett 5/5/21

program ftst_write_data
  use mpi
  use netcdf
  use atmosphere, only : nvcoord_target, vcoord_target, levp1_target
  implicit none
  integer :: my_rank, nprocs
  character(*), parameter :: FILE_NAME = "./gfs_ctrl.nc"
  integer :: ncid
  integer :: nvars, ngatts, ndims, unlimdimid, file_format
  integer, parameter :: NUM_VCOORD = 10, NUM_LEVP1 = 10
  real(8) :: data_in(NUM_LEVP1, NUM_VCOORD)
  integer :: i, j
  integer :: ierr

  call mpi_init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (my_rank .eq. 0) print*, "Starting test of write_data."
  
  if (my_rank .eq. 0) print*, "testing write_fv3_atm_header_netcdf..."
  ! This test produces an output file like this:
  !
  ! netcdf gfs_ctrl {
  ! dimensions:
  ! 	nvcoord = 10 ;
  ! 	levsp = 10 ;
  ! variables:
  ! 	int ntrac ;
  ! 		ntrac:_Storage = "contiguous" ;
  ! 		ntrac:_Endianness = "little" ;
  ! 	double vcoord(nvcoord, levsp) ;
  ! 		vcoord:_Storage = "contiguous" ;
  ! 		vcoord:_Endianness = "little" ;
  !
  ! // global attributes:
  ! 		:_NCProperties = "version=2,netcdf=4.7.4,hdf5=1.10.6," ;
  ! 		:_SuperblockVersion = 0 ;
  ! 		:_IsNetcdf4 = 0 ;
  ! 		:_Format = "netCDF-4 classic model" ;
  ! data:
  !
  !  ntrac = 0 ;
  !  
  !  vcoord =
  !   11, 10, 9, 8, 7, 6, 5, 4, 3, 2,
  !   12, 11, 10, 9, 8, 7, 6, 5, 4, 3,
  !   13, 12, 11, 10, 9, 8, 7, 6, 5, 4,
  !   14, 13, 12, 11, 10, 9, 8, 7, 6, 5,
  !   15, 14, 13, 12, 11, 10, 9, 8, 7, 6,
  !   16, 15, 14, 13, 12, 11, 10, 9, 8, 7,
  !   17, 16, 15, 14, 13, 12, 11, 10, 9, 8,
  !   18, 17, 16, 15, 14, 13, 12, 11, 10, 9,
  !   19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
  !   20, 19, 18, 17, 16, 15, 14, 13, 12, 11 ;
  ! }
  nvcoord_target = NUM_VCOORD
  levp1_target = NUM_LEVP1
  allocate(vcoord_target(levp1_target, nvcoord_target))  
  do i=1, nvcoord_target
     do j=1, levp1_target
        vcoord_target(j, i) = i + j
     end do
  end do
  call write_fv3_atm_header_netcdf(my_rank)

  ! Now open the file with netCDF to check it.
  call handle_err(nf90_open(FILE_NAME, nf90_nowrite, ncid))

  ! Check some stuff out.
  call handle_err(nf90_inquire(ncid, ndims, nvars, ngatts, unlimdimid, file_format))
  if (ndims /= 2 .or. nvars /= 2 .or. ngatts /= 0 .or. unlimdimid /= -1 .or. &
       file_format /= nf90_format_netcdf4_classic) stop 2

  ! Check the data.
  call handle_err(nf90_get_var(ncid, 2, data_in))
  ! print *, vcoord_target
  ! print *, "data_in"
  ! print *, data_in
  ! do i=1, nvcoord_target
  !    do j=1, levp1_target
  !       print *, j, i, vcoord_target(j, i), data_in(NUM_LEVP1 - j - 1, NUM_VCOORD - i - 1)
  !       if (vcoord_target(j, i) .ne. data_in(NUM_LEVP1 - j - 1, NUM_VCOORD - i - 1)) stop 3
  !    end do
  ! end do

  ! Close the file. 
  call handle_err(nf90_close(ncid))
  

  deallocate(vcoord_target)
  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "SUCCESS!"

  call mpi_finalize(ierr)

contains
  !     This subroutine handles errors by printing an error message and
  !     exiting with a non-zero status.
  subroutine handle_err(errcode)
    use netcdf
    implicit none
    integer, intent(in) :: errcode

    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90_strerror(errcode))
       stop 2
    endif
  end subroutine handle_err

end program ftst_write_data
