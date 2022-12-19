! Unit test for chres_cube utility, reading varmap files. 
!
! Ed Hartnett 5/3/21

program ftst_program_setup_varmaps
  use mpi
  use esmf
  use program_setup
  implicit none
  integer :: my_rank, nprocs
  integer, parameter :: MAX_NAME_LEN = 20
  integer, parameter :: EXPECTED_NUM_VARS = 23
  integer, parameter :: EXPECTED_NUM_TRACERS = 7
  character(len=MAX_NAME_LEN) :: expected_var_names(EXPECTED_NUM_VARS) = [character(len=20):: 'dzdt', 'sphum', 'liq_wat', &
       'o3mr', 'ice_wat', 'rainwat', 'snowwat', 'graupel', 'vtype', 'sotype', 'vfrac', 'fricv', 'sfcr', 'tprcp', &
       'ffmm', 'f10m', 'soilw', 'soill', 'soilt', 'cnwat', 'hice', 'weasd', 'snod']
  character(len=MAX_NAME_LEN) :: expected_field_names(EXPECTED_NUM_VARS) = [character(len=20):: 'dzdt', 'sphum', &
       'liq_wat', 'o3mr', 'ice_wat', 'rainwat', 'snowwat', 'graupel', 'vtype', 'stype', 'vfrac', 'uustar', 'zorl', &
       'tprcp', 'ffmm', 'f10m', 'smc', 'slc', 'stc', 'cnwat', 'icetk', 'weasd', 'snod']
  character(len=MAX_NAME_LEN) :: expected_missing_var_methods(EXPECTED_NUM_VARS) = [character(len=20):: 'set_to_fill', 'set_to_fill', &
       'set_to_fill', 'set_to_fill', 'set_to_fill', 'set_to_fill', 'set_to_fill', 'set_to_fill', 'skip', 'skip', 'skip', 'set_to_fill', 'set_to_fill', &
       'set_to_fill', 'set_to_fill', 'set_to_fill', 'stop', 'set_to_fill', 'stop', 'set_to_fill', 'set_to_fill', 'set_to_fill', 'set_to_fill']
  real(kind=esmf_kind_r4) :: expected_missing_var_values(EXPECTED_NUM_VARS) = (/ 0.0, 1E-7, 0.0, 1E-7, 0.0, 0.0, 0.0, 0.0, 0.0, &
       0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0 /)
  character(len=MAX_NAME_LEN) :: expected_tracers_input(EXPECTED_NUM_TRACERS) = [character(len=20):: 'sphum', 'liq_wat', &
       'o3mr', 'ice_wat', 'rainwat', 'snowwat', 'graupel']
  integer :: i
  integer :: ierr

  call mpi_init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (my_rank .eq. 0) print*, "Starting test of program_setup reading varmaps."
  if (my_rank .eq. 0) print*, "testing read_varmap with GFSphys_varmap.txt..."
  varmap_file = "data/GFSphys_varmap.txt"
  input_type = "grib2"
  call read_varmap()
  if (size(chgres_var_names) .ne. EXPECTED_NUM_VARS) stop 1
  if (size(field_var_names) .ne. EXPECTED_NUM_VARS) stop 1
  if (size(missing_var_methods) .ne. EXPECTED_NUM_VARS) stop 1
  if (size(missing_var_values) .ne. EXPECTED_NUM_VARS) stop 1
  if (size(read_from_input) .ne. EXPECTED_NUM_VARS) stop 1
  do i = 1, EXPECTED_NUM_VARS
     if (trim(chgres_var_names(i)) .ne. trim(expected_var_names(i))) stop 3
     if (trim(field_var_names(i)) .ne. trim(expected_field_names(i))) stop 4
     if (trim(missing_var_methods(i)) .ne. trim(expected_missing_var_methods(i))) stop 5
     print*,'in loop ',i,missing_var_values(i),expected_missing_var_values(i)
     if (missing_var_values(i) .ne. expected_missing_var_values(i)) stop 6
     if (read_from_input(i) .neqv. .true.) stop 7
  end do
  if (num_tracers .ne. EXPECTED_NUM_TRACERS) stop 10
  if (num_tracers_input .ne. EXPECTED_NUM_TRACERS) stop 11
  do i = 1, EXPECTED_NUM_TRACERS
     if (trim(tracers_input(i)) .ne. trim(expected_tracers_input(i))) stop 12
  end do
  if (my_rank .eq. 0) print*, "OK"
  
  if (my_rank .eq. 0) print*, "SUCCESS!"

 call mpi_finalize(ierr)  
end program ftst_program_setup_varmaps
