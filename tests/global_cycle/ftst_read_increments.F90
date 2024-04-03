! Unit test for global_cycle routine "read_data".
!
! Reads a sample soil increment file from file soil_xainc.$NNN
! NOTE: $NNN corresponds to (mpi rank + 1)
! Total number of processes= 6 (corresponds to 6 tiles)
! Each process is assigned one increment file on one tile to read
! If any portion of the variable does not match expected values,
! the test fails.
!
! Author: Yuan Xue, 07/17/2023

module chdir_mod
  implicit none
  interface
    integer function c_chdir(path) bind(C,name="chdir")
      use iso_c_binding
      character(kind=c_char) :: path(*)
    end function
  end interface
contains
  subroutine chdir(path, err)
    use iso_c_binding
    character(*) :: path
    integer, optional, intent(out) :: err
    integer :: loc_err

    loc_err = c_chdir(path//c_null_char)
    if (present(err)) err = loc_err
  end subroutine
end module chdir_mod

 program read_increments

 use read_write_data, only  :  read_data
 use mpi
 use chdir_mod

 implicit none

 integer:: ierr,my_rank
 integer:: lsoil, l
 integer:: lensfc, ij_input

 integer, parameter:: NUM_VALUES=4
 real, parameter :: EPSILON=0.0001

 real, allocatable:: SLCINC(:,:),STCINC(:,:)
 real:: stc_inc_expected_values_tile1(NUM_VALUES)
 real:: stc_inc_expected_values_tile2(NUM_VALUES)
 real:: stc_inc_expected_values_tile3(NUM_VALUES)
 real:: stc_inc_expected_values_tile4(NUM_VALUES)
 real:: stc_inc_expected_values_tile5(NUM_VALUES)
 real:: stc_inc_expected_values_tile6(NUM_VALUES)
 real:: slc_inc_expected_values_tile1(NUM_VALUES)
 real:: slc_inc_expected_values_tile2(NUM_VALUES)
 real:: slc_inc_expected_values_tile3(NUM_VALUES)
 real:: slc_inc_expected_values_tile4(NUM_VALUES)
 real:: slc_inc_expected_values_tile5(NUM_VALUES)
 real:: slc_inc_expected_values_tile6(NUM_VALUES)

 !expected values were extracted from MATLAB, which directly reads in xainc file
 !each tile is examined separately here
 data stc_inc_expected_values_tile1 / 3.1428, 2.9983, 2.9786, 2.9634 /
 data stc_inc_expected_values_tile2 / 2.9330, 2.9121, 2.9103, 2.9069 /
 data stc_inc_expected_values_tile3 / 2.7236, 2.7308, 2.7315, 2.7295 /
 data stc_inc_expected_values_tile4 / 3.0229, 3.0229, 3.0229, 3.0229 /
 data stc_inc_expected_values_tile5 / 2.8595, 2.8825, 2.8878, 2.8948 /
 data stc_inc_expected_values_tile6 / 2.7238, 2.7238, 2.7238, 2.7238 /
 data slc_inc_expected_values_tile1 / 0.0007, 0.0018, 0.0018, 0.0018 /
 data slc_inc_expected_values_tile2 / 0.0034, 0.0031, 0.0029, 0.0029 /
 data slc_inc_expected_values_tile3 / 0.0003, 0.0005, 0.0011, 0.0008 /
 data slc_inc_expected_values_tile4 / 0.01, 0.01, 0.01, 0.01 / 
 data slc_inc_expected_values_tile5 / 0.0019, 0.0019, 0.0020, 0.0024 /
 data slc_inc_expected_values_tile6 / 0.01, 0.01, 0.01, 0.01 /

 call mpi_init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

 lsoil=4
 lensfc=36864
 ij_input = 28541

 allocate(SLCINC(lensfc,lsoil))
 allocate(STCINC(lensfc,lsoil))

 if (my_rank .eq. 0) print*,"Starting test of global_cycle routine read_data."
 if (my_rank .eq. 0) print*,"Call routine read_data"

 call chdir("./data")
 call read_data(lsoil,lensfc,.false.,.false.,.true.,.true.,STCINC=STCINC,SLCINC=SLCINC)

 if (my_rank .eq. 0) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile1(l))&
                                                    > EPSILON) stop 20
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile1(l))&
                                                    > EPSILON) stop 40
   enddo 
   print*, "tile#", my_rank+1, "reads OK"
 endif

 if (my_rank .eq. 1) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile2(l))&
                                                    > EPSILON) stop 21
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile2(l))&
                                                    > EPSILON) stop 41
   enddo
   print*, "tile#", my_rank+1, "reads OK"
 endif

  if (my_rank .eq. 2) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile3(l))&
                                                    > EPSILON) stop 22
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile3(l))&
                                                    > EPSILON) stop 42
   enddo
   print*, "tile#", my_rank+1, "reads OK"
 endif

  if (my_rank .eq. 3) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile4(l))&
                                                    > EPSILON) stop 23
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile4(l))&
                                                    > EPSILON) stop 43
   enddo
   print*, "tile#", my_rank+1, "reads OK"
 endif

  if (my_rank .eq. 4) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile5(l))&
                                                    > EPSILON) stop 24
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile5(l))&
                                                    > EPSILON) stop 44
   enddo
   print*, "tile#", my_rank+1, "reads OK"
 endif

  if (my_rank .eq. 5) then
   do l = 1,4
     if (abs(STCINC(ij_input,l) - stc_inc_expected_values_tile6(l))&
                                                    > EPSILON) stop 25
     if (abs(SLCINC(ij_input,l) - slc_inc_expected_values_tile6(l))&
                                                    > EPSILON) stop 45
   enddo
   print*, "tile#", my_rank+1, "reads OK"
 endif

 call MPI_Barrier(MPI_COMM_WORLD, ierr)

 if (my_rank .eq. 0) print*, "ALL is OK"
 if (my_rank .eq. 0) print*, "SUCCESS!"

 deallocate(SLCINC,STCINC)
 call mpi_finalize(ierr)

 end program read_increments

