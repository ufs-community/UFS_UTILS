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

 program read_increments

 use read_write_data, only  :  read_data
 use mpi

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
 data stc_inc_expected_values_tile1 / -0.6302, -0.1116, 0.0341, 0.0 /
 data stc_inc_expected_values_tile2 / 0.0825,  0.0071, -0.0255, 0.0 /
 data stc_inc_expected_values_tile3 / 0.2070, 0.0608, 0.0001, 0.0 /
 data stc_inc_expected_values_tile4 / 0.0, 0.0, 0.0, 0.0 /
 data stc_inc_expected_values_tile5 / -0.1031, -0.0386, -0.0356, 0.0/
 data stc_inc_expected_values_tile6 / 0.0, 0.0, 0.0, 0.0 /
 data slc_inc_expected_values_tile1 / -0.0007285, 0.0000055, 0.0000003, 0.0 /
 data slc_inc_expected_values_tile2 / -0.0006, -0.0059, -0.0087, 0.0 /
 data slc_inc_expected_values_tile3 / -0.0015, 0.0030, -0.0007, 0.0 /
 data slc_inc_expected_values_tile4 / 0.0, 0.0, 0.0, 0.0 / 
 data slc_inc_expected_values_tile5 / 0.0014, 0.0012, 0.0006, 0.0 /
 data slc_inc_expected_values_tile6 / 0.0, 0.0, 0.0, 0.0 /

 call mpi_init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

 lsoil=4
 lensfc=36864
 ij_input = 28541

 allocate(SLCINC(lensfc,lsoil))
 allocate(STCINC(lensfc,lsoil))

 if (my_rank .eq. 0) print*,"Starting test of global_cycle routine read_data."
 if (my_rank .eq. 0) print*,"Call routine read_data"

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

