 program readnst

! Unit test for the read_input_nst_netcdf_file routine. 
!
! Reads a 9x9 version of the GFS v16 NetCDF surface history 
! file. This smaller version of the full file was created
! using the 'ncks' utility. The full file is too big for
! Github. The data read from the file is compared to
! expected values as determined from the 'ncdump' utility.
!
! Author George Gayno

 use esmf

 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        num_tiles_input_grid

 use input_data, only : read_input_nst_data, &
                        c_d_input_grid, &
                        c_0_input_grid, &
                        d_conv_input_grid, &
                        dt_cool_input_grid, &
                        ifd_input_grid, &
                        qrain_input_grid, &
                        tref_input_grid, &
                        w_d_input_grid, &
                        w_0_input_grid, &
                        xs_input_grid, &
                        xt_input_grid, &
                        xu_input_grid, &
                        xv_input_grid, &
                        xz_input_grid, &
                        xtts_input_grid, &
                        xzts_input_grid, &
                        z_c_input_grid, &
                        zm_input_grid

 use program_setup, only : input_type, &
                           data_dir_input_grid, &
                           sfc_files_input_grid

 implicit none

 integer, parameter :: NUM_VALUES=2

 real, parameter :: EPSILON=0.00001

 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_vm)                :: vm

 integer :: rc, localpet, npets

 real(esmf_kind_r8), allocatable    :: data_one_tile(:,:)

! The expected values were determined by the checking
! the input NetCDF history file using 'ncdump'.
! The two values are the first and last corner points.

 real :: c_d_expected_values(NUM_VALUES)
 real :: c_0_expected_values(NUM_VALUES)
 real :: d_conv_expected_values(NUM_VALUES)
 real :: dt_cool_expected_values(NUM_VALUES)
 real :: ifd_expected_values(NUM_VALUES)
 real :: qrain_expected_values(NUM_VALUES)
 real :: tref_expected_values(NUM_VALUES)
 real :: wd_expected_values(NUM_VALUES)
 real :: w0_expected_values(NUM_VALUES)
 real :: xs_expected_values(NUM_VALUES)
 real :: xt_expected_values(NUM_VALUES)
 real :: xtts_expected_values(NUM_VALUES)
 real :: xu_expected_values(NUM_VALUES)
 real :: xv_expected_values(NUM_VALUES)
 real :: xz_expected_values(NUM_VALUES)
 real :: xzts_expected_values(NUM_VALUES)
 real :: zc_expected_values(NUM_VALUES)
 real :: zm_expected_values(NUM_VALUES)

 data c_d_expected_values / -76.08778, -83.7614/
 data c_0_expected_values / 0.0638366, 0.0646863 /
 data dt_cool_expected_values / 0.412696, 0.3639118 /
 data d_conv_expected_values / 0.0, 0.0/
 data ifd_expected_values / 1.0, 1.0 /
 data qrain_expected_values / 0.0, 0.003217 /
 data tref_expected_values / 302.315, 302.54385 /
 data wd_expected_values / 0.001898, 0.0 /
 data w0_expected_values / -0.000331, 0.0 /
 data xs_expected_values / 0.02394, 0.0 /
 data xt_expected_values / 0.466102, 0.0 /
 data xtts_expected_values / 0.05393, 0.0 /
 data xu_expected_values / -0.22587, 0.0 /
 data xv_expected_values / -0.094777, 0.0 /
 data xz_expected_values / 7.71445, 30.0 /
 data xzts_expected_values / 0.913822, 0.0 /
 data zc_expected_values / 0.0008382, 0.0007719 /
 data zm_expected_values / 0.0, 0.0 /

 print*,"Starting test of read_input_nst_netcdf_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 i_input = 9
 j_input = 9

 input_type = "gaussian_netcdf"
 num_tiles_input_grid = 1
 data_dir_input_grid = "data/"
 sfc_files_input_grid(1) = "gfs.v16.sfc.history.nc"

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

! Call the chgres_cube read routine.

 call read_input_nst_data(localpet)
 
 allocate(data_one_tile(i_input,j_input))

 call ESMF_FieldGather(c_d_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - c_d_expected_values(1)) > EPSILON) stop 2
   if (abs(data_one_tile(i_input,j_input) - c_d_expected_values(2)) > EPSILON) stop 4
 endif
 call ESMF_FieldGather(c_0_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - c_0_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile(i_input,j_input) - c_0_expected_values(2)) > EPSILON) stop 8
 endif
 call ESMF_FieldGather(d_conv_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - d_conv_expected_values(1)) > EPSILON) stop 10
   if (abs(data_one_tile(i_input,j_input) - d_conv_expected_values(2)) > EPSILON) stop 12
 endif
 call ESMF_FieldGather(dt_cool_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - dt_cool_expected_values(1)) > EPSILON) stop 14
   if (abs(data_one_tile(i_input,j_input) - dt_cool_expected_values(2)) > EPSILON) stop 16
 endif
 call ESMF_FieldGather(ifd_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - ifd_expected_values(1)) > EPSILON) stop 18
   if (abs(data_one_tile(i_input,j_input) - ifd_expected_values(2)) > EPSILON) stop 20
 endif
 call ESMF_FieldGather(qrain_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - qrain_expected_values(1)) > EPSILON) stop 22
   if (abs(data_one_tile(i_input,j_input) - qrain_expected_values(2)) > EPSILON) stop 24
 endif
 call ESMF_FieldGather(tref_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - tref_expected_values(1)) > EPSILON) stop 26
   if (abs(data_one_tile(i_input,j_input) - tref_expected_values(2)) > EPSILON) stop 28
 endif
 call ESMF_FieldGather(w_d_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - wd_expected_values(1)) > EPSILON) stop 30
   if (abs(data_one_tile(i_input,j_input) - wd_expected_values(2)) > EPSILON) stop 32
 endif
 call ESMF_FieldGather(w_0_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - w0_expected_values(1)) > EPSILON) stop 34
   if (abs(data_one_tile(i_input,j_input) - w0_expected_values(2)) > EPSILON) stop 36
 endif
 call ESMF_FieldGather(xs_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xs_expected_values(1)) > EPSILON) stop 38
   if (abs(data_one_tile(i_input,j_input) - xs_expected_values(2)) > EPSILON) stop 40
 endif
 call ESMF_FieldGather(xt_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xt_expected_values(1)) > EPSILON) stop 42
   if (abs(data_one_tile(i_input,j_input) - xt_expected_values(2)) > EPSILON) stop 44
 endif
 call ESMF_FieldGather(xu_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xu_expected_values(1)) > EPSILON) stop 46
   if (abs(data_one_tile(i_input,j_input) - xu_expected_values(2)) > EPSILON) stop 48
 endif
 call ESMF_FieldGather(xv_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xv_expected_values(1)) > EPSILON) stop 50
   if (abs(data_one_tile(i_input,j_input) - xv_expected_values(2)) > EPSILON) stop 52
 endif
 call ESMF_FieldGather(xz_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xz_expected_values(1)) > EPSILON) stop 54
   if (abs(data_one_tile(i_input,j_input) - xz_expected_values(2)) > EPSILON) stop 56
 endif
 call ESMF_FieldGather(xtts_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xtts_expected_values(1)) > EPSILON) stop 58
   if (abs(data_one_tile(i_input,j_input) - xtts_expected_values(2)) > EPSILON) stop 60
 endif
 call ESMF_FieldGather(xzts_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - xzts_expected_values(1)) > EPSILON) stop 62
   if (abs(data_one_tile(i_input,j_input) - xzts_expected_values(2)) > EPSILON) stop 64
 endif
 call ESMF_FieldGather(z_c_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - zc_expected_values(1)) > EPSILON) stop 66
   if (abs(data_one_tile(i_input,j_input) - zc_expected_values(2)) > EPSILON) stop 68
 endif
 call ESMF_FieldGather(zm_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - zm_expected_values(1)) > EPSILON) stop 70
   if (abs(data_one_tile(i_input,j_input) - zm_expected_values(2)) > EPSILON) stop 72
 endif

 print*,"OK"

 deallocate(data_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program readnst
