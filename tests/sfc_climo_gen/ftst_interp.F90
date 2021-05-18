 program test_input

! Unit test for the 'interp' routine. Read and 
! interpolate a field of vegetation type to a 
! small FV3 regional grid. Check the result against
! known values.

 use program_setup, only : mosaic_file_mdl, &
                           orog_dir_mdl, &
                           orog_files_mdl
 use model_grid
 use source_grid
 use esmf

 implicit none

 integer, parameter :: NUM_VALUES=4

 character(len=150) :: input_vegetation_type_file

 integer :: localpet, npets, rc

 real(esmf_kind_r4), allocatable    :: vegt_mdl_one_tile(:,:)
 real(esmf_kind_r4) :: expected_values(NUM_VALUES)

 type(esmf_regridmethod_flag) :: method
 type(esmf_vm)                :: vm

 data expected_values /10.0, 8.0, 12.0, 14.0/

 print*,"Starting test of routine interp."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 mosaic_file_mdl="data/C450_mosaic.nc"
 orog_dir_mdl="data"
 orog_files_mdl="C450_oro_data.tile7.nc"

 call define_model_grid(localpet, npets)

 input_vegetation_type_file="data/vegetation_type.viirs.igbp.0.1.nc"
 call define_source_grid(localpet, npets, input_vegetation_type_file)

! Call sfc_climo_gen routine.

 method=ESMF_REGRIDMETHOD_NEAREST_STOD
 call interp(localpet, method, input_vegetation_type_file)

 allocate(vegt_mdl_one_tile(i_mdl,j_mdl))
 call ESMF_FieldGather(vegt_field_mdl, vegt_mdl_one_tile, rootPet=0, tile=1, rc=rc)

! Check against expected values.

 if (vegt_mdl_one_tile(10,24) /= expected_values(1)) stop 2
 if (vegt_mdl_one_tile(28,27) /= expected_values(2)) stop 3
 if (vegt_mdl_one_tile(17,4) /= expected_values(3)) stop 4
 if (vegt_mdl_one_tile(27,13) /= expected_values(4)) stop 5

 call source_grid_cleanup
 call model_grid_cleanup

 deallocate(vegt_mdl_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program test_input
