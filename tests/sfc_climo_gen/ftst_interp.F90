 program test_input

! Unit test for the 'interp' routine. Read and 
! interpolate a field of vegetation type and
! check against known values.

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

 print*,"- INITIALIZE ESMF"
 call ESMF_Initialize(rc=rc)

 print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=rc)

 print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 mosaic_file_mdl="/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/tests/sfc_climo_gen/data/C450_mosaic.nc"
 orog_dir_mdl="/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/tests/sfc_climo_gen/data"
 orog_files_mdl="C450_oro_data.tile7.nc"

 call define_model_grid(localpet, npets)

 input_vegetation_type_file="/scratch1/NCEPDEV/da/George.Gayno/noscrub/viirs.veg/vegetation_type.viirs.igbp.0.1.nc"
 call define_source_grid(localpet, npets, input_vegetation_type_file)

 method=ESMF_REGRIDMETHOD_NEAREST_STOD
 call interp(localpet, method, input_vegetation_type_file)

 allocate(vegt_mdl_one_tile(i_mdl,j_mdl))
 call ESMF_FieldGather(vegt_field_mdl, vegt_mdl_one_tile, rootPet=0, tile=1, rc=rc)

 if (vegt_mdl_one_tile(10,24) /= expected_values(1)) stop 2
 if (vegt_mdl_one_tile(28,27) /= expected_values(2)) stop 3
 if (vegt_mdl_one_tile(17,4) /= expected_values(3)) stop 4
 if (vegt_mdl_one_tile(27,13) /= expected_values(4)) stop 5

 deallocate(vegt_mdl_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program test_input
