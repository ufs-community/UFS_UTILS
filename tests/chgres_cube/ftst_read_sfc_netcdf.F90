 program readsfc

! Unit test for the read_input_sfc_netcdf_file routine. 
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

 use input_data, only : read_input_sfc_data, &
                        landsea_mask_input_grid, &
                        terrain_input_grid

 use program_setup, only : input_type, &
                           data_dir_input_grid, &
                           sfc_files_input_grid

 implicit none

 integer, parameter :: NUM_VALUES=2

 real, parameter :: EPSILON=0.0001

 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_vm)                :: vm

 integer :: rc, localpet, npets

 real(esmf_kind_r8), allocatable    :: data_one_tile(:,:)

! The expected values were determined by the checking
! the input NetCDF history file using 'ncdump'.
! The two values are the first and last corner points.

 real :: landsea_mask_expected_values(NUM_VALUES)
 real :: terrain_expected_values(NUM_VALUES)

 data landsea_mask_expected_values /1.0, 1.0/
 data terrain_expected_values /319.5036, 405.1155/

 print*,"Starting test of read_input_sfc_netcdf_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 i_input = 9
 j_input = 9

 input_type = "gaussian_netcdf"
 num_tiles_input_grid = 1
 data_dir_input_grid = "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/tests/chgres_cube/data/"
 sfc_files_input_grid(1) = "gfs.v16.sfc2.history.nc"

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

 call read_input_sfc_data(localpet)
 
 allocate(data_one_tile(i_input,j_input))

 call ESMF_FieldGather(landsea_mask_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'mask ',data_one_tile
   if (abs(data_one_tile(1,1) - landsea_mask_expected_values(1)) > EPSILON) stop 2
   if (abs(data_one_tile(i_input,j_input) - landsea_mask_expected_values(2)) > EPSILON) stop 4
 endif

 call ESMF_FieldGather(terrain_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'terrain ',data_one_tile
   if (abs(data_one_tile(1,1) - terrain_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile(i_input,j_input) - terrain_expected_values(2)) > EPSILON) stop 8
 endif

 print*,"OK"

 deallocate(data_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program readsfc
