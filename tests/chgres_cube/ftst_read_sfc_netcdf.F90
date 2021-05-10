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
                        lsoil_input, &
                        landsea_mask_input_grid, &
                        terrain_input_grid, &
                        soilm_liq_input_grid, &
                        soilm_tot_input_grid, &
                        soil_temp_input_grid, &
                        seaice_fract_input_grid, &
                        seaice_depth_input_grid, &
                        seaice_skin_temp_input_grid, &
                        snow_liq_equiv_input_grid, &
                        snow_depth_input_grid, &
                        veg_type_input_grid, &
                        soil_type_input_grid, &
                        t2m_input_grid, &
                        q2m_input_grid, &
                        tprcp_input_grid, &
                        f10m_input_grid, &
                        ffmm_input_grid, &
                        ustar_input_grid, &
                        srflag_input_grid, &
                        skin_temp_input_grid

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
 real(esmf_kind_r8), allocatable    :: data_one_tile_3d(:,:,:)

! The expected values were determined by the checking
! the input NetCDF history file using 'ncdump'.
! The two values are the first and last corner points.

 real :: landsea_mask_expected_values(NUM_VALUES)
 real :: terrain_expected_values(NUM_VALUES)
 real :: soilm_liq1_expected_values(NUM_VALUES)
 real :: soilm_liq2_expected_values(NUM_VALUES)
 real :: soilm_liq3_expected_values(NUM_VALUES)
 real :: soilm_liq4_expected_values(NUM_VALUES)

 data landsea_mask_expected_values /1.0, 1.0/
 data terrain_expected_values /319.5036, 405.1155/
 data soilm_liq1_expected_values / 0.33374, 0.33583/
 data soilm_liq2_expected_values / 0.33563, 0.33659/
 data soilm_liq3_expected_values / 0.33697, 0.33554/
 data soilm_liq4_expected_values / 0.33997, 0.33849/

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
 allocate(data_one_tile_3d(i_input,j_input,lsoil_input))

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

 call ESMF_FieldGather(soilm_liq_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'slc 1 ',data_one_tile_3d(1,1,:)
   print*,'slc last ',data_one_tile_3d(i_input,j_input,:)
   if (abs(data_one_tile_3d(1,1,1) - soilm_liq1_expected_values(1)) > EPSILON) stop 10
   if (abs(data_one_tile_3d(i_input,j_input,1) - soilm_liq1_expected_values(2)) > EPSILON) stop 11
   if (abs(data_one_tile_3d(1,1,2) - soilm_liq2_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile_3d(i_input,j_input,2) - soilm_liq2_expected_values(2)) > EPSILON) stop 12
   if (abs(data_one_tile_3d(1,1,3) - soilm_liq3_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile_3d(i_input,j_input,3) - soilm_liq3_expected_values(2)) > EPSILON) stop 13
   if (abs(data_one_tile_3d(1,1,4) - soilm_liq4_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile_3d(i_input,j_input,4) - soilm_liq4_expected_values(2)) > EPSILON) stop 14
 endif

 call ESMF_FieldGather(soilm_tot_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'slm 1 ',data_one_tile_3d(1,1,:)
   print*,'slm last ',data_one_tile_3d(i_input,j_input,:)
 endif

 call ESMF_FieldGather(soil_temp_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'stc 1 ',data_one_tile_3d(1,1,:)
   print*,'stc last ',data_one_tile_3d(i_input,j_input,:)
 endif

 call ESMF_FieldGather(seaice_fract_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'icef ',data_one_tile
 endif

 call ESMF_FieldGather(seaice_depth_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'iced ',data_one_tile
 endif

 call ESMF_FieldGather(seaice_skin_temp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'icet ',data_one_tile
 endif

 call ESMF_FieldGather(snow_liq_equiv_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'weasd ',data_one_tile
 endif

 call ESMF_FieldGather(snow_depth_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'snod ',data_one_tile
 endif

 call ESMF_FieldGather(veg_type_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'vegt ',data_one_tile
 endif

 call ESMF_FieldGather(soil_type_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'soilt ',data_one_tile
 endif

 call ESMF_FieldGather(t2m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'t2m ',data_one_tile
 endif

 call ESMF_FieldGather(q2m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'q2m ',data_one_tile
 endif

 call ESMF_FieldGather(tprcp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'tprcp ',data_one_tile
 endif

 call ESMF_FieldGather(f10m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'f10m ',data_one_tile
 endif

 call ESMF_FieldGather(ffmm_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'ffmm ',data_one_tile
 endif

 call ESMF_FieldGather(ustar_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'ustar ',data_one_tile
 endif

 call ESMF_FieldGather(srflag_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'srflag ',data_one_tile
 endif

 call ESMF_FieldGather(skin_temp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   print*,'skin t ',data_one_tile
 endif

 print*,"OK"

 deallocate(data_one_tile, data_one_tile_3d)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program readsfc
