 program readsfcnetcdf

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

 use atm_input_data, only : terrain_input_grid

 use sfc_input_data, only : read_input_sfc_data, &
                        lsoil_input, &
                        landsea_mask_input_grid, &
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
                        skin_temp_input_grid, &
                        canopy_mc_input_grid, &
                        z0_input_grid

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

 real :: landsea_mask_expected_values(NUM_VALUES) ! land-sea mask
 real :: terrain_expected_values(NUM_VALUES) ! terrain height
 real :: soilm_liq1_expected_values(NUM_VALUES) ! layer 1 liquid soil moisture
 real :: soilm_liq2_expected_values(NUM_VALUES) ! layer 2 liquid soil moisture
 real :: soilm_liq3_expected_values(NUM_VALUES) ! layer 3 liquid soil moisture
 real :: soilm_liq4_expected_values(NUM_VALUES) ! layer 4 liquid soil moisture
 real :: soilm_tot1_expected_values(NUM_VALUES) ! layer 1 total soil moisture
 real :: soilm_tot2_expected_values(NUM_VALUES) ! layer 2 total soil moisture
 real :: soilm_tot3_expected_values(NUM_VALUES) ! layer 3 total soil moisture
 real :: soilm_tot4_expected_values(NUM_VALUES) ! layer 4 total soil moisture
 real :: soil_temp1_expected_values(NUM_VALUES) ! layer 1 soil temperature
 real :: soil_temp2_expected_values(NUM_VALUES) ! layer 2 soil temperature
 real :: soil_temp3_expected_values(NUM_VALUES) ! layer 3 soil temperature
 real :: soil_temp4_expected_values(NUM_VALUES) ! layer 4 soil temperature
 real :: seaice_fract_expected_values(NUM_VALUES) ! sea ice fraction
 real :: seaice_depth_expected_values(NUM_VALUES) ! sea ice depth
 real :: seaice_skin_temp_expected_values(NUM_VALUES) ! sea ice skin temperature
 real :: snow_liq_equiv_expected_values(NUM_VALUES) ! liquid equivalent snow depth
 real :: snow_depth_expected_values(NUM_VALUES) ! physical snow depth
 real :: veg_type_expected_values(NUM_VALUES) ! vegetation type
 real :: soil_type_expected_values(NUM_VALUES) ! soil type
 real :: t2m_expected_values(NUM_VALUES) ! two-meter temperature
 real :: q2m_expected_values(NUM_VALUES) ! two-meter specific humidity
 real :: tprcp_expected_values(NUM_VALUES) ! precipitation
 real :: f10m_expected_values(NUM_VALUES) ! log((z0+10)*l/z0)
 real :: ffmm_expected_values(NUM_VALUES) ! log((z0+z1)*l/z0)
 real :: ustar_expected_values(NUM_VALUES) ! friction velocity
 real :: srflag_expected_values(NUM_VALUES) ! snow/rain flag
 real :: skin_temp_expected_values(NUM_VALUES) ! skin temperature
 real :: canopy_mc_expected_values(NUM_VALUES) ! canopy moisture content
 real :: z0_expected_values(NUM_VALUES) ! roughness length

 data landsea_mask_expected_values /1.0, 1.0/
 data terrain_expected_values /319.5036, 405.1155/
 data soilm_liq1_expected_values / 0.33374, 0.33583/
 data soilm_liq2_expected_values / 0.33563, 0.33659/
 data soilm_liq3_expected_values / 0.33697, 0.33554/
 data soilm_liq4_expected_values / 0.33997, 0.33849/
 data soilm_tot1_expected_values / 0.33374, 0.33583/
 data soilm_tot2_expected_values / 0.33563, 0.33659/
 data soilm_tot3_expected_values / 0.33697, 0.33554/
 data soilm_tot4_expected_values / 0.33997, 0.33849/
 data soil_temp1_expected_values / 278.8724, 278.7153/
 data soil_temp2_expected_values / 279.4426, 280.144/
 data soil_temp3_expected_values / 281.495, 282.1981/
 data soil_temp4_expected_values / 284.2445, 285.1065/
 data seaice_fract_expected_values / 0.0, 0.0/
 data seaice_depth_expected_values / 0.0, 0.0/
 data seaice_skin_temp_expected_values / 277.2036, 277.8394/
 data snow_liq_equiv_expected_values / 0.0, 0.0 /
 data snow_depth_expected_values / 0.0, 0.0/
 data veg_type_expected_values / 14.0, 4.0 /
 data soil_type_expected_values / 4.0, 4.0 /
 data t2m_expected_values / 277.7575, 278.1006 /
 data q2m_expected_values / 0.0043968, 0.0050449/
 data tprcp_expected_values / 1.627e-06, 2.596e-07/
 data f10m_expected_values / 0.99733, 0.99539/
 data ffmm_expected_values / 3.90207, 2.97655/
 data ustar_expected_values / 0.44228, 0.27333/
 data srflag_expected_values / 0.0, 0.0/
 data skin_temp_expected_values / 277.2036, 277.8394/
 data canopy_mc_expected_values / 0.25392, 0.44883/
 data z0_expected_values / 0.25, 0.826/

 print*,"Starting test of read_input_sfc_netcdf_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 i_input = 9
 j_input = 9

 input_type = "gaussian_netcdf"
 num_tiles_input_grid = 1
 data_dir_input_grid = "data/"
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

 call ESMF_FieldGather(terrain_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - terrain_expected_values(1)) > EPSILON) stop 6
   if (abs(data_one_tile(i_input,j_input) - terrain_expected_values(2)) > EPSILON) stop 8
 endif

 call ESMF_FieldGather(soilm_liq_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile_3d(1,1,1) - soilm_liq1_expected_values(1)) > EPSILON) stop 10
   if (abs(data_one_tile_3d(i_input,j_input,1) - soilm_liq1_expected_values(2)) > EPSILON) stop 11
   if (abs(data_one_tile_3d(1,1,2) - soilm_liq2_expected_values(1)) > EPSILON) stop 12
   if (abs(data_one_tile_3d(i_input,j_input,2) - soilm_liq2_expected_values(2)) > EPSILON) stop 13
   if (abs(data_one_tile_3d(1,1,3) - soilm_liq3_expected_values(1)) > EPSILON) stop 14
   if (abs(data_one_tile_3d(i_input,j_input,3) - soilm_liq3_expected_values(2)) > EPSILON) stop 15
   if (abs(data_one_tile_3d(1,1,4) - soilm_liq4_expected_values(1)) > EPSILON) stop 16
   if (abs(data_one_tile_3d(i_input,j_input,4) - soilm_liq4_expected_values(2)) > EPSILON) stop 17
 endif

 call ESMF_FieldGather(soilm_tot_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile_3d(1,1,1) - soilm_tot1_expected_values(1)) > EPSILON) stop 20
   if (abs(data_one_tile_3d(i_input,j_input,1) - soilm_tot1_expected_values(2)) > EPSILON) stop 21
   if (abs(data_one_tile_3d(1,1,2) - soilm_tot2_expected_values(1)) > EPSILON) stop 22
   if (abs(data_one_tile_3d(i_input,j_input,2) - soilm_tot2_expected_values(2)) > EPSILON) stop 23
   if (abs(data_one_tile_3d(1,1,3) - soilm_tot3_expected_values(1)) > EPSILON) stop 24
   if (abs(data_one_tile_3d(i_input,j_input,3) - soilm_tot3_expected_values(2)) > EPSILON) stop 25
   if (abs(data_one_tile_3d(1,1,4) - soilm_tot4_expected_values(1)) > EPSILON) stop 26
   if (abs(data_one_tile_3d(i_input,j_input,4) - soilm_tot4_expected_values(2)) > EPSILON) stop 27
 endif

 call ESMF_FieldGather(soil_temp_input_grid, data_one_tile_3d, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile_3d(1,1,1) - soil_temp1_expected_values(1)) > EPSILON) stop 30
   if (abs(data_one_tile_3d(i_input,j_input,1) - soil_temp1_expected_values(2)) > EPSILON) stop 31
   if (abs(data_one_tile_3d(1,1,2) - soil_temp2_expected_values(1)) > EPSILON) stop 32
   if (abs(data_one_tile_3d(i_input,j_input,2) - soil_temp2_expected_values(2)) > EPSILON) stop 33
   if (abs(data_one_tile_3d(1,1,3) - soil_temp3_expected_values(1)) > EPSILON) stop 34
   if (abs(data_one_tile_3d(i_input,j_input,3) - soil_temp3_expected_values(2)) > EPSILON) stop 35
   if (abs(data_one_tile_3d(1,1,4) - soil_temp4_expected_values(1)) > EPSILON) stop 36
   if (abs(data_one_tile_3d(i_input,j_input,4) - soil_temp4_expected_values(2)) > EPSILON) stop 37
 endif

 call ESMF_FieldGather(landsea_mask_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - landsea_mask_expected_values(1)) > EPSILON) stop 39
   if (abs(data_one_tile(i_input,j_input) - landsea_mask_expected_values(2)) > EPSILON) stop 40
 endif

 call ESMF_FieldGather(seaice_fract_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - seaice_fract_expected_values(1)) > EPSILON) stop 41
   if (abs(data_one_tile(i_input,j_input) - seaice_fract_expected_values(2)) > EPSILON) stop 42
 endif

 call ESMF_FieldGather(seaice_depth_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - seaice_depth_expected_values(1)) > EPSILON) stop 43
   if (abs(data_one_tile(i_input,j_input) - seaice_depth_expected_values(2)) > EPSILON) stop 44
 endif

 call ESMF_FieldGather(seaice_skin_temp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - seaice_skin_temp_expected_values(1)) > EPSILON) stop 45
   if (abs(data_one_tile(i_input,j_input) - seaice_skin_temp_expected_values(2)) > EPSILON) stop 46
 endif

 call ESMF_FieldGather(snow_liq_equiv_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - snow_liq_equiv_expected_values(1)) > EPSILON) stop 47
   if (abs(data_one_tile(i_input,j_input) - snow_liq_equiv_expected_values(2)) > EPSILON) stop 48
 endif

 call ESMF_FieldGather(snow_depth_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - snow_depth_expected_values(1)) > EPSILON) stop 49
   if (abs(data_one_tile(i_input,j_input) - snow_depth_expected_values(2)) > EPSILON) stop 50
 endif

 call ESMF_FieldGather(veg_type_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - veg_type_expected_values(1)) > EPSILON) stop 51
   if (abs(data_one_tile(i_input,j_input) - veg_type_expected_values(2)) > EPSILON) stop 52
 endif

 call ESMF_FieldGather(soil_type_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - soil_type_expected_values(1)) > EPSILON) stop 53
   if (abs(data_one_tile(i_input,j_input) - soil_type_expected_values(2)) > EPSILON) stop 54
 endif

 call ESMF_FieldGather(t2m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - t2m_expected_values(1)) > EPSILON) stop 55
   if (abs(data_one_tile(i_input,j_input) - t2m_expected_values(2)) > EPSILON) stop 56
 endif

 call ESMF_FieldGather(q2m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - q2m_expected_values(1)) > EPSILON) stop 57
   if (abs(data_one_tile(i_input,j_input) - q2m_expected_values(2)) > EPSILON) stop 58
 endif

 call ESMF_FieldGather(tprcp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - tprcp_expected_values(1)) > EPSILON) stop 59
   if (abs(data_one_tile(i_input,j_input) - tprcp_expected_values(2)) > EPSILON) stop 60
 endif

 call ESMF_FieldGather(f10m_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - f10m_expected_values(1)) > EPSILON) stop 61
   if (abs(data_one_tile(i_input,j_input) - f10m_expected_values(2)) > EPSILON) stop 62
 endif

 call ESMF_FieldGather(ffmm_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - ffmm_expected_values(1)) > EPSILON) stop 63
   if (abs(data_one_tile(i_input,j_input) - ffmm_expected_values(2)) > EPSILON) stop 64
 endif

 call ESMF_FieldGather(ustar_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - ustar_expected_values(1)) > EPSILON) stop 65
   if (abs(data_one_tile(i_input,j_input) - ustar_expected_values(2)) > EPSILON) stop 66
 endif

 call ESMF_FieldGather(srflag_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - srflag_expected_values(1)) > EPSILON) stop 67
   if (abs(data_one_tile(i_input,j_input) - srflag_expected_values(2)) > EPSILON) stop 68
 endif

 call ESMF_FieldGather(skin_temp_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - skin_temp_expected_values(1)) > EPSILON) stop 69
   if (abs(data_one_tile(i_input,j_input) - skin_temp_expected_values(2)) > EPSILON) stop 70
 endif

 call ESMF_FieldGather(canopy_mc_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - canopy_mc_expected_values(1)) > EPSILON) stop 71
   if (abs(data_one_tile(i_input,j_input) - canopy_mc_expected_values(2)) > EPSILON) stop 72
 endif

 call ESMF_FieldGather(z0_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (localpet == 0) then
   if (abs(data_one_tile(1,1) - z0_expected_values(1)) > EPSILON) stop 73
   if (abs(data_one_tile(i_input,j_input) - z0_expected_values(2)) > EPSILON) stop 74
 endif

 print*,"OK"

 deallocate(data_one_tile, data_one_tile_3d)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program readsfcnetcdf
