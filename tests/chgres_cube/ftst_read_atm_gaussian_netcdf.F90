 program read_atm_gaussian_netcdf

! Unit test for the "read_input_atm_gaussian_netcdf_file"
! routine of chgres_cube.
!
! Reads a 9x9 version of the GFS v16 NetCDF atmospheric history 
! file. This smaller version of the full file was created
! using the 'ncks' utility. The full file is too big for
! Github. The data read from the file is compared to
! expected values as determined from the 'ncdump' utility.
!
! Author George Gayno

 use esmf

 use model_grid, only : i_input, j_input, &
                        input_grid, &
                        num_tiles_input_grid, &
                        latitude_input_grid, &
                        longitude_input_grid

 use atm_input_data, only : read_input_atm_data, &
                        lev_input, &
                        levp1_input, &
                        temp_input_grid, &
                        dzdt_input_grid, &
                        ps_input_grid, &
                        pres_input_grid, &
                        xwind_input_grid, &
                        ywind_input_grid, &
                        zwind_input_grid, &
                        terrain_input_grid, &
                        tracers_input_grid

 use program_setup, only : input_type, &
                           data_dir_input_grid, &
                           num_tracers_input, &
                           tracers_input, &
                           num_tracers, &
                           atm_files_input_grid

 implicit none

 integer, parameter :: EXPECTED_LEV_INPUT=127 ! Number of vertical layers.
 integer, parameter :: EXPECTED_LEVP1_INPUT=128 ! Number of vertical layer
                                                ! interfaces.

 integer, parameter :: NUM_VALUES=2

 real, parameter :: EPSILON=0.0001
 real, parameter :: EPSILON_SMALL=0.0000001

 real :: expected_values_tmp(NUM_VALUES)
 real :: expected_values_spfh(NUM_VALUES)
 real :: expected_values_clwmr(NUM_VALUES)
 real :: expected_values_o3mr(NUM_VALUES)
 real :: expected_values_dzdt(NUM_VALUES)
 real :: expected_values_ps(NUM_VALUES)
 real :: expected_values_pres(NUM_VALUES)
 real :: expected_values_xwind(NUM_VALUES)
 real :: expected_values_ywind(NUM_VALUES)
 real :: expected_values_zwind(NUM_VALUES)
 real :: expected_values_terrain(NUM_VALUES)

 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_vm)                :: vm

 integer :: n, rc, localpet, npets

 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable  :: data3d_one_tile(:,:,:)
 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)

! The expected values were determined by the checking
! the input NetCDF history file using 'ncdump'.

 data expected_values_tmp / 301.6022, 182.6277 / ! Temperature
 data expected_values_spfh / 0.01331, 3.754e-06 / ! Specific humidity
 data expected_values_clwmr / 0.0, 0.0 / ! Cloud liquid water
 data expected_values_o3mr / 6.2015e-08, 2.632e-07 / ! Ozone
 data expected_values_dzdt / 0.0, 0.0 / ! Vertical velocity
 data expected_values_ps / 100971.7454, 100941.9350 / ! Surface pressure in
                                                      ! Pascals.
 data expected_values_pres / 100846.9917, 1.30199 / ! 3-d pressure in Pascals.
 data expected_values_xwind / -5.9860, -10.2035 / ! 'x' component wind
 data expected_values_ywind / -2.8133, -0.99638 / ! 'y' component wind
 data expected_values_zwind / 0.0, 0.0 / ! 'z' component wind
 data expected_values_terrain / 0.0, 0.0 / ! Terrain height.

 print*,"Starting test of read_input_atm_gaussian_netcdf_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 i_input = 9
 j_input = 9

 num_tracers=3
 num_tracers_input=num_tracers
 tracers_input(1)="spfh"
 tracers_input(2)="clwmr"
 tracers_input(3)="o3mr"

 input_type = "gaussian_netcdf"
 num_tiles_input_grid = 1
 data_dir_input_grid = "data/"
 atm_files_input_grid(1) = "gfs.v16.atm.history.nc"

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)

 latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_latitude", &
                                   rc=rc)

 longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="input_grid_longitude", &
                                   rc=rc)

! Using the north pole/greenwich simplifies the comparison
! of u/v winds with the x/y/z components.

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))
 latitude=90.0
 longitude=0.0

 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)

! Call the chgres_cube read routine.

 call read_input_atm_data(localpet)

 if (lev_input /= EXPECTED_LEV_INPUT) stop 2
 if (levp1_input /= EXPECTED_LEVP1_INPUT) stop 3

 allocate(data3d_one_tile(i_input,j_input,lev_input))
 allocate(data_one_tile(i_input,j_input))

 call ESMF_FieldGather(temp_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_tmp(1)) > EPSILON) stop 4
 if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_tmp(2)) > EPSILON) stop 5

 do n = 1, num_tracers
   call ESMF_FieldGather(tracers_input_grid(n), data3d_one_tile, rootPet=0, rc=rc)
   if (trim(tracers_input(n)) == 'spfh') then
     if (abs(data3d_one_tile(1,1,1) - expected_values_spfh(1)) > EPSILON) stop 6
     if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_spfh(2)) > EPSILON_SMALL) stop 7
   endif
   if (trim(tracers_input(n)) == 'clwmr') then
     if (abs(data3d_one_tile(1,1,1) - expected_values_clwmr(1)) > EPSILON) stop 8
     if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_clwmr(2)) > EPSILON) stop 9
   endif
   if (trim(tracers_input(n)) == 'o3mr') then
     if (abs(data3d_one_tile(1,1,1) - expected_values_o3mr(1)) > EPSILON_SMALL) stop 10
     if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_o3mr(2)) > EPSILON_SMALL) stop 11
   endif
 enddo

 call ESMF_FieldGather(dzdt_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_dzdt(1)) > EPSILON) stop 14
 if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_dzdt(2)) > EPSILON) stop 15

 call ESMF_FieldGather(ps_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (abs(data_one_tile(1,1) - expected_values_ps(1)) > EPSILON) stop 16
 if (abs(data_one_tile(i_input,j_input) - expected_values_ps(2)) > EPSILON) stop 17

 call ESMF_FieldGather(pres_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_pres(1)) > EPSILON) stop 18
 if (abs(data3d_one_tile(i_input,j_input,lev_input) - expected_values_pres(2)) > EPSILON) stop 19

 call ESMF_FieldGather(xwind_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_xwind(1)) > EPSILON) stop 20
 if (abs(data3d_one_tile(1,1,lev_input) - expected_values_xwind(2)) > EPSILON) stop 21
 call ESMF_FieldGather(ywind_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_ywind(1)) > EPSILON) stop 22
 if (abs(data3d_one_tile(1,1,lev_input) - expected_values_ywind(2)) > EPSILON) stop 23
 call ESMF_FieldGather(zwind_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(1,1,1) - expected_values_zwind(1)) > EPSILON) stop 24
 if (abs(data3d_one_tile(1,1,lev_input) - expected_values_zwind(2)) > EPSILON) stop 25

 call ESMF_FieldGather(terrain_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (abs(data_one_tile(1,1) - expected_values_terrain(1)) > EPSILON) stop 26
 if (abs(data_one_tile(i_input,j_input) - expected_values_terrain(2)) > EPSILON) stop 27

 print*,"OK"

 deallocate(latitude, longitude, data3d_one_tile, data_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program read_atm_gaussian_netcdf
