 program read_atm_grib2

 use esmf

 use input_data, only    :  read_input_atm_data, &
                            lev_input, &
                            levp1_input, &
                            temp_input_grid, tracers_input_grid, &
                            dzdt_input_grid, pres_input_grid, &
                            ps_input_grid

 use program_setup, only :  input_type, data_dir_input_grid, &
                            grib2_file_input_grid, &
                            read_varmap, varmap_file, &
                            num_tracers, num_tracers_input, &
                            tracers, tracers_input, external_model

 use model_grid, only    : input_grid, &
                           i_input, j_input, &
                           latitude_input_grid, &
                           longitude_input_grid

 implicit none

 type(esmf_vm)                :: vm

 type(esmf_polekind_flag)     :: polekindflag(2)

 integer, parameter :: EXPECTED_LEV_INPUT=31 ! Number of vertical layers.
 integer, parameter :: EXPECTED_LEVP1_INPUT=32 ! Number of vertical layer
                                                ! interfaces.

 integer, parameter :: NUM_VALUES=2

 real, parameter    :: EPSILON=0.0001
 real, parameter    :: EPSILON_SMALL=0.0000001

 integer            :: rc, localpet, n, npets
 integer            :: i_check(NUM_VALUES), j_check(NUM_VALUES), k_check(NUM_VALUES)

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)
 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)
 real(esmf_kind_r8), allocatable  :: data3d_one_tile(:,:,:)

 real :: expected_values_tmp(NUM_VALUES)
 real :: expected_values_sphum(NUM_VALUES)
 real :: expected_values_liq_wat(NUM_VALUES)
 real :: expected_values_o3mr(NUM_VALUES)
 real :: expected_values_icewat(NUM_VALUES)
 real :: expected_values_rainwat(NUM_VALUES)
 real :: expected_values_snowwat(NUM_VALUES)
 real :: expected_values_graupel(NUM_VALUES)
 real :: expected_values_dzdt(NUM_VALUES)
 real :: expected_values_pres(NUM_VALUES)
 real :: expected_values_ps(NUM_VALUES)

! The expected values were determined by the checking
! the input GRIB2 file using stand-alone 'wgrib2'.

 data expected_values_tmp / 300.5728, 262.8000 / ! Temperature
 data expected_values_sphum / 0.01659, 0.0 / ! Specific humidity
 data expected_values_liq_wat / 0.0, 0.0 / ! Cloud liquid water
 data expected_values_o3mr / 6.69e-08, 4.94e-06 / ! Ozone
 data expected_values_icewat / 0.0, 0.0 / ! Ice water
 data expected_values_rainwat / 0.0, 0.0 / ! Rain water
 data expected_values_snowwat / 0.0, 0.0 / ! Snow water
 data expected_values_graupel / 0.0, 0.0 / ! Graupel
 data expected_values_dzdt / 8.48e-03, 0.0 / ! Vertical velocity
 data expected_values_pres / 100000.0, 100.0 / ! 3-d pressure
 data expected_values_ps / 100662.9, 100662.9 / ! Surface pressure

 print*,"Starting test of read_atm_grib2_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 external_model="GFS"
 input_type="grib2"
 data_dir_input_grid = "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/gfs.grib2"
 grib2_file_input_grid = "gfs.t00z.pgrb2.0p50.f000"

 i_input = 720
 j_input = 361

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)


 num_tracers=7
 num_tracers_input=num_tracers

 tracers_input(1)="spfh"
 tracers_input(2)="clwmr"
 tracers_input(3)="o3mr"
 tracers_input(4)="icmr"
 tracers_input(5)="rwmr"
 tracers_input(6)="snmr"
 tracers_input(7)="grle"

 tracers(1)="sphum"
 tracers(2)="liq_wat"
 tracers(3)="o3mr"
 tracers(4)="ice_wat"
 tracers(5)="rainwat"
 tracers(6)="snowwat"
 tracers(7)="graupel"

 varmap_file ="/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/parm/varmap_tables/GFSphys_var_map.txt"
 call read_varmap

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
! of u/v winds.

 allocate(latitude(i_input,j_input))
 allocate(longitude(i_input,j_input))
 latitude=90.0
 longitude=0.0

 call ESMF_FieldScatter(longitude_input_grid, longitude, rootpet=0, rc=rc)
 call ESMF_FieldScatter(latitude_input_grid, latitude, rootpet=0, rc=rc)

 call read_input_atm_data(localpet)

 if (lev_input /= EXPECTED_LEV_INPUT) stop 2
 if (levp1_input /= EXPECTED_LEVP1_INPUT) stop 3

 allocate(data3d_one_tile(i_input,j_input,lev_input))
 allocate(data_one_tile(i_input,j_input))

 i_check(1) = i_input/2
 j_check(1) = j_input/2
 k_check(1) = 1

 i_check(2) = i_input/2
 j_check(2) = j_input/2
 k_check(2) = lev_input

 call ESMF_FieldGather(temp_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_tmp(1)) > EPSILON) stop 4
 if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_tmp(2)) > EPSILON) stop 5

 do n = 1, num_tracers
   call ESMF_FieldGather(tracers_input_grid(n), data3d_one_tile, rootPet=0, rc=rc)
   if (trim(tracers(n)) == 'sphum') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_sphum(1)) > EPSILON) stop 6
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_sphum(2)) > EPSILON_SMALL) stop 7
   endif
   if (trim(tracers(n)) == 'liq_wat') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_liq_wat(1)) > EPSILON_SMALL) stop 8
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_liq_wat(2)) > EPSILON_SMALL) stop 9
   endif
   if (trim(tracers(n)) == 'o3mr') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_o3mr(1)) > EPSILON_SMALL) stop 10
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_o3mr(2)) > EPSILON_SMALL) stop 11
   endif
   if (trim(tracers(n)) == 'ice_wat') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_icewat(1)) > EPSILON_SMALL) stop 12
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_icewat(2)) > EPSILON_SMALL) stop 13
   endif
   if (trim(tracers(n)) == 'rainwat') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_rainwat(1)) > EPSILON_SMALL) stop 14
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_rainwat(2)) > EPSILON_SMALL) stop 15
   endif
   if (trim(tracers(n)) == 'snowwat') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_snowwat(1)) > EPSILON_SMALL) stop 16
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_snowwat(2)) > EPSILON_SMALL) stop 17
   endif
   if (trim(tracers(n)) == 'graupel') then
     if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_graupel(1)) > EPSILON_SMALL) stop 18
     if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_graupel(2)) > EPSILON_SMALL) stop 19
   endif
 enddo

 call ESMF_FieldGather(dzdt_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_dzdt(1)) > EPSILON) stop 20
 if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_dzdt(2)) > EPSILON) stop 21

 call ESMF_FieldGather(pres_input_grid, data3d_one_tile, rootPet=0, rc=rc)
 if (abs(data3d_one_tile(i_check(1),j_check(1),k_check(1)) - expected_values_pres(1)) > EPSILON) stop 22
 if (abs(data3d_one_tile(i_check(2),j_check(2),k_check(2)) - expected_values_pres(2)) > EPSILON) stop 23

 call ESMF_FieldGather(ps_input_grid, data_one_tile, rootPet=0, rc=rc)
 if (abs(data_one_tile(i_check(1),j_check(1)) - expected_values_ps(1)) > EPSILON) stop 24

 deallocate(latitude, longitude, data3d_one_tile, data_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program read_atm_grib2
