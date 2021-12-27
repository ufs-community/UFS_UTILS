 program read_atm_grib2

 use esmf

 use input_data, only    :  read_input_atm_data, &
                            lev_input, &
                            levp1_input

 use program_setup, only :  input_type, data_dir_input_grid, &
                            grib2_file_input_grid, &
                            read_varmap, varmap_file, &
                            num_tracers, num_tracers_input, &
                            tracers_input, external_model

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

 integer :: rc, localpet, npets

 real(esmf_kind_r8), allocatable  :: latitude(:,:)
 real(esmf_kind_r8), allocatable  :: longitude(:,:)

 print*,"Starting test of read_atm_grib2_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 external_model="GFS"
 input_type="grib2"
 data_dir_input_grid = "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/gfs.grib2"
 grib2_file_input_grid = "gfs.t00z.pgrb2.0p50.f000"

 varmap_file ="/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/parm/varmap_tables/GFSphys_var_map.txt"
 call read_varmap

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


 num_tracers=3
 num_tracers_input=num_tracers
 tracers_input(1)="spfh"
 tracers_input(2)="clwmr"
 tracers_input(3)="o3mr"

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

 deallocate(latitude, longitude)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program read_atm_grib2
