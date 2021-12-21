 program read_atm_grib2

 use esmf

 use input_data, only    :  read_input_atm_data

 use program_setup, only :  input_type, data_dir_input_grid, &
                            grib2_file_input_grid, &
                            read_varmap, varmap_file

 use model_grid, only    : input_grid, &
                           i_input, j_input

 implicit none

 type(esmf_vm)                :: vm

 type(esmf_polekind_flag)     :: polekindflag(2)

 integer :: rc, localpet, npets

 print*,"Starting test of read_atm_grib2_file."

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

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

 call read_input_atm_data(localpet)


 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program read_atm_grib2
