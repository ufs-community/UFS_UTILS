 program readnst

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

 type(esmf_polekind_flag)     :: polekindflag(2)
 type(esmf_vm)                :: vm

 integer :: rc, localpet, npets

 real(esmf_kind_r8), allocatable    :: data_one_tile(:,:)

 call mpi_init(rc)

 call ESMF_Initialize(rc=rc)

 call ESMF_VMGetGlobal(vm, rc=rc)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=rc)

 i_input = 9
 j_input = 9

 input_type = "gaussian_netcdf"
 num_tiles_input_grid = 1
 data_dir_input_grid = "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/tests/chgres_cube/data/"
 sfc_files_input_grid(1) = "sfc.nc"

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
 input_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                   maxIndex=(/i_input,j_input/), &
                                   polekindflag=polekindflag, &
                                   periodicDim=1, &
                                   poleDim=2,  &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, &
                                   regDecomp=(/1,npets/),  &
                                   indexflag=ESMF_INDEX_GLOBAL, rc=rc)
 print*,'rc ',rc

 call read_input_nst_data(localpet)
 
 allocate(data_one_tile(i_input,j_input))

 call ESMF_FieldGather(c_d_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'c_d ',data_one_tile
 endif
 call ESMF_FieldGather(c_0_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'c_0 ',data_one_tile
 endif
 call ESMF_FieldGather(d_conv_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'d_conv ',data_one_tile
 endif
 call ESMF_FieldGather(dt_cool_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'dt_cool ',data_one_tile
 endif
 call ESMF_FieldGather(ifd_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'ifd ',data_one_tile
 endif
 call ESMF_FieldGather(qrain_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'qrain ',data_one_tile
 endif
 call ESMF_FieldGather(tref_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'tref ',data_one_tile
 endif
 call ESMF_FieldGather(w_d_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'w_d ',data_one_tile
 endif
 call ESMF_FieldGather(w_0_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'w_0 ',data_one_tile
 endif
 call ESMF_FieldGather(xs_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xs ',data_one_tile
 endif
 call ESMF_FieldGather(xt_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xt ',data_one_tile
 endif
 call ESMF_FieldGather(xu_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xu ',data_one_tile
 endif
 call ESMF_FieldGather(xv_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xv ',data_one_tile
 endif
 call ESMF_FieldGather(xz_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xz ',data_one_tile
 endif
 call ESMF_FieldGather(xtts_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xtts ',data_one_tile
 endif
 call ESMF_FieldGather(xzts_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'xzts ',data_one_tile
 endif
 call ESMF_FieldGather(z_c_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'z_c ',data_one_tile
 endif
 call ESMF_FieldGather(zm_input_grid, data_one_tile, rootPet=0, rc=rc)
 print*,'rc ',rc
 if (localpet == 0) then
   print*,'zm ',data_one_tile
 endif

 deallocate(data_one_tile)

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI)

 call mpi_finalize(rc)

 print*,"SUCCESS!"

 end program readnst
