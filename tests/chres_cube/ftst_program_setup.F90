! Unit test for chres_cube utility, input_data.F90.
!
! Ed Hartnett 2/16/21

program ftst_program_setup
  use esmf
  use netcdf
  use program_setup
  implicit none
  
  print*, "Starting test of program_setup."

  print*, "testing read_setup_namelist..."
  call read_setup_namelist()
  if (cycle_mon .ne. 7 .or. cycle_day .ne. 4 .or. cycle_hour .ne. 12) stop 4
  if (.not. convert_atm .or. .not. convert_sfc .or. .not. convert_nst) stop 5
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 6
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/C96_mosaic.nc") stop 7
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/fix_sfc") stop 8
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/") stop 9
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l64.txt") stop 10
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/fv3.nemsio") stop 11
  if (trim(atm_files_input_grid(1)) .ne. 'gfs.t12z.atmf000.nemsio') stop 12
  if (trim(sfc_files_input_grid(1)) .ne. 'gfs.t12z.sfcf000.nemsio') stop 13
  if (varmap_file .ne. "NULL") stop 14
  if (thomp_mp_climo_file .ne. "NULL") stop 16
  if (trim(cres_target_grid) .ne. "C96") stop 17
  if (atm_weight_file .ne. "NULL") stop 18
  if (trim(input_type) .ne. "gaussian_nemsio") stop 19
  if (trim(external_model) .ne. "GFS") stop 20
  if (num_tracers .ne. 7) stop 21
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "liq_wat" .or. tracers(3) .ne. "o3mr" .or. &
       tracers(4) .ne. "ice_wat" .or. tracers(5) .ne. "rainwat" .or. tracers(6) .ne. "snowwat" .or. &
       tracers(7) .ne. "graupel") stop 22
  if (tracers_input(1) .ne. "spfh" .or. tracers_input(2) .ne. "clwmr" .or. &
       tracers_input(3) .ne. "o3mr" .or. tracers_input(4) .ne. "icmr" .or. &
       tracers_input(5) .ne. "rwmr" .or. tracers_input(6) .ne. "snmr" .or. &
       tracers_input(7) .ne. "grle") stop 23
  print*, "OK"
  
  print*, "SUCCESS!"
end program ftst_program_setup
