! Unit test for chres_cube utility, input_data.F90.
!
! Ed Hartnett 2/16/21

program ftst_program_setup
  use mpi
  use esmf
  use netcdf
  use program_setup
  implicit none
  integer :: is
  integer :: my_rank, nprocs
  integer :: ierr

  call mpi_init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (my_rank .eq. 0) print*, "Starting test of program_setup."
  if (my_rank .eq. 0) print*, "testing read_setup_namelist with file fort.41..."
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

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "testing read_setup_namelist with config_fv3_tiled..."
  call read_setup_namelist("config_fv3_tiled.nml")
  if (cycle_mon .ne. 10 .or. cycle_day .ne. 3 .or. cycle_hour .ne. 0) stop 34
  if (.not. convert_atm .or. .not. convert_sfc .or. .not. convert_nst) stop 35
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 36
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/C192_mosaic.nc") stop 37
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/fix_sfc") stop 38
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/") stop 39
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l64.txt") stop 40
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/fv3.history") stop 41
  if (trim(atm_files_input_grid(1)) .ne. 'dynf000.tile1.nc') stop 42
  if (trim(sfc_files_input_grid(1)) .ne. 'phyf000.tile1.nc') stop 43
  if (varmap_file .ne. "NULL") stop 44
  if (thomp_mp_climo_file .ne. "NULL") stop 46
  if (trim(cres_target_grid) .ne. "C192") stop 47
  if (atm_weight_file .ne. "NULL") stop 48
  if (trim(input_type) .ne. "history") stop 49
  if (trim(external_model) .ne. "GFS") stop 50
  if (num_tracers .ne. 3) stop 51
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "liq_wat" .or. tracers(3) .ne. "o3mr") stop 52
  if (tracers_input(1) .ne. "spfh" .or. tracers_input(2) .ne. "clwmr" .or. &
       tracers_input(3) .ne. "o3mr") stop 53

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  if (my_rank .eq. 0) print*, "OK"

  print*, "testing read_setup_namelist with config_fv3_tiled_warm_restart..."
  call read_setup_namelist("config_fv3_tiled_warm_restart.nml")
  if (cycle_mon .ne. 7 .or. cycle_day .ne. 6 .or. cycle_hour .ne. 12) stop 64
  if (.not. convert_atm .or. .not. convert_sfc .or. .not. convert_nst) stop 65
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 66  
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/C96_mosaic.nc") stop 67
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/fix_sfc") stop 68
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/") stop 69
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l64.txt") stop 70
  if (trim(mosaic_file_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C384/C384_mosaic.nc") stop 71
  if (trim(orog_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C384/") stop 72
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/fv3.restart") stop 73
  if (trim(orog_files_target_grid(1)) .ne. "C96_oro_data.tile1.nc") stop 176
  if (trim(orog_files_target_grid(2)) .ne. "C96_oro_data.tile2.nc") stop 176
  if (trim(orog_files_target_grid(3)) .ne. "C96_oro_data.tile3.nc") stop 176
  if (trim(orog_files_target_grid(4)) .ne. "C96_oro_data.tile4.nc") stop 176
  if (trim(orog_files_target_grid(5)) .ne. "C96_oro_data.tile5.nc") stop 176
  if (trim(orog_files_target_grid(6)) .ne. "C96_oro_data.tile6.nc") stop 176
  if (trim(orog_files_input_grid(1)) .ne. "C384_oro_data.tile1.nc") stop 177
  if (trim(orog_files_input_grid(2)) .ne. "C384_oro_data.tile2.nc") stop 177
  if (trim(orog_files_input_grid(3)) .ne. "C384_oro_data.tile3.nc") stop 177
  if (trim(orog_files_input_grid(4)) .ne. "C384_oro_data.tile4.nc") stop 177
  if (trim(orog_files_input_grid(5)) .ne. "C384_oro_data.tile5.nc") stop 177
  if (trim(orog_files_input_grid(6)) .ne. "C384_oro_data.tile6.nc") stop 177
  if (trim(atm_core_files_input_grid(1)) .ne. "20190706.120000.fv_core.res.tile1.nc") stop 178
  if (trim(atm_core_files_input_grid(2)) .ne. "20190706.120000.fv_core.res.tile2.nc") stop 178
  if (trim(atm_core_files_input_grid(3)) .ne. "20190706.120000.fv_core.res.tile3.nc") stop 178
  if (trim(atm_core_files_input_grid(4)) .ne. "20190706.120000.fv_core.res.tile4.nc") stop 178
  if (trim(atm_core_files_input_grid(5)) .ne. "20190706.120000.fv_core.res.tile5.nc") stop 178
  if (trim(atm_core_files_input_grid(6)) .ne. "20190706.120000.fv_core.res.tile6.nc") stop 178
  if (trim(atm_tracer_files_input_grid(1)) .ne. "20190706.120000.fv_tracer.res.tile1.nc") stop 179
  if (trim(atm_tracer_files_input_grid(2)) .ne. "20190706.120000.fv_tracer.res.tile2.nc") stop 179
  if (trim(atm_tracer_files_input_grid(3)) .ne. "20190706.120000.fv_tracer.res.tile3.nc") stop 179
  if (trim(atm_tracer_files_input_grid(4)) .ne. "20190706.120000.fv_tracer.res.tile4.nc") stop 179
  if (trim(atm_tracer_files_input_grid(5)) .ne. "20190706.120000.fv_tracer.res.tile5.nc") stop 179
  if (trim(atm_tracer_files_input_grid(6)) .ne. "20190706.120000.fv_tracer.res.tile6.nc") stop 179
  if (trim(sfc_files_input_grid(1)) .ne. "20190706.120000.sfc_data.tile1.nc") stop 180
  if (trim(sfc_files_input_grid(2)) .ne. "20190706.120000.sfc_data.tile2.nc") stop 180
  if (trim(sfc_files_input_grid(3)) .ne. "20190706.120000.sfc_data.tile3.nc") stop 180
  if (trim(sfc_files_input_grid(4)) .ne. "20190706.120000.sfc_data.tile4.nc") stop 180
  if (trim(sfc_files_input_grid(5)) .ne. "20190706.120000.sfc_data.tile5.nc") stop 180
  if (trim(sfc_files_input_grid(6)) .ne. "20190706.120000.sfc_data.tile6.nc") stop 180
  if (trim(input_type) .ne. "restart") stop 89
  if (num_tracers .ne. 7) stop 173
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "liq_wat" .or. tracers(3) .ne. "o3mr" .or. &
       tracers(4) .ne. "ice_wat" .or. tracers(5) .ne. "rainwat" .or. tracers(6) .ne. "snowwat" .or. &
       tracers(7) .ne. "graupel") stop 174
  if (tracers_input(1) .ne. "sphum" .or. tracers_input(2) .ne. "liq_wat" .or. &
       tracers_input(3) .ne. "o3mr" .or. tracers_input(4) .ne. "ice_wat" .or. &
       tracers_input(5) .ne. "rainwat" .or. tracers_input(6) .ne. "snowwat" .or. &
       tracers_input(7) .ne. "graupel") stop 175

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  print*, "OK"

  if (my_rank .eq. 0) print*, "testing read_setup_namelist with config_gaussian_nemsio..."
  call read_setup_namelist("config_gaussian_nemsio.nml")
  if (cycle_mon .ne. 7 .or. cycle_day .ne. 4 .or. cycle_hour .ne. 12) stop 74
  if (.not. convert_atm .or. .not. convert_sfc .or. .not. convert_nst) stop 75
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 76
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/C96_mosaic.nc") stop 77
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/fix_sfc") stop 78
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/") stop 79
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l64.txt") stop 80
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/fv3.nemsio") stop 81
  if (trim(atm_files_input_grid(1)) .ne. 'gfs.t12z.atmf000.nemsio') stop 82
  if (trim(sfc_files_input_grid(1)) .ne. 'gfs.t12z.sfcf000.nemsio') stop 83
  if (varmap_file .ne. "NULL") stop 84
  if (thomp_mp_climo_file .ne. "NULL") stop 86
  if (trim(cres_target_grid) .ne. "C96") stop 87
  if (atm_weight_file .ne. "NULL") stop 88
  if (trim(input_type) .ne. "gaussian_nemsio") stop 89
  if (trim(external_model) .ne. "GFS") stop 90
  if (num_tracers .ne. 7) stop 21
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "liq_wat" .or. tracers(3) .ne. "o3mr" .or. &
       tracers(4) .ne. "ice_wat" .or. tracers(5) .ne. "rainwat" .or. tracers(6) .ne. "snowwat" .or. &
       tracers(7) .ne. "graupel") stop 22
  if (tracers_input(1) .ne. "spfh" .or. tracers_input(2) .ne. "clwmr" .or. &
       tracers_input(3) .ne. "o3mr" .or. tracers_input(4) .ne. "icmr" .or. &
       tracers_input(5) .ne. "rwmr" .or. tracers_input(6) .ne. "snmr" .or. &
       tracers_input(7) .ne. "grle") stop 23

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "testing read_setup_namelist with config_spectral_sigio..."
  call read_setup_namelist("config_spectral_sigio.nml")
  if (cycle_mon .ne. 7 .or. cycle_day .ne. 17 .or. cycle_hour .ne. 0) stop 114
  if (.not. convert_atm .or. .not. convert_sfc .or. convert_nst) stop 115
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 116
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/C96_mosaic.nc") stop 117
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/fix_sfc") stop 118
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C96/") stop 119
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l64.txt") stop 120
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/gfs.sigio") stop 121
  if (trim(atm_files_input_grid(1)) .ne. 'gdas.t00z.sanl') stop 122
  if (trim(sfc_files_input_grid(1)) .ne. 'gdas.t00z.sfcanl') stop 123
  if (varmap_file .ne. "NULL") stop 124
  if (thomp_mp_climo_file .ne. "NULL") stop 126
  if (trim(cres_target_grid) .ne. "C96") stop 127
  if (atm_weight_file .ne. "NULL") stop 128
  if (trim(input_type) .ne. "gfs_sigio") stop 129
  if (trim(external_model) .ne. "GFS") stop 130
  if (num_tracers .ne. 3) stop 131
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "o3mr" .or. tracers(3) .ne. "liq_wat") stop 132
  if (tracers_input(1) .ne. "spfh" .or. tracers_input(2) .ne. "o3mr" .or. &
       tracers_input(3) .ne. "clwmr") stop 133

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "testing read_setup_namelist with config_gfs_grib2..."
  call read_setup_namelist("config_gfs_grib2.nml")
  if (cycle_mon .ne. 11 .or. cycle_day .ne. 4 .or. cycle_hour .ne. 0) stop 94
  if (.not. convert_atm .or. .not. convert_sfc .or. convert_nst) stop 95
  if (regional .ne. 0 .or. halo_bndy .ne. 0 .or. halo_blend .ne. 0) stop 96
  if (trim(mosaic_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/C192_mosaic.nc") stop 97
  if (trim(fix_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/fix_sfc") stop 98
  if (trim(orog_dir_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/fix/C192/") stop 99
  if (trim(vcoord_file_target_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../fix/fix_am/global_hyblev.l65.txt") stop 100
  if (trim(data_dir_input_grid) .ne. "/scratch1/NCEPDEV/da/George.Gayno/noscrub/reg_tests/chgres_cube/input_data/gfs.grib2") stop 101
  if (trim(grib2_file_input_grid) .ne. 'gfs.t00z.pgrb2.0p50.f000') stop 102
  if (trim(sfc_files_input_grid(1)) .ne. 'gdas.t00z.sfcanl') stop 103
  if (varmap_file .ne. "/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/reg_tests/chgres_cube/../../parm/varmap_tables/GFSphys_var_map.txt") stop 1010
  if (thomp_mp_climo_file .ne. "NULL") stop 106
  if (trim(cres_target_grid) .ne. "C192") stop 107
  if (atm_weight_file .ne. "NULL") stop 108
  if (trim(input_type) .ne. "grib2") stop 109
  if (trim(external_model) .ne. "GFS") stop 110
  if (num_tracers .ne. 7) stop 111
  if (tracers(1) .ne. "sphum" .or. tracers(2) .ne. "liq_wat" .or. tracers(3) .ne. "o3mr" .or. &
       tracers(4) .ne. "ice_wat" .or. tracers(5) .ne. "rainwat" .or. tracers(6) .ne. "snowwat" .or. &
       tracers(7) .ne. "graupel") stop 22
  if (tracers_input(1) .ne. "spfh" .or. tracers_input(2) .ne. "clwmr" .or. &
       tracers_input(3) .ne. "o3mr" .or. tracers_input(4) .ne. "icmr" .or. &
       tracers_input(5) .ne. "rwmr" .or. tracers_input(6) .ne. "snmr" .or. &
       tracers_input(7) .ne. "grle") stop 23

  ! Reset the tracers array.
  do is = 1, max_tracers
     tracers(is) = "NULL"
     tracers_input(is) = "NULL"
  enddo
  if (my_rank .eq. 0) print*, "OK"

  if (my_rank .eq. 0) print*, "SUCCESS!"

 call mpi_finalize(ierr)  
end program ftst_program_setup
