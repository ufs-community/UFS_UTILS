program ftst_program_setup

! Unit test for emc_snow2mdl utility, program_setup.
!
! Reads the program namelist and compares each
! variable to expected values.
!
! Author: George Gayno

 use program_setup

 implicit none

 print*, "Starting test of program_setup."
 print*, "testing read_setup_namelist with file fort.41..."

 call read_config_nml

 if (trim(autosnow_file) /= "autosnow.grb") stop 2
 if (trim(nesdis_snow_file) /= "imssnow96.grb") stop 3
 if (trim(nesdis_lsmask_file) /= "mask.grb") stop 4
 if (trim(afwa_snow_global_file) /= "global_snow.grb") stop 5
 if (trim(afwa_snow_nh_file) /= "NPR.SNWN.SP.S1200.MESH16") stop 6
 if (trim(afwa_snow_sh_file) /= "NPR.SNWS.SP.S1200.MESH16") stop 7
 if (trim(afwa_lsmask_nh_file) /= "afwa_mask.nh.bin") stop 8
 if (trim(afwa_lsmask_sh_file) /= "afwa_mask.sh.bin") stop 9

 if (trim(climo_qc_file) /= "emcsfc_snow_cover_climo.grib2") stop 10

 if (trim(model_lat_file) /= "global_latitudes.t1534.3072.1536.grb") stop 11
 if (trim(model_lon_file) /= "global_longitudes.t1534.3072.1536.grb") stop 12
 if (trim(model_lsmask_file) /= "global_slmask.t1534.3072.1536.grb") stop 13
 if (trim(gfs_lpl_file) /= "global_lonsperlat.t1534.3072.1536.txt") stop 14

 if (trim(model_snow_file) /= "snogrb_model") stop 15
 if (output_grib2) stop 16

 if (grib_year /= 12) stop 17
 if (grib_month /= 10) stop 18
 if (grib_day /= 29) stop 19
 if (grib_hour /= 0) stop 20

 if (lat_threshold /= 55.0) stop 22
 if (min_snow_depth /= 0.05) stop 23
 if (snow_cvr_threshold /= 50.0) stop 24

 print*, "OK"
 print*, "SUCCESS!"

end program ftst_program_setup
