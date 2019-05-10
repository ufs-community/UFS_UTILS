 module program_setup

!--------------------------------------------------------------------------
! Module program_setup
!
! Abstract: Set up program execution
!
! Public Subroutines:
! -------------------
! read_setup_namelist          Reads configuration namelist
! calc_soil_params_driver      Computes soil parameters
!
! Public variables:
! -----------------
! atm_file_input_grid             File names of input atmospheric data.
!                                 History or gaussian input type only.
! atm_core_files_input_grid       File names of input atmospheric restart
!                                 core files.
! atm_tracer_files_input_grid     File names of input atmospheric restart
!                                 tracer files.
! atm_weight_file                 File containing pre-computed weights
!                                 to horizontally interpolate
!                                 atmospheric fields.
! bb_target                       Soil 'b' parameter, target grid
! convert_atm                     Convert atmospheric data when true.
! convert_nst                     Convert nst data when true.
! convert_sfc                     Convert sfc data when true.
! cres_target_grid                Target grid resolution, i.e., C768.
! cycle_mon/day/hour              Cycle month/day/hour
! data_dir_input_grid             Directory containing input atm or sfc
!                                 files.
! drysmc_input/target             Air dry soil moisture content input/
!                                 target grids.
! fix_dir_target_grid             Directory containing target grid
!                                 pre-computed fixed data (ex: soil type)
! halo_blend                      Number of row/cols of blending halo,
!                                 where model tendencies and lateral
!                                 boundary tendencies are applied.
!                                 Regional target grids only.
! halo_bndy                       Number of row/cols of lateral halo,
!                                 where pure lateral bndy conditions are 
!                                 applied (regional target grids).
! input_type                      Input data type: "restart" for fv3
!                                 tiled restart files; "history" for fv3
!                                 tiled history files; "gaussian"
!                                 for fv3 gaussian nemsio files;
!                                 "gfs_gaussian" for spectral gfs gaussian
!                                 nemsio files.
! max_tracers                     Maximum number of atmospheric tracers
!                                 processed
! maxsmc_input/target             Maximum soil moisture content input/
!                                 target grids
! mosaic_file_input_grid          Input grid mosaic file.  Not used
!                                 with "gaussian" or "gfs_gaussian"
!                                 input type.
! mosaic_file_target_grid         Target grid mosaic file
! nst_files_input_grid            File name of input nst data.  Only
!                                 used for input_type "gfs_gaussian".
! num_tracers                     Number of atmospheric tracers to
!                                 be processed.
! orog_dir_input_grid             Directory containing the input grid
!                                 orography files.  Not used for "gaussian"
!                                 or "gfs_gaussian" input types.
! orog_files_input_grid           Input grid orography files.  Not used
!                                 for "gaussian" or "gfs_gaussian"
!                                 input types.
! orog_dir_target_grid            Directory containing the target grid
!                                 orography files.
! orog_files_target_grid          Target grid orography files.
! refsmc_input/target             Reference soil moisture content input/
!                                 target grids (onset of soil moisture
!                                 stress).
! regional                        For regional target grids.  When '1'
!                                 remove boundary halo region from
!                                 atmospheric/surface data and
!                                 output atmospheric boundary file.
!                                 When '2' output boundary file only.
!                                 Default is '0' (global grids).
! satpsi_target                   Saturated soil potential, target grid
! sfc_files_input_grid            File names containing input surface data.
! tracers                         Name of each atmos tracer to be processed.
!                                 These names will be used to identify
!                                 the tracer records in the output files.
!                                 Follows the convention in the field table.
! tracers_input                   Name of each atmos tracer record in 
!                                 the input file.  May be different from
!                                 value in 'tracers'. 
! vcoord_file_target_grid         Vertical coordinate definition file
! wltsmc_input/target             Wilting point soil moisture content
!                                 input/target grids
! 
!--------------------------------------------------------------------------

 implicit none

 private
 
 character(len=500), public      :: base_install_dir = "NULL"
 character(len=500), public      :: atm_files_input_grid(6) = "NULL"
 character(len=500), public      :: atm_core_files_input_grid(7) = "NULL"
 character(len=500), public      :: atm_tracer_files_input_grid(6) = "NULL"
 character(len=500), public      :: data_dir_input_grid = "NULL"
 character(len=500), public      :: varmap_tables_dir = "parm/varmap_tables"
 character(len=500), public      :: fix_dir_target_grid = "NULL"
 character(len=500), public      :: mosaic_file_input_grid = "NULL"
 character(len=500), public      :: mosaic_file_target_grid = "NULL"
 character(len=500), public      :: nst_files_input_grid = "NULL"
 character(len=500), public      :: grib2_file_input_grid = "NULL"
 character(len=500), public      :: geogrid_file_input_grid="NULL"
 character(len=500), public      :: orog_dir_input_grid = "NULL"
 character(len=500), public      :: orog_files_input_grid(6) = "NULL"
 character(len=500), public      :: orog_dir_target_grid = "NULL"
 character(len=500), public      :: orog_files_target_grid(6) = "NULL"
 character(len=500), public      :: sfc_files_input_grid(6) = "NULL"
 character(len=500), public      :: vcoord_file_target_grid = "NULL"
 character(len=6),   public      :: cres_target_grid = "NULL"
 character(len=500), public      :: atm_weight_file="NULL"
 character(len=20),  public      :: input_type="restart"
 character(len=20),  public      :: external_model="GFS"  !Default assume gfs data
 character(len=20), public       :: phys_suite="GFS"      !Default to gfs physics suite
 character(len=1000), public     :: wgrib2_path="wgrib2"

 integer, parameter, public      :: max_tracers=100
 integer, public                 :: num_tracers
 
 logical, allocatable, public    :: read_from_input(:)
 
 character(len=20), public       :: tracers(max_tracers)="NULL"
 character(len=20), public       :: tracers_input(max_tracers)="NULL"
 character(len=20), allocatable, public      :: missing_var_methods(:)
 character(len=20), allocatable, public      :: chgres_var_names(:)
 character(len=20), allocatable, public      :: field_var_names(:)
 
 
 integer, public                 :: cycle_mon = -999
 integer, public                 :: cycle_day = -999
 integer, public                 :: cycle_hour = -999
 integer, public                 :: regional = 0
 integer, public                 :: halo_bndy = 4
 integer, public                 :: halo_blend = 0

 logical, public                 :: convert_atm = .false.
 logical, public                 :: convert_nst = .false.
 logical, public                 :: convert_sfc = .false.
 logical, public                 :: replace_o3mr = .true.
 logical, public                 :: replace_clwmr = .true.

 real, allocatable, public       :: drysmc_input(:), drysmc_target(:)
 real, allocatable, public       :: maxsmc_input(:), maxsmc_target(:)
 real, allocatable, public       :: refsmc_input(:), refsmc_target(:)
 real, allocatable, public       :: wltsmc_input(:), wltsmc_target(:)
 real, allocatable, public       :: bb_target(:),    satpsi_target(:)
 real, allocatable, public       :: missing_var_values(:)
 

 public :: read_setup_namelist
 public :: calc_soil_params_driver
 public :: read_varmap
 public :: get_var_cond

 contains

 subroutine read_setup_namelist

 implicit none

 integer                     :: is, ie, ierr


 namelist /config/ base_install_dir, &
                   varmap_tables_dir, &
                   mosaic_file_target_grid, &
                   fix_dir_target_grid,     &
                   orog_dir_target_grid,    &
                   orog_files_target_grid,  &
                   mosaic_file_input_grid,  &
                   orog_dir_input_grid,     &
                   orog_files_input_grid,   &
                   nst_files_input_grid,    &
                   sfc_files_input_grid,    &
                   atm_files_input_grid,    &
                   atm_core_files_input_grid,    &
                   atm_tracer_files_input_grid,    &
                   grib2_file_input_grid, &
                   geogrid_file_input_grid, &
                   data_dir_input_grid,     &
                   vcoord_file_target_grid, &
                   cycle_mon, cycle_day,    &
                   cycle_hour, convert_atm, &
                   convert_nst, convert_sfc, &
                   regional, input_type, external_model, &
                   atm_weight_file, tracers, &
                   tracers_input, replace_o3mr, &
                   replace_clwmr, phys_suite, wgrib2_path, &
                   halo_bndy, halo_blend

 print*,"- READ SETUP NAMELIST"

 open(41, file="./fort.41", iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(41, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
 close (41)
 
 orog_dir_target_grid = trim(orog_dir_target_grid) // '/'
 orog_dir_input_grid = trim(orog_dir_input_grid) // '/'

!-------------------------------------------------------------------------
! Determine CRES of target grid from the name of the mosaic file.
!-------------------------------------------------------------------------

 is = index(mosaic_file_target_grid, "/", .true.)
 ie = index(mosaic_file_target_grid, "_mosaic")

 if (is == 0 .or. ie == 0) then
   call error_handler("CANT DETERMINE CRES FROM MOSAIC FILE.", 1)
 endif
   
 cres_target_grid = mosaic_file_target_grid(is+1:ie-1)

 print *, convert_sfc, convert_atm
 if (.not. convert_sfc .and. .not. convert_atm) then
   call error_handler("MUST CONVERT EITHER AN ATM OR SFC FILE.", 1)
 endif

!-------------------------------------------------------------------------
! Flag for processing stand-alone regional grid.  When '1', 
! remove halo from atmospheric and surface data and output
! atmospheric lateral boundary condition file. When '2',
! create lateral boundary file only.  When '0' (the default),
! process normally as a global grid.
!-------------------------------------------------------------------------

 if (regional > 0) then
   !halo_bndy = 4
   !halo_blend = 0
   print*,"- PROCESSING A REGIONAL NEST WITH A BOUNDARY HALO OF ",halo_bndy
   print*,"- PROCESSING A REGIONAL NEST WITH A BLENDING HALO OF ",halo_blend
 endif

 num_tracers = 0
 do is = 1, max_tracers
   if (trim(tracers(is)) == "NULL") exit
   num_tracers = num_tracers + 1
   print*,"- WILL PROCESS TRACER ", trim(tracers(is))
 enddo

!-------------------------------------------------------------------------
! Ensure program recognizes the input data type.  
!-------------------------------------------------------------------------

 select case (trim(input_type))
   case ("restart")
     print*,'- INPUT DATA FROM FV3 TILED RESTART FILES.'
   case ("history")
     print*,'- INPUT DATA FROM FV3 TILED HISTORY FILES.'
   case ("gaussian")
     print*,'- INPUT DATA FROM FV3 GAUSSIAN NEMSIO FILE.'
   case ("gfs_gaussian")
     print*,'- INPUT DATA FROM SPECTRAL GFS GAUSSIAN NEMSIO FILE.'
   case ("gfs_spectral")
     print*,'- INPUT DATA FROM SPECTRAL GFS SIGIO/SFCIO FILE.'
   case ("grib2")
     print*,'- INPUT DATA FROM A GRIB2 FILE'
   case default
     call error_handler("UNRECOGNIZED INPUT DATA TYPE.", 1)
 end select

 return

 end subroutine read_setup_namelist

subroutine read_varmap

 implicit none

 integer                    :: istat, k, nvars
 character(len=500)         :: varmap_table_file
 
 varmap_table_file = trim(base_install_dir) // "/" // trim(varmap_tables_dir) // "/" &
                    // trim(phys_suite) // "phys_var_map.txt"
 

 print*,"OPEN VARIABLE MAPPING FILE: ", trim(varmap_table_file)
 open(14, file=trim(varmap_table_file), form='formatted', iostat=istat)
 if (istat /= 0) then
   call error_handler("OPENING VARIABLE MAPPING FILE", istat)
 endif

 read(14, *, iostat=istat) nvars
 if (istat /= 0) call error_handler("READING VARIABLE MAPPING FILE", istat)


 allocate(chgres_var_names(nvars))
 allocate(field_var_names(nvars))
 allocate(missing_var_methods(nvars))
 allocate(missing_var_values(nvars))
 allocate(read_from_input(nvars))

 
 read_from_input(:) = .true.

 
 do k = 1,nvars
  read(14, *, iostat=istat) chgres_var_names(k), field_var_names(k) , & 
                           missing_var_methods(k), missing_var_values(k)
 enddo
 
 if (istat /= 0) call error_handler("READING VARIABLE MAPPING FILE", istat)


 print*
 do k = 1, nvars
   print*,'VAR MISSING METHOD FOR', chgres_var_names(k), 'IS: ', missing_var_methods(k)
 enddo

 close(14)

 end subroutine read_varmap

! ----------------------------------------------------------------------------------------
! Find conditions for handling missing variables from varmap arrays
! ----------------------------------------------------------------------------------------

subroutine get_var_cond(var_name,this_miss_var_method,this_miss_var_value, &
                            this_field_var_name, loc)
  use esmf
  
  implicit none
  character(len=20)         :: var_name
  
  character(len=20), optional, intent(out) :: this_miss_var_method, &
                                              this_field_var_name
  real(esmf_kind_r4), optional, intent(out):: this_miss_var_value                                           
  
  integer, optional, intent(out)        :: loc
  
  integer                               :: i, tmp(size(missing_var_methods))
  
  i=0
  
  tmp(:)=0
  where(chgres_var_names==var_name) tmp=1
  
  i = maxloc(merge(1.,0.,chgres_var_names == var_name),dim=1) !findloc(chgres_var_names,var_name)
  print*, i
  if (maxval(tmp).eq.0) then
    print*, "WARNING: NO ENTRY FOR ", trim(var_name), " IN VARMAP TABLE. WILL SKIP " // &
            "VARIABLE IF NOT FOUND IN EXTERNAL MODEL FILE"
            
    if(present(this_miss_var_method)) this_miss_var_method = "skip"
    if(present(this_miss_var_value)) this_miss_var_value = -9999.9_esmf_kind_r4
    if(present(this_field_var_name)) this_field_var_name = "NULL"
    if(present(loc)) loc = 9999
  else
    if(present(this_miss_var_method)) this_miss_var_method = missing_var_methods(i)
    if(present(this_miss_var_value)) this_miss_var_value = missing_var_values(i)
    if(present(this_field_var_name)) this_field_var_name = field_var_names(i)
    if(present(loc)) loc = i
  endif
  
end subroutine get_var_cond

 subroutine calc_soil_params_driver(localpet)

 implicit none

 integer, intent(in)       :: localpet

 integer, parameter        :: num_statsgo = 16
 real, parameter           :: smlow_statsgo = 0.5
 real, parameter           :: smhigh_statsgo = 6.0

! zobler soil type used by spectral gfs prior to June 2017.
 integer, parameter        :: num_zobler = 9
 real, parameter           :: smlow_zobler = 0.5
 real, parameter           :: smhigh_zobler = 6.0

 integer                   :: num_soil_cats

 real                      :: bb_statsgo(num_statsgo)
 real                      :: maxsmc_statsgo(num_statsgo)
 real                      :: satdk_statsgo(num_statsgo)
 real                      :: satpsi_statsgo(num_statsgo)

 real                      :: bb_zobler(num_zobler)
 real                      :: maxsmc_zobler(num_zobler)
 real                      :: satdk_zobler(num_zobler)
 real                      :: satpsi_zobler(num_zobler)

 real, allocatable         :: bb(:)
 real                      :: smlow, smhigh
 real, allocatable         :: f11(:)
 real, allocatable         :: satdk(:)
 real, allocatable         :: satpsi(:)
 real, allocatable         :: satdw(:)

 data bb_statsgo /4.05, 4.26, 4.74, 5.33, 5.33, 5.25, &
            6.77, 8.72, 8.17, 10.73, 10.39, 11.55, &
            5.25, -9.99, 4.05, 4.26/

 data maxsmc_statsgo /0.395, 0.421, 0.434, 0.476, 0.476, 0.439, &
              0.404, 0.464, 0.465, 0.406, 0.468, 0.457, &
              0.464, -9.99, 0.200, 0.421/

 data satdk_statsgo /1.7600e-4, 1.4078e-5, 5.2304e-6, 2.8089e-6, 2.8089e-6, &
             3.3770e-6, 4.4518e-6, 2.0348e-6, 2.4464e-6, 7.2199e-6, &
             1.3444e-6, 9.7384e-7, 3.3770e-6,     -9.99, 1.4078e-5, &
             1.4078e-5/

 data satpsi_statsgo /0.0350, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548, &
              0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677, &
              0.3548, -9.99,  0.0350, 0.0363/

 data bb_zobler /4.26,  8.72, 11.55,  4.74, 10.73,  8.17, &
                 6.77,  5.25,  4.26/

 data maxsmc_zobler /0.421, 0.464, 0.468, 0.434, 0.406, 0.465, &
                     0.404, 0.439, 0.421/

 data satdk_zobler /1.41e-5, 0.20e-5, 0.10e-5, 0.52e-5, 0.72e-5, &
                    0.25e-5, 0.45e-5, 0.34e-5, 1.41e-5/

 data satpsi_zobler /0.040, 0.620, 0.470, 0.140, 0.100, 0.260,  &
                     0.140, 0.360, 0.040/

!-------------------------------------------------------------------------
! Compute soil parameters for the input grid.
!-------------------------------------------------------------------------

 select case (trim(input_type))
   case ("gfs_spectral")
     print*,'- INPUT GRID USED ZOBLER SOIL TYPES.'
     num_soil_cats = num_zobler
   case default
     print*,'- INPUT GRID USED STATSGO SOIL TYPES.'
     num_soil_cats = num_statsgo
 end select

 allocate(maxsmc_input(num_soil_cats))
 allocate(wltsmc_input(num_soil_cats))
 allocate(drysmc_input(num_soil_cats))
 allocate(refsmc_input(num_soil_cats))
 allocate(bb(num_soil_cats))
 allocate(satdk(num_soil_cats))
 allocate(satpsi(num_soil_cats))
 allocate(satdw(num_soil_cats))
 allocate(f11(num_soil_cats))

 select case (trim(input_type))
   case ("gfs_spectral")
     smlow  = smlow_zobler
     smhigh = smhigh_zobler
     maxsmc_input = maxsmc_zobler
     bb     = bb_zobler
     satdk  = satdk_zobler
     satpsi = satpsi_zobler
   case default
     smlow  = smlow_statsgo
     smhigh = smhigh_statsgo
     maxsmc_input = maxsmc_statsgo
     bb     = bb_statsgo
     satdk  = satdk_statsgo
     satpsi = satpsi_statsgo
 end select

 call calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc_input, &
                       bb, satpsi, satdw, f11, refsmc_input, drysmc_input, wltsmc_input)

 deallocate(bb, satdk, satpsi, satdw, f11)

 if (localpet == 0) print*,'maxsmc input grid ',maxsmc_input
 if (localpet == 0) print*,'wltsmc input grid ',wltsmc_input

!-------------------------------------------------------------------------
! Compute soil parameters for the target grid.
!-------------------------------------------------------------------------

 print*,'- TARGET GRID USEING STATSGO SOIL TYPES.'

 num_soil_cats = num_statsgo

 allocate(maxsmc_target(num_soil_cats))
 allocate(wltsmc_target(num_soil_cats))
 allocate(drysmc_target(num_soil_cats))
 allocate(refsmc_target(num_soil_cats))
 allocate(bb_target(num_soil_cats))
 allocate(satpsi_target(num_soil_cats))
 allocate(satdk(num_soil_cats))
 allocate(satdw(num_soil_cats))
 allocate(f11(num_soil_cats))

 smlow  = smlow_statsgo
 smhigh = smhigh_statsgo
 maxsmc_target = maxsmc_statsgo
 bb_target     = bb_statsgo
 satdk  = satdk_statsgo
 satpsi_target = satpsi_statsgo

 call calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc_target, &
                       bb_target, satpsi_target, satdw, f11, refsmc_target, drysmc_target, wltsmc_target)

 deallocate(satdk, satdw, f11)

 if (localpet == 0) print*,'maxsmc target grid ',maxsmc_target
 if (localpet == 0) print*,'wltsmc input grid ',wltsmc_target

 end subroutine calc_soil_params_driver

 subroutine calc_soil_params(num_soil_cats, smlow, smhigh, satdk,  &
            maxsmc, bb, satpsi, satdw, f11, refsmc, drysmc, wltsmc)

 implicit none

 integer, intent(in)            :: num_soil_cats

 real, intent(in)               :: smlow, smhigh
 real, intent(in)               :: bb(num_soil_cats)
 real, intent(in)               :: maxsmc(num_soil_cats)
 real, intent(in)               :: satdk(num_soil_cats)
 real, intent(in)               :: satpsi(num_soil_cats)

 real, intent(out)              :: f11(num_soil_cats)
 real, intent(out)              :: satdw(num_soil_cats)
 real, intent(out)              :: refsmc(num_soil_cats)
 real, intent(out)              :: drysmc(num_soil_cats)
 real, intent(out)              :: wltsmc(num_soil_cats)

 integer                        :: i

 real                      :: refsmc1
 real                      :: wltsmc1

 satdw = 0.0
 f11   = 0.0
 refsmc = 0.0
 wltsmc = 0.0
 drysmc = 0.0

 do i = 1, num_soil_cats

   if (maxsmc(i) > 0.0) then

   SATDW(I)  = BB(I)*SATDK(I)*(SATPSI(I)/MAXSMC(I))
   F11(I) = ALOG10(SATPSI(I)) + BB(I)*ALOG10(MAXSMC(I)) + 2.0
   REFSMC1 = MAXSMC(I)*(5.79E-9/SATDK(I)) **(1.0/(2.0*BB(I)+3.0))
   REFSMC(I) = REFSMC1 + (MAXSMC(I)-REFSMC1) / SMHIGH
   WLTSMC1 = MAXSMC(I) * (200.0/SATPSI(I))**(-1.0/BB(I))
   WLTSMC(I) = WLTSMC1 - SMLOW * WLTSMC1

!----------------------------------------------------------------------
!  CURRENT VERSION DRYSMC VALUES THAT EQUATE TO WLTSMC.
!  FUTURE VERSION COULD LET DRYSMC BE INDEPENDENTLY SET VIA NAMELIST.
!----------------------------------------------------------------------

   DRYSMC(I) = WLTSMC(I)

   end if

 END DO

 end subroutine calc_soil_params

 end module program_setup
