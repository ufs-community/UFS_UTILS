!> @file
!! @brief Set up program execution.
!! @author George Gayno NCEP/EMC

!> This module contains code to read the setup namelist file, handle
!! the varmap file for GRIB2 data, and calculate the soil parameters.
!!
!! @author George Gayno NCEP/EMC
 module program_setup

 use utilities, only                    : error_handler, to_lower

 implicit none

 private
 
 character(len=500), public      :: varmap_file = "NULL" !< REQUIRED. Full path of the relevant varmap file.
 character(len=500), public      :: atm_files_input_grid(6) = "NULL" !< File names of input
                                                                     !< atmospheric data. Not used
                                                                     !< for "grib2" or "restart"
                                                                     !< input types.
 character(len=500), public      :: atm_core_files_input_grid(7) = "NULL" !<  File names of input atmospheric restart core files.  Only used for 'restart' input type.
 character(len=500), public      :: atm_tracer_files_input_grid(6) = "NULL" !< File names of input atmospheric restart tracer files.  Only used for 'restart' input type.
 character(len=500), public      :: data_dir_input_grid = "NULL"  !< Directory containing input atm or sfc files.
 character(len=500), public      :: fix_dir_target_grid = "NULL" !< Directory containing target grid pre-computed fixed data (ex: soil type).
 character(len=500), public      :: mosaic_file_input_grid = "NULL" !< Input grid mosaic file.  Only used for "restart" or "history" input type.
 character(len=500), public      :: mosaic_file_target_grid = "NULL" !< Target grid mosaic file.
 character(len=500), public      :: nst_files_input_grid = "NULL" !< File name of input nst data.  Only used for input_type "gfs_gaussian_nemsio".
 character(len=500), public      :: grib2_file_input_grid = "NULL" !<  REQUIRED. File name of grib2 input data. Assumes atmospheric and surface data are in a single file. 
 character(len=500), public      :: geogrid_file_input_grid = "NULL" !< Name of "geogrid" file, which contains static
                                                                     !! surface fields on the input grid.  GRIB2 option
                                                                     !! only.
 character(len=500), public      :: orog_dir_input_grid = "NULL" !<  Directory containing the input grid orography files.  Only used for "restart" or "history" input types.
 character(len=500), public      :: orog_files_input_grid(6) = "NULL" !<  Input grid orography files.  Only used for "restart" or "history" input types.
 character(len=500), public      :: orog_dir_target_grid = "NULL" !<  Directory containing the target grid orography files.
 character(len=500), public      :: orog_files_target_grid(6) = "NULL" !<  Target grid orography files.
 character(len=500), public      :: sfc_files_input_grid(6) = "NULL" !<  File names containing input surface data. Not used for 'grib2' input type.
 character(len=500), public      :: vcoord_file_target_grid = "NULL" !<  Vertical coordinate definition file.
 character(len=500), public      :: thomp_mp_climo_file= "NULL" !<  Path/name to the Thompson MP climatology file.
 character(len=6),   public      :: cres_target_grid = "NULL" !<  Target grid resolution, i.e., C768.
 character(len=500), public      :: atm_weight_file="NULL" !<  File containing pre-computed weights to horizontally interpolate atmospheric fields.
 character(len=25),  public      :: input_type="restart" !< Input data type: 
!!                                 - "restart" for fv3 tiled warm restart
!!                                    files (netcdf).
!!                                 - "history" for fv3 tiled history files
!!                                    (netcdf).
!!                                 - "gaussian_nemsio" for fv3 gaussian
!!                                    nemsio files;
!!                                 - "gaussian_netcdf" for fv3 gaussian
!!                                    netcdf files.
!!                                 - "grib2" for grib2 files.
!!                                 - "gfs_gaussian_nemsio" for spectral gfs
!!                                    gaussian nemsio files
!!                                 - "gfs_sigio" for spectral gfs
!!                                    gfs sigio/sfcio files.
 character(len=20),  public      :: external_model="GFS"  !< The model that the input data is derived from. Current supported options are: "GFS", "HRRR", "NAM", "RAP". Default: "GFS"
 
 integer, parameter, public      :: max_tracers=100 !< Maximum number of atmospheric tracers processed.
 integer, public                 :: num_tracers !< Number of atmospheric tracers to be processed.
 integer, public                 :: num_tracers_input !< Number of atmospheric tracers in input file.
 
 logical, allocatable, public    :: read_from_input(:) !< When false, variable was not read from GRIB2 
                                                       !! input file.
 
 character(len=20), public       :: tracers(max_tracers)="NULL" !< Name of each atmos tracer to be processed.
                                                                !! These names will be used to identify
                                                                !! the tracer records in the output files.
                                                                !! Follows the convention in the field table.
                                                                !! FOR GRIB2 FILES: Not used. Tracers instead taken
                                                                !! from the varmap file. 
 character(len=20), public       :: tracers_input(max_tracers)="NULL" !<  Name of each atmos tracer record in 
                                                                      !! the input file.  May be different from
                                                                      !! value in 'tracers'. 
                                                                      !! FOR GRIB2 FILES: Not used. Tracers instead taken
                                                                      !! from the varmap file. 
 character(len=20), allocatable, public      :: missing_var_methods(:) !< Method to replace missing GRIB2 input
                                                                       !! records.
 character(len=20), allocatable, public      :: chgres_var_names(:) !< Varmap table variable name as recognized
                                                                    !! by this program.
 character(len=20), allocatable, public      :: field_var_names(:)  !< The GRIB2 variable name in the varmap table.
 
 
 integer, public                 :: cycle_year = -999 !< Cycle year.
 integer, public                 :: cycle_mon = -999 !< Cycle month.
 integer, public                 :: cycle_day = -999 !< Cycle day.
 integer, public                 :: cycle_hour = -999 !< Cycle hour.
 integer, public                 :: regional = 0 !<  For regional target grids.  When '1' remove boundary halo region from atmospheric/surface data and
                                                 !! output atmospheric boundary file. When '2' output boundary file only. Default is '0' (global grids).
 integer, public                 :: halo_bndy = 0 !< Number of row/cols of lateral halo, where pure lateral bndy conditions are applied (regional target grids).
 integer, public                 :: halo_blend = 0 !< Number of row/cols of blending halo, where model tendencies and lateral boundary tendencies are applied. Regional target grids only.
 integer, public                 :: nsoill_out = 4 !<  Number of soil levels desired in the output data. chgres_cube can interpolate from 9 input to 4 output levels. DEFAULT: 4.

 logical, public                 :: convert_atm = .false. !< Convert atmospheric data when true.
 logical, public                 :: convert_nst = .false. !< Convert nst data when true.
 logical, public                 :: convert_sfc = .false. !< Convert sfc data when true.
 logical, public                 :: wam_cold_start = .false. !< When true, cold start for whole atmosphere model.
 
 ! Options for replacing vegetation/soil type, veg fraction, and lai with data from the grib2 file
 ! Default is to use climatology instead
 logical, public                 :: vgtyp_from_climo = .true. !<  If false, interpolate vegetation type from the input 
                                                              !! data to the target grid instead of using data from 
                                                              !! static data. Use with caution as vegetation categories
                                                              !! can vary. Default: True.
 logical, public                 :: sotyp_from_climo = .true. !<  If false, interpolate soil type from the input 
                                                              !! data to the target grid instead of using data from 
                                                              !! static data. Use with caution as the code assumes
                                                              !! input soil type use STATSGO soil categories.
                                                              !! Default: True.
 logical, public                 :: vgfrc_from_climo = .true. !<  If false, interpolate vegetation fraction from the input 
                                                              !! data to the target grid instead of using data from 
                                                              !! static data. Use with caution as vegetation categories
                                                              !! can vary. Default: True.

 logical, public                 :: minmax_vgfrc_from_climo = .true. !<  If false, interpolate min/max vegetation fraction from 
                                                                     !! the input data to the target grid instead of using data
                                                                     !! from static data. Use with caution as vegetation
                                                                     !! categories can vary. Default: True.
 logical, public                 :: lai_from_climo = .true. !< If false, interpolate leaf area index from the input 
                                                            !! data to the target grid instead of using data from 
                                                            !! static data. Default: True.
 logical, public                 :: tg3_from_soil = .false. !<  If false, use lowest level soil temperature for the
                                                            !! base soil temperature instead of using data from 
                                                            !! static data. Default: False.
 logical, public                 :: use_thomp_mp_climo=.false. !<  When true, read and process Thompson MP climatological tracers.  False, when 'thomp_mp_climo_file' is NULL.

 real, allocatable, public       :: drysmc_input(:)   !< Air dry soil moisture content input grid.
 real, allocatable, public       :: drysmc_target(:)  !< Air dry soil moisture content target grid.
 real, allocatable, public       :: maxsmc_input(:) !< Maximum soil moisture content input grid.
 real, allocatable, public       :: maxsmc_target(:) !< Maximum soil moisture content target grid.
 real, allocatable, public       :: refsmc_input(:) !<  Reference soil moisture content input grid (onset of soil moisture stress).
 real, allocatable, public       :: refsmc_target(:) !<  Reference soil moisture content target grid (onset of soil moisture stress).
 real, allocatable, public       :: wltsmc_input(:) !< Plant wilting point soil moisture content input grid.
 real, allocatable, public       :: wltsmc_target(:) !< Plant wilting point soil moisture content target grid.
 real, allocatable, public       :: bb_target(:)  !<  Soil 'b' parameter, target grid
 real, allocatable, public       :: satpsi_target(:) !<   Saturated soil potential, target grid
 real, allocatable, public       :: missing_var_values(:) !< If input GRIB2 record is missing, the variable
                                                          !! is set to this value.
 

 public :: read_setup_namelist
 public :: calc_soil_params_driver
 public :: read_varmap
 public :: get_var_cond

 contains

!> Reads program configuration namelist.
!!
!! @param filename the name of the configuration file (defaults to
!! ./fort.41).
!! @author George Gayno NCEP/EMC
 subroutine read_setup_namelist(filename)
 implicit none

 character(len=*), intent(in), optional :: filename
 character(:), allocatable :: filename_to_use
 

 integer                     :: is, ie, ierr


 namelist /config/ varmap_file, &
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
                   cycle_year, cycle_mon, cycle_day,    &
                   cycle_hour, convert_atm, &
                   convert_nst, convert_sfc, &
                   wam_cold_start, &
                   vgtyp_from_climo, &
                   sotyp_from_climo, &
                   vgfrc_from_climo, &
                   minmax_vgfrc_from_climo, &
                   lai_from_climo, tg3_from_soil, &
                   regional, input_type, &
                   external_model, &
                   atm_weight_file, tracers, &
                   tracers_input, &
                   halo_bndy, & 
                   halo_blend, &
                   nsoill_out, &
                   thomp_mp_climo_file

 print*,"- READ SETUP NAMELIST"

 if (present(filename)) then
    filename_to_use = filename
 else
    filename_to_use = "./fort.41"
 endif

 open(41, file=filename_to_use, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(41, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
 close (41)
 
 call to_lower(input_type)
 
 orog_dir_target_grid = trim(orog_dir_target_grid) // '/'
 orog_dir_input_grid = trim(orog_dir_input_grid) // '/'

!-------------------------------------------------------------------------
! Determine CRES of target grid from the name of the mosaic file.
!-------------------------------------------------------------------------

 is = index(mosaic_file_target_grid, "/", .true.)
 ie = index(mosaic_file_target_grid, "mosaic") - 1

 if (is == 0 .or. ie == 0) then
   call error_handler("CANT DETERMINE CRES FROM MOSAIC FILE.", 1)
 endif
   
 cres_target_grid = mosaic_file_target_grid(is+1:ie-1)

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
   print*,"- PROCESSING A REGIONAL NEST WITH A BOUNDARY HALO OF ",halo_bndy
   print*,"- PROCESSING A REGIONAL NEST WITH A BLENDING HALO OF ",halo_blend
 else
   halo_bndy = 0
   halo_blend = 0
 endif

 num_tracers = 0
 do is = 1, max_tracers
   if (trim(tracers(is)) == "NULL") exit
   num_tracers = num_tracers + 1
   print*,"- TRACER NAME IN OUTPUT FILE ", trim(tracers(is))
 enddo
 
 num_tracers_input = 0
 do is = 1, max_tracers
   if (trim(tracers_input(is)) == "NULL") exit
   num_tracers_input = num_tracers_input + 1
   print*,"- TRACER NAME IN INPUT FILE ", trim(tracers_input(is))
 enddo

!-------------------------------------------------------------------------
! Ensure spo, spo2, and spo3 in tracers list if wam ic is on
!-------------------------------------------------------------------------

 if( wam_cold_start ) then
    ierr=3
    do is = 1, num_tracers
      if(trim(tracers(is)) == "spo"  ) ierr = ierr - 1
      if(trim(tracers(is)) == "spo2" ) ierr = ierr - 1
      if(trim(tracers(is)) == "spo3" ) ierr = ierr - 1
    enddo
    if (ierr /= 0) then
      print*,"-ERROR: spo, spo2, and spo3 should be in tracers namelist"
      call error_handler("WAM TRACER NAMELIST.", ierr)
    endif
    print*,"- WAM COLDSTART OPTION IS TURNED ON."
 endif  

!-------------------------------------------------------------------------
! Ensure program recognizes the input data type.  
!-------------------------------------------------------------------------

 select case (trim(input_type))
   case ("restart")
     print*,'- INPUT DATA FROM FV3 TILED RESTART FILES.'
   case ("history")
     print*,'- INPUT DATA FROM FV3 TILED HISTORY FILES.'
   case ("gaussian_nemsio")
     print*,'- INPUT DATA FROM FV3 GAUSSIAN NEMSIO FILE.'
   case ("gfs_gaussian_nemsio")
     print*,'- INPUT DATA FROM SPECTRAL GFS GAUSSIAN NEMSIO FILE.'
   case ("gfs_sigio")
     print*,'- INPUT DATA FROM SPECTRAL GFS SIGIO/SFCIO FILE.'
   case ("gaussian_netcdf")
     print*,'- INPUT DATA FROM FV3 GAUSSIAN NETCDF FILE.'
   case ("grib2")
     print*,'- INPUT DATA FROM A GRIB2 FILE'
   case default
     call error_handler("UNRECOGNIZED INPUT DATA TYPE.", 1)
 end select

!-------------------------------------------------------------------------
! Ensure proper file variable provided for grib2 input  
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2") then
	 if (trim(grib2_file_input_grid) == "NULL" .or. trim(grib2_file_input_grid) == "") then
		 call error_handler("FOR GRIB2 DATA, PLEASE PROVIDE GRIB2_FILE_INPUT_GRID", 1)
	 endif
 endif

 !-------------------------------------------------------------------------
! For grib2 input, warn about possibly unsupported external model types
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2") then
	 if (.not. any((/character(4)::"GFS","NAM","RAP","HRRR"/)==trim(external_model))) then
		 call error_handler( "KNOWN SUPPORTED external_model INPUTS ARE GFS, NAM, RAP, AND HRRR. " // &
		 "IF YOU WISH TO PROCESS GRIB2 DATA FROM ANOTHER MODEL, YOU MAY ATTEMPT TO DO SO AT YOUR OWN RISK. " // &
		 "ONE WAY TO DO THIS IS PROVIDE NAM FOR external_model AS IT IS A RELATIVELY STRAIGHT-" // &
		 "FORWARD REGIONAL GRIB2 FILE. YOU MAY ALSO COMMENT OUT THIS ERROR MESSAGE IN " // &
		 "program_setup.f90 LINE 389. NO GUARANTEE IS PROVIDED THAT THE CODE WILL WORK OR "// &
		 "THAT THE RESULTING DATA WILL BE CORRECT OR WORK WITH THE ATMOSPHERIC MODEL.", 1)
	 endif
 endif

!-------------------------------------------------------------------------
! For grib2 hrrr input without geogrid file input, warn that soil moisture interpolation
! will be less accurate
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2" .and. trim(external_model)=="HRRR") then
	 if (trim(geogrid_file_input_grid) == "NULL" .or. trim(grib2_file_input_grid) == "") then
		 print*, "HRRR DATA DOES NOT CONTAIN SOIL TYPE INFORMATION. WITHOUT &
			GEOGRID_FILE_INPUT_GRID SPECIFIED, SOIL MOISTURE INTERPOLATION MAY BE LESS &
			ACCURATE. "
	 endif
 endif
 
 if (trim(thomp_mp_climo_file) /= "NULL") then
   use_thomp_mp_climo=.true.
   print*,"- WILL PROCESS CLIMO THOMPSON MP TRACERS FROM FILE: ", trim(thomp_mp_climo_file)
 endif

 return

 end subroutine read_setup_namelist

!> Reads the variable mapping table, which is required for
!! initializing with GRIB2 data.
!!
!! The varmap files has entries that look like this:
!!
!! <pre>dzdt dzdt set_to_fill 0 D</pre>
!!
!! These are the chgres_var_name, field_var_name, missing_var_method,
!! missing_var_value, var_type.
!!
!! The missing_var_method is one of:
!! * set_to_fill
!! * skip
!! * stop
!!
!! The var_type is one of:
!! * T - tracer.
!! * D - variables processed by atmosphere subroutine that are not
!! tracers.
!! * S - variables processed by surface subroutine that are not
!! tracers.
!!
!! @author Larissa Reames
!! @author Jeff Beck
subroutine read_varmap

 implicit none

 integer                    :: istat, k, nvars
 character(len=500)         :: line
 character(len=20),allocatable  :: var_type(:)

 if (trim(input_type) == "grib2") then 

   print*,"OPEN VARIABLE MAPPING FILE: ", trim(varmap_file)
   open(14, file=trim(varmap_file), form='formatted', iostat=istat)
   if (istat /= 0) then
     call error_handler("OPENING VARIABLE MAPPING FILE", istat)
   endif

   num_tracers_input = 0
   nvars = 0

   !Loop over lines of file to count the number of variables
   do
     read(14, '(A)', iostat=istat) line !chgres_var_names_tmp(k)!, field_var_names(k) , &
                          ! missing_var_methods(k), missing_var_values(k), var_type(k)
     if (istat/=0) exit
     if ( trim(line) .eq. '' ) cycle
     nvars = nvars+1
   enddo
   if ( nvars == 0) call error_handler("VARMAP FILE IS EMPTY.", -1)

   allocate(chgres_var_names(nvars))
   allocate(field_var_names(nvars))
   allocate(missing_var_methods(nvars))
   allocate(missing_var_values(nvars))
   allocate(read_from_input(nvars))
   allocate(var_type(nvars))

   read_from_input(:) = .true.
   rewind(14)
    do k = 1,nvars
      read(14, *, iostat=istat) chgres_var_names(k), field_var_names(k) , &
                           missing_var_methods(k), missing_var_values(k), var_type(k)
     if (istat /= 0) call error_handler("READING VARIABLE MAPPING FILE", istat)
     if(trim(var_type(k))=='T') then
       num_tracers_input = num_tracers_input + 1
       tracers_input(num_tracers_input)=chgres_var_names(k)
       if ((trim(chgres_var_names(k)) == "ice_aero" .or. trim(chgres_var_names(k)) == "liq_aero") .and. &
           trim(thomp_mp_climo_file) .ne. "NULL" .and. trim(input_type) == "grib2") then
           call error_handler("VARMAP TABLE CONTAINS TRACER ENTRIES FOR THOMPSON AEROSOLS liq_aero or "// &
           "ice_aero. REMOVE THESE ENTRIES OR REMOVE THE NAMELIST ENTRY FOR "// &
           "thomp_mp_climo_file AND TRY AGAIN.",1)
       endif
     endif
    enddo
   close(14)
   num_tracers = num_tracers_input
 endif
end subroutine read_varmap

!> Search the variable mapping table to find conditions for handling 
!! missing variables.  Only applicable when using GRIB2 data as
!! input.
!!
!! @param [in] var_name  table variable name to search for
!! @param [out] this_miss_var_method  the method used to replace missing data
!! @param [out] this_miss_var_value  the value used to replace missing data
!! @param [out] this_field_var_name  name of variable in output file. not
!!                                   currently implemented.
!! @param [out] loc  variable table location index
!! @author Larissa Reames
!! @author Jeff Beck
subroutine get_var_cond(var_name,this_miss_var_method,this_miss_var_value, &
                            this_field_var_name, loc)
  use esmf
  
  implicit none
  character(len=20), intent(in) :: var_name
  
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

!> Driver routine to compute soil parameters for each
!! soil type. Works for Zobler and STATSGO soil categories.
!!
!! The calculations are those used in the Noah Land Surface Model. For
!! more information see <a
!! href="https://doi.org/10.1029/2002JD003296">Implementation of Noah
!! land surface model advances in the National Centers for
!! Environmental Prediction operational mesoscale Eta model</a>.
!!
!! For more details about the soil parameters/properties see <a
!! href="https://journals.ametsoc.org/view/journals/mwre/129/4/1520-0493_2001_129_0569_caalsh_2.0.co_2.xml">Coupling
!! an Advanced Land Surface–Hydrology Model with the Penn State–NCAR
!! MM5 Modeling System. Part I: Model Implementation and
!! Sensitivity</a>.
!!
!! The original source for soil properties is here:
!!
!! Cosby, B. J., G. M. Hornberger, R. B. Clapp, and T. R. Ginn, 1984:
!! <a
!! href="https://agupubs.onlinelibrary.wiley.com/doi/10.1029/WR020i006p00682">A
!! statistical exploration of the relationships of soil moisture
!! characteristics to the physical properties of soils</a>. Water
!! Resour. Res.,20, 682–690.
!!
!! The parameters in this subroutine were copied from
!! https://github.com/HelinWei-NOAA/ccpp-physics/blob/master/physics/set_soilveg.f
!! values need to be kept in sync with set_soilveg.f.
!!
!! For more information about these parameters see
!! https://github.com/HelinWei-NOAA/ccpp-physics/blob/master/physics/sflx.f.
!!
!! @param [in] localpet  ESMF local persistent execution thread
!! @author George Gayno NCEP/EMC
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
 real, allocatable         :: satdk(:)
 real, allocatable         :: satpsi(:)
 real, allocatable         :: satdw(:)

! using stasgo table
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
   case ("gfs_sigio")
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

 select case (trim(input_type))
   case ("gfs_sigio")
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
                       bb, satpsi, satdw, refsmc_input, drysmc_input, wltsmc_input)

 deallocate(bb, satdk, satpsi, satdw)

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

 smlow  = smlow_statsgo
 smhigh = smhigh_statsgo
 maxsmc_target = maxsmc_statsgo
 bb_target     = bb_statsgo
 satdk  = satdk_statsgo
 satpsi_target = satpsi_statsgo

 call calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc_target, &
                       bb_target, satpsi_target, satdw, refsmc_target, drysmc_target, wltsmc_target)

 deallocate(satdk, satdw)

 if (localpet == 0) print*,'maxsmc target grid ',maxsmc_target
 if (localpet == 0) print*,'wltsmc input grid ',wltsmc_target

 end subroutine calc_soil_params_driver

!> Compute soil parameters. These will be used to rescale soil
!! moisture differences in soil type between the input grid and target
!! model grid.
!!
!! @param [in] num_soil_cats  number of soil type categories
!! @param [in] smlow  reference parameter for wltsmc
!! @param [in] smhigh reference parameter for refsmc
!! @param [in] satdk  saturated soil moisture hydraulic conductivity
!! @param [in] maxsmc maximum soil moisture (porosity)
!! @param [in] bb  soil 'b' parameter
!! @param [in] satpsi saturated soil potential
!! @param [out] satdw  saturated soil diffusivity/conductivity coefficient
!! @param [out] refsmc  onset of soil moisture stress (field capacity)
!! @param [out] drysmc  air dry soil moisture limit
!! @param [out] wltsmc  plant soil moisture wilting point
!! @author George Gayno NCEP/EMC
 subroutine calc_soil_params(num_soil_cats, smlow, smhigh, satdk,  &
            maxsmc, bb, satpsi, satdw, refsmc, drysmc, wltsmc)

 implicit none

 integer, intent(in)            :: num_soil_cats

 real, intent(in)               :: smlow, smhigh
 real, intent(in)               :: bb(num_soil_cats)
 real, intent(in)               :: maxsmc(num_soil_cats)
 real, intent(in)               :: satdk(num_soil_cats)
 real, intent(in)               :: satpsi(num_soil_cats)

 real, intent(out)              :: satdw(num_soil_cats)
 real, intent(out)              :: refsmc(num_soil_cats)
 real, intent(out)              :: drysmc(num_soil_cats)
 real, intent(out)              :: wltsmc(num_soil_cats)

 integer                        :: i

 real                      :: refsmc1
 real                      :: wltsmc1

 satdw = 0.0
 refsmc = 0.0
 wltsmc = 0.0
 drysmc = 0.0

 do i = 1, num_soil_cats

   if (maxsmc(i) > 0.0) then

   SATDW(I)  = BB(I)*SATDK(I)*(SATPSI(I)/MAXSMC(I))
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
