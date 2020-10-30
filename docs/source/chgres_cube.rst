.. _chgres_cube:

************************
Introduction
************************

The chgres_cube program creates initial condition files to “coldstart” the forecast model.  The initial conditions are created from either Global Forecast System (GFS) gridded binary version 2 (GRIB2), NOAA Environmental Modeling System Input/Output (NEMSIO) data, or Network Common Data Form (NetCDF) data.

************************************************
Where to find GFS GRIB2, NEMSIO and NetCDF data
************************************************

**GRIB2:**

      * 0.25-degree data (last 10 days only) - Use the **gfs.tHHz.pgrb2.0p25.f000** files in subdirectory gfs.YYYYMMDD/HH `here <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>`_.`

      * 0.5-degree data - Use the **gfs_4_YYYYMMDD_00HH_000.grb2** file, under **GFS Forecasts 004 (0.5-deg)** here: `NCDC - Global Forecast System <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs>`__.  Note: *Tests were not done with the AVN, MRF or analysis data*.

      * 1.0-degree data - Use the **gfs_3_YYYYMMDD_00HH_000.grb2 file**, under **GFS Forecasts 003 (1-deg)** here: `NCDC - Global Forecast System <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs>`__.  Note: *Tests were not done with the AVN, MRF or analysis data*.

**NEMSIO:**

      * T1534 gaussian (last 10 days only) - Use the **gfs.tHHz.atmanl.nemsio** (atmospheric fields) and **gfs.tHHz.sfcanl.nemsio** (surface fields) files in subdirectory gfs.YYYYMMDD/HH `here <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>`_.

**NetCDF:**

      * T1534 gaussian (don't have any more details at this time).

************************************************
Initializing with GRIB2 data - some caveats
************************************************


Keep this in mind when using GFS GRIB2 data for model initialization:

      * GRIB2 data does not contain the fields needed for the Near Sea Surface Temperature (NSST) scheme.  See the next section for options on running the forecast model in this situation.
      * Data is coarse (in vertical and horizontal) compared to the NCEP operational GFS .  May not provide a good initialization (especially for the surface).  Recommendations:

	      * C96 - use 0.25, 0.5 or 1.0-degree GRIB2 data
	      * C192 - use 0.25 or 0.5-degree GRIB2 data
	      * C384 - use 0.25-degree GRIB2 data
	      * C768 - try the 0.25-degree GRIB2 data.  But it may not work well.
      * Sea/lake ice thickness and column temperatures are not available.  So, nominal values of 1.5 m and 265 K are used.
      * Soil moisture in the GRIB2 files is created using bilinear interpolation and, therefore, may be a mixture of values from different soil types.  Could result in poor latent/sensible heat fluxes.
      * Ozone is not available at all isobaric levels.  Missing levels are set to a nominal value defined in the variable mapping (VARMAP) file (1E-07).
      * Only tested with GRIB2 data from GFS v14 and v15 (from 12z July 19, 2017 to current).  May not work with older GFS data.  Will not work with GRIB2 data from other models.

**Near Sea Surface Temperature (NSST) data and GRIB2 initialization**

The issue with not having NSST data is important.  In GFS we use the foundation temperature (Tref) and add a diurnal warming/cooling layer using NSST. This is the surface temperature that is passed to the atmospheric boundary layer. This is a critical feature, especially when we are doing Data Assimilation.

When using NEMSIO or NetCDF data to initialize the model, both the foundation and surface temperature are available and the atmospheric model should be run using the NSST option as this will properly account for in the forward run of the model.

In GRIB2 files only the Tsfc is stored and that is set as foundation temperature as well. So the diurnal heating / cooling is baked into the initial condition for the extent of the run. This can be critical if the model is being initialized when the ocean is warm and initialization is occuring at the peak of the diurnal warming. That warm ocean will be baked in for the extent of the run and may spawn off a number of fake hurricanes. The user has two options -- either to use a spin up cycle to spin up NSST (set **nstf_name** = [2,1,0,0,0] in **input.nml** of the model namelist file. This will create a diurnal cycle after 24 hours of spin up), or to run the model without any NSST option ( set **nstf_name** = [0,0,0,0,0] in **input.nml** of the model namelist file. The user will also have to choose one of the no NSST physics suite options in **input.nml**).

Note, that neither of these two options will get rid of the underlying baked in heating/cooling in the surface temperature fields. For most cases this may not be an issue, but where it is then the user will either have to initialize the model with NEMSIO or NetCDF data, or replace the surface temperature in the GRIB2 fields with independently obtained foundation temperature.

************************************************
chgres_cube namelist options
************************************************

Namelist variables with “input” in their name refer to data input to chgres_cube.  Namelist variables with “target” in their name refer to the FV3 horizontal and vertical grid (i.e., the target grid chgres_cube is mapping to).

When using GRIB2 data as input to chgres_cube, set namelist as follows:

      * fix_dir_target_grid - Path to the tiled FV3 surface climatological files (such as albedo).
      * mosaic_file_target_grid - Path and name of the FV3 mosaic file.
      * orog_dir_target_grid - directory containing the tiled FV3 orography and grid files (NetCDF).
      * orog_files_target_grid - names of the six tiled FV3 orography files.
      * vcoord_file_target_grid - path and name of the model vertical coordinate definition file (“global_hyblev.l$LEVS.txt).
      * data_dir_input_grid - directory containing the GRIB2 initial conditions data
      * grib2_file_input_grid - name of the GRIB2 input data file
      * varmap_file - path and name of the variable mapping (VARMAP) table.  See below for details on this table.
      * input_type - input data type. Set to ‘grib2’
      * cycle_mon/day/hour - month/day/hour of your model initialization
      * convert_atm - set to ‘true’ to process the atmospheric fields
      * convert_sfc - set to ‘true’ to process the surface fields

When using NEMSIO data as input to chgres_cube, set namelist as follows:

      * fix_dir_target_grid - Path to the tiled FV3 surface climatological files (such as albedo).
      * mosaic_file_target_grid - Path and name of the FV3 mosaic file.
      * orog_dir_target_grid - directory containing the tiled FV3 orography and grid files (NetCDF).
      * orog_files_target_grid - names of the six tiled FV3 orography files.
      * vcoord_file_target_grid - path and name of the model vertical coordinate definition file (“global_hyblev.l$LEVS.txt).
      * data_dir_input_grid - directory containing the NEMSIO input data
      * atm_files_input_grid - name of the NEMSIO input atmospheric data file
      * sfc_files_input_grid - name of the NEMSIO input surface/Near Sea Surface Temperature (NSST) data file
      * input_type - input data type. Set to ‘gaussian_nemsio’.
      * cycle_mon/day/hour - month/day/hour of your model run
      * convert_atm - set to ‘true’ to process the atmospheric fields
      * convert_sfc - set to ‘true’ to process the surface fields
      * convert_nst - set to ‘true’ to process NSST fields
      * tracers_input - names of tracer records in input file.  For GFDL microphysics, set to “spfh”,”clwmr”,”o3mr”,”icmr”,”rwmr”,”snmr”,”grle”.
      * tracers - names of tracer records in output file expected by model.  For GFDL microphysics, set to “sphum”,”liq_wat”,”o3mr”,”ice_wat”,”rainwat”,”snowwat”,”graupel”.

When using NetCDF data as input to chgres_cube, set namelist as follows:

      * fix_dir_target_grid - Path to the tiled FV3 surface climatological files (such as albedo).
      * mosaic_file_target_grid - Path and name of the FV3 mosaic file.
      * orog_dir_target_grid - directory containing the tiled FV3 orography and grid files (NetCDF).
      * orog_files_target_grid - names of the six tiled FV3 orography files.
      * vcoord_file_target_grid - path and name of the model vertical coordinate definition file (“global_hyblev.l$LEVS.txt).
      * data_dir_input_grid - directory containing the NetCDF input data
      * atm_files_input_grid - name of the NetCDF input atmospheric data file
      * sfc_files_input_grid - name of the NetCDF input surface/Near Sea Surface Temperature (NSST) data file
      * input_type - input data type. Set to ‘gaussian_netcdf’.
      * cycle_mon/day/hour - month/day/hour of your model run
      * convert_atm - set to ‘true’ to process the atmospheric fields
      * convert_sfc - set to ‘true’ to process the surface fields
      * convert_nst - set to ‘true’ to process NSST fields
      * tracers_input - names of tracer records in input file.  For GFDL microphysics, set to “spfh”,”clwmr”,”o3mr”,”icmr”,”rwmr”,”snmr”,”grle”.
      * tracers - names of tracer records in output file expected by model.  For GFDL microphysics, set to “sphum”,”liq_wat”,”o3mr”,”ice_wat”,”rainwat”,”snowwat”,”graupel”.


************************
Compiling the program
************************

chgres_cube requires cmake 3.12 or higher. It can be built as part of the NCEPLIBS unified build system that includes two separate build systems -- one for the third party libraries that are needed and the other for the libraries and utilities themselves. See https://github.com/NOAA-EMC/NCEPLIBS-external/wiki for more detailed information.

If the NCEPLIBS have been installed and the user wants to compile chgres_cube again

      * make sure paths are set to hdf5, compiler, mpi and cmake
      * In a bash environment run
	      * cd /path/to/nceplibs/installed
	      * source bin/setenv_nceplibs.sh (this will set all necessary environments)
      * set cmake compiler - export FC=ifort (if ifort is the compiler chosen)
      * cd to where you checked out the UFS_Utils
      * mkdir build and cd build
      * cmake .. -DCMAKE_INSTALL_PREFIX=/path/where/you/want/the/code/installed -DCMAKE_PREFIX_PATH=/path/to/nceplibs/installed
      * make -j x (where x is a number that can be chosen to speed up the make, usually 8)
      * make install
      * if you do get errors that cmake cannot find "FindNETCDF" or "FindESMF", run: git submodule update --init --recursive

************************************************
Program inputs and outputs
************************************************

**Inputs**

The following four sets of files are located here: https://ftp.emc.ncep.noaa.gov/EIB/UFS/global/fix/fix_fv3_gmted2010.v20191213/

      * FV3 mosaic file - (NetCDF format)
	      * CRES_mosaic.nc

      * FV3 grid files - (NetCDF format)
	      * CRES_grid.tile1.nc
	      * CRES_grid.tile2.nc
	      * CRES_grid.tile3.nc
	      * CRES_grid.tile4.nc
	      * CRES_grid.tile5.nc
	      * CRES_grid.tile6.nc

      * FV3 orography files - (NetCDF format)
	      * CRES_oro_data.tile1.nc
	      * CRES_oro_data.tile2.nc
	      * CRES_oro_data.tile3.nc
	      * CRES_oro_data.tile4.nc
	      * CRES_oro_data.tile5.nc
	      * CRES_oro_data.tile6.nc

      * FV3 surface climatological files - Located under the ./fix_sfc sub-directory.  One file for each tile.  NetCDF format.
	      * CRES.facsf.tileX.nc (fractional coverage for strong/weak zenith angle dependent albedo)
	      * CRES.maximum_snow_albedo.tileX.nc (maximum snow albedo)
	      * CRES.slope_type.tileX.nc (slope type)
	      * CRES.snowfree_albedo.tileX.nc (snow-free albedo)
	      * CRES.soil_type.tileX.nc (soil type)
	      * CRES.subtrate_temperature.tileX.nc (soil substrate temperature)
	      * CRES.vegetation_greenness.tileX.nc (vegetation greenness)
	      * CRES.vegetation_type.tileX.nc (vegetation type)

      * FV3 vertical coordinate file.  Text file.  Located here: https://ftp.emc.ncep.noaa.gov/EIB/UFS/global/fix/fix_am.v20191213/
	      * global_hyblev.l$LEVS.txt

      * Input data files.  GRIB2, NEMSIO or NetCDF.  See above section for how to find this data.

**Outputs**

      * Atmospheric “coldstart” files.  NetCDF.
	      * out.atm.tile1.nc
	      * out.atm.tile2.nc
	      * out.atm.tile3.nc
	      * out.atm.tile4.nc
	      * out.atm.tile5.nc
	      * out.atm.tile6.nc

      * Surface/Near Sea Surface Temperature (NSST) “coldstart” files.  NetCDF
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc

************************************************
Running the program stand alone
************************************************

      * Locate your input files.  See above for a list.
      * Set the namelist for your experiment.  See above for an explanation of the namelist entries.
      * Link the namelist to Fortran unit number 41, i.e.”
	      * ln -fs your-namelist-file  ./fort.41
      * Load any required runtime libraries.  For example, you may need to load libraries for NetCDF and/or your Fortran compiler.
      * Run the program with an MPI task count that is a multiple of six.  This is an ESMF library requirement when processing a six-tiled global grid.

************************
Code structure
************************

Note on variable names: “input” refers to the data input to the program (i.e., GRIB2, NEMSIO, NetCDF).  “Target” refers to the target or FV3 model grid.  See routine doc blocks for more details.

      * chgres.F90 - This is the main driver routine.
      * program_setup.F90 - Sets up the program execution.
	      * Reads program namelist
	      * Computes required soil parameters
	      * Reads the variable mapping (VARMAP) table.
      * model_grid.F90 - Sets up the ESMF grid objects for the input data grid and target FV3 grid.
      * static_data.F90 - Reads static surface climatological data for the target FV3 grid (such as soil type and vegetation type).  Time interpolates time-varying fields, such as monthly plant greenness, to the model run time.  Data for each target FV3 resolution resides in the ‘fixed’ directory.  Set path via the fix_dir_target_grid namelist variable.
      * write_data.F90 - Writes the tiled and header files expected by the forecast model.
      * input_data.F90 - Contains routines to read atmospheric and surface data from GRIB2, NEMSIO and NetCDF files.
      * utils.f90 - Contains utility routines, such as error handling.
      * grib2_util.F90 -  Routines to (1) convert from RH to specific humidity; (2) convert from omega to dzdt.  Required for GRIB2 input data.
      * atmosphere.F90 - Process atmospheric fields.  Horizontally interpolate from input to target FV3 grid using ESMF regridding.  Adjust surface pressure according to terrain differences between input and target grid.  Vertically interpolate to target FV3 grid vertical levels.  Description of main routines:
	      * read_vcoord_info - Reads model vertical coordinate definition file (as specified by namelist variable vcoord_file_target_grid).
	      * newps - computes adjusted surface pressure given a new terrain height.
	      * newpr1 - computes 3-D pressure given an adjusted surface pressure.
	      * vintg - vertically interpolate atmospheric fields to target FV3 grid.
      * surface.F90 - process land, sea/lake ice, open water/Near Sea Surface Temperature (NSST) fields.  Assumes the input land data are Noah LSM-based, and the fv3 run will use the Noah LSM.   NSST fields are not available when using GRIB2 input data.  Description of main routines:
	      * interp - horizontally interpolate fields from input to target FV3 grid.
	      * calc_liq_soil_moisture - compute liquid portion of total soil moisture.
	      * adjust_soilt_for_terrain - adjust soil temperature for large differences between input and target FV3 grids.
	      * rescale_soil_moisture - adjust total soil moisture for differences between soil type on input and target FV3 grids.  Required to preserve latent/sensible heat fluxes.  Assumes Noah LSM.
	      * roughness - set roughness length at land and sea/lake ice.  At land, a vegetation type-based lookup table is used.
	      * qc_check - some consistency checks.
      * search_util.f90 - searches for the nearest valid land/non-land data where the input and target fv3 land-mask differ.  Example: when the target FV3 grid depicts an island that is not resolved by the input data.  If nearby valid data is not found, a default value is used.

************************************************
Making changes to the chgres_cube program
************************************************

chgres_cube is part of the UFS_UTILS repository (https://github.com/NOAA-EMC/UFS_UTILS). When wanting to contribute to this repository developers shall follow the Gitflow software development process

      * Developers shall create their own fork of the UFS_UTILS repository
      * Developers shall create a ‘feature’ branch off ‘develop’ in their fork for all changes.
      * Developers shall open an issue and reference it in all commits.

For more details, see the UFS_UTILS wiki page: https://github.com/NOAA-EMC/UFS_UTILS/wiki

Changes that support current or future NCEP operations will be given priority for inclusion into the authoritative repository.

************************************************
Variable Mapping (VARMAP) table
************************************************

The VARMAP table, set in the chgres_cube namelist (variable varmap_file), controls how chgres_cube handles variables that might be missing from the GRIB2 files. Since there are so many different versions of GRIB2 files, it's often uncertain what fields are available even if you know what source model the data is coming from.  Each file contains the following:  (Note, only the GFS physics suite is currently supported.)

Column 1: Name the code searches for in the table. Do not change.  Some definitions:

      * dzdt - vertical velocity
      * sphum - specific humidity
      * liq_wat - liquid water mixing ratio
      * o3mr - ozone mixing ratio
      * ice_wat - ice water mixing ratio
      * rainwat - rain water mixing ratio
      * snowwat - snow water mixing ratio
      * graupel - graupel mixing ratio
      * vtype - vegetation type
      * sotype - soil type
      * vfrac - plant greenness fraction
      * fricv - friction velocity
      * sfcr - roughness length
      * tprcp - precipitation rate
      * ffmm - surface exchange coefficient for momentum
      * ffhh - surface exchange coefficient for heat
      * f10m - log((sfcr+10)/sfcr)
      * soilw - total volumetric soil moisture
      * soill - liquid volumetric soil moisture
      * soilt - soil column temperature
      * cnwat - plant canopy water content
      * hice - sea/lake ice thickness
      * weasd - snow liquid equivalent
      * snod - physical snow depth

Column 2: Name of the variable in the output “coldstart” files. Unimplemented.

Column 3: Behavior when the code can't find the variable in the input file. Options are:

      * "skip": Don't write to the output file.
      * "set_to_fill": Set to user-specified field value (see column 4).
      * "stop": Force an exception and stop code execution. Use this if you absolutely require a field to be present.

Column 4: If column 3 = "set_to_fill", then this value is used to fill in all points in the input field. These values may be overwritten by the code before output depending on the variable (especially for surface variables).

Column 5: Variable type descriptor. Variable names designated as tracers are used to populate the list of tracers to read from the GRIB2 file and write to output, so make sure all tracers you wish to read have an entry. Note that if you wish to add a tracer name that is not already included in the appropriate VARMAP file, this will require modification of the chgres_cube code. Valid choices are:

      * “T”: 3-dimensional tracer array
      * “D”: 3-dimensional non-tracer array
      * “S”: 2-dimensional surface array
