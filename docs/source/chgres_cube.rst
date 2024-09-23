.. _chgres_cube:

Introduction
------------

The chgres_cube program creates initial condition files to coldstart the forecast model.  The initial conditions are created from either Finite-Volume Sphere (FV3) Global Forecast System (GFS), North American Mesoscale Forecast System (NAM), Rapid Refresh (RAP), or High Resolution Rapid Refresh (HRRR) gridded binary version 2 (GRIB2) data.

Code structure
--------------

Note on variable names: “input” refers to the data input to the program (i.e., GRIB2, NEMSIO, NetCDF).  “Target” refers to the target or FV3 model grid.  See routine doc blocks for more details.

The program assumes Noah/Noah-MP LSM coefficients for certain soil thresholds. In the future, an option will be added to use RUC LSM thresholds.

      * chgres.F90 - This is the main driver routine.
      * program_setup.F90 - Sets up the program execution.

            * Reads program namelist
            * Computes required soil parameters
            * Reads the variable mapping (VARMAP) table.
      * model_grid.F90 - Sets up the ESMF grid objects for the input data grid and target FV3 grid.
      * static_data.F90 - Reads static surface climatological data for the target FV3 grid (such as soil type and vegetation type).  Time interpolates time-varying fields, such as monthly plant greenness, to the model run time. Set path to these files via the fix_dir_target_grid namelist variable.
      * write_data.F90 - Writes the tiled and header files expected by the forecast model.
      * atm_input_data.F90 - Contains routines to read input atmospheric data from GRIB2, NEMSIO and NetCDF files.
      * nst_input_data.F90 - Contains routines to read input NSST data from NEMSIO and NetCDF files.
      * sfc_input_data.F90 - Contains routines to read input surface data from GRIB2, NEMSIO and NetCDF files.
      * utils.F90 - Contains utility routines, such as error handling.
      * grib2_util.F90 -  Routines to (1) convert from RH to specific humidity; (2) convert from omega to dzdt.  Required for GRIB2 input data.
      * atmosphere.F90 - Process atmospheric fields.  Horizontally interpolate from input to target FV3 grid using ESMF regridding.  Adjust surface pressure according to terrain differences between input and target grid.  Vertically interpolate to target FV3 grid vertical levels.  Description of main routines:

            * read_vcoord_info - Reads model vertical coordinate definition file (as specified by namelist variable vcoord_file_target_grid).
            * newps - Computes adjusted surface pressure given a new terrain height.
            * newpr1 - Computes 3-D pressure given an adjusted surface pressure.
            * vintg - vertically interpolate atmospheric fields to target FV3 grid.
            * vintg_wam - vertically interpolate atmospheric fields to the thermosphere. Supports the Whole Atmosphere Model.
      * atmosphere_target_data.F90 - Holds the target grid atmospheric ESMF fields.
      * surface.F90 - process land, sea/lake ice, open water/Near Sea Surface Temperature (NSST) fields.  NSST fields are not available when using GRIB2 input data.  Description of main routines:

            * interp - horizontally interpolate fields from input to target FV3 grid.
            * calc_liq_soil_moisture - compute liquid portion of total soil moisture.
            * adjust_soilt_for_terrain - adjust soil temperature for large differences between input and target FV3 grids.
            * rescale_soil_moisture - adjust total soil moisture for differences between soil type on input and target FV3 grids.  Required to preserve latent/sensible heat fluxes.
            * roughness - set roughness length at land and sea/lake ice.  At land, a vegetation type-based lookup table is used.
            * qc_check - some consistency checks.
      * surface_target_data.F90 - Holds the target grid surface ESMF fields.
      * search_util.F90 - searches for the nearest valid land/non-land data where the input and target fv3 land-mask differ.  Example: when the target FV3 grid depicts an island that is not resolved by the input data.  If nearby valid data is not found, a default value is used.
      * thompson_mp_climo_data.F90 - Processes climatological Thompson micro-physics fields.
      * wam_climo_data.f90 - Process vertical profile climatological data for the Whole Atmosphere Model.

Configuring and using chgres_cube for global applications
---------------------------------------------------------

Program inputs and outputs for global applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Inputs**

Users may create their own global grids, or use the pre-defined files located in the `./CRES directories <https://noaa-nws-global-pds.s3.amazonaws.com/index.html#fix/orog/20231027/>`_. (where CRES is the atmospheric resolution and mxRES is the ocean resolution).

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
	      * CRES.mxRES_oro_data.tile1.nc
	      * CRES.mxRES_oro_data.tile2.nc
	      * CRES.mxRES_oro_data.tile3.nc
	      * CRES.mxRES_oro_data.tile4.nc
	      * CRES.mxRES_oro_data.tile5.nc
	      * CRES.mxRES_oro_data.tile6.nc

      * FV3 surface climatological files - Located under the `./CRES/sfc <https://noaa-nws-global-pds.s3.amazonaws.com/index.html#fix/orog/20231027/>`_ subdirectories. One file for each tile.  NetCDF format.
	      * CRES.mxRES.facsf.tileX.nc (fractional coverage for strong/weak zenith angle dependent albedo)
	      * CRES.mxRES.maximum_snow_albedo.tileX.nc (maximum snow albedo)
	      * CRES.mxRES.slope_type.tileX.nc (slope type)
	      * CRES.mxRES.snowfree_albedo.tileX.nc (snow-free albedo)
	      * CRES.mxRES.soil_type.tileX.nc (soil type)
	      * CRES.mxRES.subtrate_temperature.tileX.nc (soil substrate temperature)
	      * CRES.mxRES.vegetation_greenness.tileX.nc (vegetation greenness)
	      * CRES.mxRES.vegetation_type.tileX.nc (vegetation type)

      * FV3 vertical coordinate file.  Text file.  `Located here <https://noaa-nws-global-pds.s3.amazonaws.com/index.html#fix/am/20220805/>`_.
	      * global_hyblev.l$LEVS.txt

      * Input data files.  GRIB2, NEMSIO or NetCDF.  See the next section for how to find this data.

**Outputs**

      * Atmospheric coldstart files.  NetCDF.
	      * out.atm.tile1.nc
	      * out.atm.tile2.nc
	      * out.atm.tile3.nc
	      * out.atm.tile4.nc
	      * out.atm.tile5.nc
	      * out.atm.tile6.nc

      * Surface/Near Sea Surface Temperature (NSST) coldstart files.  NetCDF
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc
	      * out.sfc.tile1.nc


Where to find GFS GRIB2 and NetCDF data for global applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**GRIB2**

      * 0.25-degree data (last 10 days only) - Use the **gfs.tHHz.pgrb2.0p25.f000** files in subdirectory ./gfs.YYYYMMDD/HH/atmos on `NOMADS <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>`_.

      * 0.5-degree data - Use the **gfs_4_YYYYMMDD_HHHH_000.grb2** file, under **GFS Forecasts 004 (0.5-deg)** here: `NCEI - Global Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast>`_.  Note: *Tests were not done with the AVN, MRF or analysis data*.

      * 1.0-degree data - Use the **gfs_3_YYYYMMDD_HHHH_000.grb2 file**, under **GFS Forecasts 003 (1.0-deg)** here: `NCEI - Global Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast>`_.  Note: *Tests were not done with the AVN, MRF or analysis data*.

**NetCDF**

      * T1534 gaussian (last 10 days only) - Use the **gfs.tHHz.atmanl.nc** (atmospheric fields) and **gfs.tHHz.sfcanl.nc** (surface fields) files in subdirectory ./gfs.YYYYMMDD/HH/atmos on `NOMADS <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>`_.

Initializing global domains with GRIB2 data - some caveats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Keep these things in mind when using GFS GRIB2 data for model initialization.**

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
      * Note that when concatenating grib2 files for use in initialization of global simulations, it is possible to inadvertently introduce duplicate variables and levels into the subsequent grib2 files.  Chgres_cube will automatically fail with a warning message indicating that the grib2 file used contains these duplicate entries.  Prior to continuing it will be necessary to strip out duplicate entries.  Users can remove these entries through use of wgrib2, such as in the following command:
              * ``wgrib2 IN.grb -submsg 1 | unique.pl | wgrib2 -i IN.grb -GRIB OUT.grb``, where IN.grb is the original concatenated grib2 file, and OUT.grb is the resulting grib2 file, with duplicates removed.  The "unique.pl" Perl script is as follows, taken from the `Tricks for wgrib2 <https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/tricks.wgrib2>`_ website:

                      .. code-block:: console

                         ----------------------- unique.pl ------------------------
                         #!/usr/bin/perl -w
                         # print only lines where fields 3..N are different
                         # 
                         while (<STDIN>) {
                            chomp;
                            $line = $_;
                            $_ =~ s/^[0-9.]*:[0-9]*://;
                            if (! defined $inv{$_}) { 
                              $inv{$_} = 1;
                              print "$line\n";
                            }
                         }
                         --------------------- end unique.pl ----------------------

Near Sea Surface Temperature (NSST) data and GRIB2 initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The issue with not having NSST data is important.  In GFS we use the foundation temperature (Tref) and add a diurnal warming/cooling layer using NSST. This is the surface temperature that is passed to the atmospheric boundary layer. This is a critical feature, especially when we are doing Data Assimilation.

When using NEMSIO or NetCDF data to initialize the model, both the foundation and surface temperature are available and the atmospheric model should be run using the NSST option as this will properly account for in the forward run of the model.

In GRIB2 files only the Tsfc is stored and that is set as foundation temperature as well. So the diurnal heating / cooling is baked into the initial condition for the extent of the run. This can be critical if the model is being initialized when the ocean is warm and initialization is occuring at the peak of the diurnal warming. That warm ocean will be baked in for the extent of the run and may spawn off a number of fake hurricanes. The user has two options -- either to use a spin up cycle to spin up NSST (set **nstf_name** = [2,1,0,0,0] in **input.nml** of the model namelist file. This will create a diurnal cycle after 24 hours of spin up), or to run the model without any NSST option ( set **nstf_name** = [0,0,0,0,0] in **input.nml** of the model namelist file. The user will also have to choose one of the no NSST physics suite options in **input.nml**).

Note, that neither of these two options will get rid of the underlying baked in heating/cooling in the surface temperature fields. For most cases this may not be an issue, but where it is then the user will either have to initialize the model with NEMSIO or NetCDF data, or replace the surface temperature in the GRIB2 fields with independently obtained foundation temperature.

Global chgres_cube namelist options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Namelist variables with “input” in their name refer to data input to chgres_cube.  Namelist variables with “target” in their name refer to the FV3 horizontal and vertical grid (i.e., the target grid chgres_cube is mapping to).

Namelist settings for using **GRIB2** data as input in global chgres_cube applications 

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

Namelist settings for using **NEMSIO** data as input in global chgres_cube applications

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

Namelist settings for using **NetCDF** data as input in global chgres_cube applications 

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

Configuring and using chgres_cube for regional applications
----------------------------------------------------------------

Regional program inputs and outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Inputs**

The following four sets of files/directories should all be located in the same directory (orog_dir_target_grid in the namelist):

      * FV3 mosaic file - (NetCDF format)
	      * CRES_mosaic.halo4.nc

      * FV3 grid files - (NetCDF format)
	      * CRES_grid.tile7.halo4.nc 

      * FV3 orography files - (NetCDF format)
	      * CRES_oro_data.tile7.halo4.nc

      * FV3 surface climatological files - NetCDF format.  Linked without the “halo4” (e.g., CRES.facsf.tile7.halo4.nc linked as CRES.facsf.tile7.nc)
	      * CRES.facsf.tile7.halo4.nc (fractional coverage for strong/weak zenith angle dependent albedo)
	      * CRES.maximum_snow_albedo.tile7.halo4.nc (maximum snow albedo)
	      * CRES.slope_type.tile7.halo4.nc (slope type)
	      * CRES.snowfree_albedo.tile7.halo4.nc (snow-free albedo)
	      * CRES.soil_type.tile7.halo4.nc (soil type)
	      * CRES.subtrate_temperature.tile7.halo4.nc (soil substrate temperature)
	      * CRES.vegetation_greenness.tile7.halo4.nc (vegetation greenness)
	      * CRES.vegetation_type.tile7.halo4.nc (vegetation type)

      * FV3 vertical coordinate file.  Text file.  `Located here <https://noaa-nws-global-pds.s3.amazonaws.com/index.html#fix/am/20220805/>`_.
	      * global_hyblev.l$LEVS.txt

      * Input data files. GRIB2 only.  See the next section for how to find this data.

**Outputs**

      * Atmospheric coldstart file.  NetCDF.
        * out.atm.tile7.nc

      * Surface coldstart file.  NetCDF.
        * out.sfc.tile7.nc

Where to find FV3GFS, NAM, HRRR, and RAP GRIB2 data for regional applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**FV3GFS**

      * 0.25-degree data (last 10 days only) - Use the **gfs.tHHz.pgrb2.0p25.f000** files in subdirectory ./gfs.YYYYMMDD/HH/atmos on `NOMADS <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>`_.

      * 0.5-degree data - Use the **gfs_4_YYYYMMDD_HHHH_000.grb2** file, under **GFS Forecasts 004 (0.5-deg)** here: `NCEI - Global Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast>`_.  Note: *Tests were not done with the AVN, MRF or analysis data*.

      * 1.0-degree data - Use the **gfs_3_YYYYMMDD_HHHH_000.grb2 file**, under **GFS Forecasts 003 (1.0-deg)** here: `NCEI - Global Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast>`_.  Note: *Tests were not done with the AVN, MRF or analysis data*.

**NAM**

     * 12-km data from last few days (NOMADS) - Use the **nam.tHHz.conusnest.hiresfHH.tmHH.grib2** files in subdirectory nam.YYYYMMDD on `NOMADS <https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/>`_.

     * 12-km data starting 2020 - Use the **nam_218_YYYYMMDD_HHHH_000.grb2 file**, under **NAM Forecasts NAM-NMM 218 (12km) Domain** here: `NCEI - North American Mesoscale Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/north-american-mesoscale>`_.

     * 12-km archived data prior to 2020 can be requested through the Archive Information Request System `here <https://www.ncei.noaa.gov/has/HAS.FileAppRouter?datasetname=NAM218&subqueryby=STATION&applname=&outdest=FILE>`_.

**HRRR**
 
      * 3-km operational data from previous few days (NOMADS) - Use the **hrrr.tHHz.wrfnatfHH.grib2** files in the subdirectory ./hrrr.YYYYMMDD/conus `here <https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/>`_.

      * 3-km operational data from 2015 to present (AWS S3): Go `here <https://registry.opendata.aws/noaa-hrrr-pds/>`__ and click “Browse Bucket.” Type "YYYYMMDD" in to the Search bar. Use the **hrrr.t00z.wrfnatf00.grib2** files in the directory hrrr.YYYYMMDD/conus/.

      * 3-km operational data from 2015 to present (Google Cloud): Go `here <https://console.cloud.google.com/marketplace/product/noaa-public/hrrr>`__ and click “View Dataset.” Type “hrrr.YYYYMMDD” into the “Filter” box. Use the **hrrr.tHHz.wrfnatfFF.grib2** files in the hrrr.YYYYMMDD/conus directory.

      * 3-km operational data from 2016 to present (University of Utah): `Click here <http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/hrrr_download.cgi>`__.

**RAP**

      * 13-km operational data for the previous few days (NOMADS): Use the **rap.tHHz.wrfnatfHH.grib2** files in the subdirectory ./rap.YYYYMMDD `here <https://nomads.ncep.noaa.gov/pub/data/nccf/com/rap/prod/>`_.

      * 13-km isobaric level data from previous 6 months : Use the **rap_130_YYYYMMDD_HHHH_0FF.grb2** files from **RAP Forecasts - RAP 130 (13km) - Domain** at NCEI `here <https://www.ncei.noaa.gov/products/weather-climate-models/rapid-refresh-update>`_.

      * 13-km archived isobaric data older than 6 months can be requested through the Archive Information Request System `here <https://www.ncei.noaa.gov/has/HAS.FileAppRouter?datasetname=RAP130&subqueryby=STATION&applname=&outdest=FILE>`_.


Initializing regional domains with GRIB2 data - some caveats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Keep these things in mind when using FV3GFS GRIB2 data for model initialization:

      * GRIB2 data does not contain the fields needed for the Near Sea Surface Temperature (NSST) scheme.  
      * External model recommendations for pre-defined CONUS grids:

              * 3-km domain, HRRR or RAP data is recommended
              * 13-km domain: RAP or GFS data is recommended
              * 25-km domain: GFS data is recommended
      * Sea/lake ice thickness and column temperatures are not available.  So, nominal values of 1.5 m and 265 K are used.
      * For FV3GFS GRIB2 data, soil moisture is created using bilinear interpolation and, therefore, may be a mixture of values from different soil types. Could result in poor latent/sensible heat fluxes.
      * Ozone is not available at all isobaric levels. Missing levels are set to a nominal value defined in the variable mapping (VARMAP) file (1E-07).
      * Only tested with GRIB2 data from FV3GFS, RAP, NAM, and HRRR data. May not work with GRIB2 data from other models. Use these at your own risk.
      * Note that when concatenating grib2 files for use in initialization of regional simulations, it is possible to inadvertently introduce duplicate variables and levels into the subsequent grib2 files.  Chgres_cube will automatically fail with a warning message indicating that the grib2 file used contains these duplicate entries.  Prior to continuing it will be necessary to strip out duplicate entries.  Users can remove these entries through use of wgrib2, such as in the following command:
              * ``wgrib2 IN.grb -submsg 1 | unique.pl | wgrib2 -i IN.grb -GRIB OUT.grb``, where IN.grb is the original concatenated grib2 file, and OUT.grb is the resulting grib2 file, with duplicates removed.  The "unique.pl" Perl script is as follows, taken from the `Tricks for wgrib2 <https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/tricks.wgrib2>`_ website:
                      
                      .. code-block:: console
                            
                         ----------------------- unique.pl ------------------------
                         #!/usr/bin/perl -w
                         # print only lines where fields 3..N are different
                         #
                         while (<STDIN>) {
                            chomp;
                            $line = $_;
                            $_ =~ s/^[0-9.]*:[0-9]*://;
                            if (! defined $inv{$_}) {
                              $inv{$_} = 1;
                              print "$line\n";
                            }
                         }
                         --------------------- end unique.pl ----------------------

Regional chgres_cube namelist options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Namelist variables with “input” in their name refer to data input to chgres_cube.  Namelist variables with “target” in their name refer to the FV3-LAM horizontal and vertical grid (i.e., the target grid chgres_cube is mapping to).

**Required Entries**

      * fix_dir_target_grid - Path to the FV3-LAM surface climatological files (such as albedo).
      * mosaic_file_target_grid - Path and name of the FV3-LAM mosaic file.
      * orog_dir_target_grid - Directory containing the FV3-LAM orography and grid files (NetCDF).
      * orog_files_target_grid - Names of the FV3-LAM orography file.
      * vcoord_file_target_grid - Path and name of the model vertical coordinate definition file (“global_hyblev.l$LEVS.txt).
      * data_dir_input_grid - Directory containing the GRIB2 initial conditions data
      * grib2_file_input_grid - Name of the GRIB2 input data file
      * varmap_file - Path and name of the variable mapping (VARMAP) table.  See below for details on this table.
      * input_type - Input data type. Set to ‘grib2’
      * cycle_mon/day/hour - Month/day/hour of your model initialization
      * convert_atm - Set to ‘true’ to process atmospheric fields
      * convert_sfc - Set to ‘true’ to process surface fields
      * regional
 
              * Set to 0 to create initial condition atmospheric file
              * Set to 1 to create initial condition atmospheric file and zero hour boundary condition file
              * Set to 2 to create a boundary condition file. Use this option for all but the initialization time.
      * halo_blend - Integer number of row/columns to apply halo blending into the domain, where model and lateral boundary tendencies are applied.
      * halo_bndy - Integer number of rows/columns that exist within the halo, where pure lateral boundary conditions are applied.
      * external_model - Name of source model for input data. Valid options: 'GFS', 'NAM', 'RAP', 'HRRR', 'RRFS'. (Default: 'GFS')

**Optional Entries**

      * geogrid_file_input_grid - Full path to the RAP or HRRR geogrid file corresponding to the external model input data. Only used with external_model = ‘HRRR’ or ‘RAP’. 
      * nsoill_out - Number of soil levels to produce in the sfc_data.nc file (Default: 4).
      * sotyp_from_climo - Use soil type from climatology. Valid options: .true. or .false. (Default: .true.)
      * vgtyp_from_climo - Use vegetation type from climatology. Valid Options: .true. or  .false. (Default: .true.)
      * vgfrc_from_climo - Use vegetation fraction from climatology. Valid options: .true. or .false. (Default: .true.)
      * lai_from_climo - Use leaf area index from climatology. Valid options: .true. or .false. (Default: .true.)
      * minmax_vgfrc_from_climo - Use min/max vegetation fraction from climatology. Valid options: .true. or .false. (Default: .true.)
      * tg3_from_soil - Use tg3 from input soil. Valid options: .true. or .false. . Default: .false.
      * thomp_mp_climo_file - Location of Thompson aerosol climatology file. Provide only if you wish to use these aerosol variables.
      * wam_cold_start - Cold start for the Whole Atmosphere Model. Valid Options: .true. or .false. (Default: .false.)
      * use_rh - Use relative humidity instead of specific humidity when reading in external model grib2 files (Default: .false.)
      * calrh - Type of relative humidity to specific humidity calculation to use (Default: 0; use existing calculation, or 1; use calculation consistent with GFSv15/v16)

Variable Mapping (VARMAP) table
-------------------------------

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
      * "intrp": Ln(pressure) interpolation to missing levels. Linear interpolation and extrapolation are possible, but require modifying the value of "LINLOG" in input_data.F90 to anything other than 2, or to a negative number, respectively.
      * "set_to_fill": Set to user-specified field value (see column 4).
      * "stop": Force an exception and stop code execution. Use this if you absolutely require a field to be present.

Column 4: If column 3 = "set_to_fill", then this value is used to fill in all points in the input field. These values may be overwritten by the code before output depending on the variable (especially for surface variables).

Column 5: Variable type descriptor. Variable names designated as tracers are used to populate the list of tracers to read from the GRIB2 file and write to output, so make sure all tracers you wish to read have an entry. Note that if you wish to add a tracer name that is not already included in the appropriate VARMAP file, this will require modification of the chgres_cube code. Valid choices are:

      * “T”: 3-dimensional tracer array
      * “D”: 3-dimensional non-tracer array
      * “S”: 2-dimensional surface array

Running the program stand alone
-------------------------------

      * Locate your input files.  See above for a list.
      * Set the namelist for your experiment.  See above for an explanation of the namelist entries.
      * Link the namelist to Fortran unit number 41, i.e.”
        * ln -fs your-namelist-file  ./fort.41
      * Load any required runtime libraries.  For example, you may need to load libraries for NetCDF and/or your Fortran compiler.
      * Run the program with an MPI task count that is a multiple of six.  This is an ESMF library requirement when processing a six-tiled global grid.

Making changes to the chgres_cube program
-----------------------------------------

chgres_cube is part of the UFS_UTILS repository (https://github.com/ufs-community/UFS_UTILS). When wanting to contribute to this repository developers shall follow the Gitflow software development process

      * Developers shall create their own fork of the UFS_UTILS repository
      * Developers shall create a ‘feature’ branch off ‘develop’ in their fork for all changes.
      * Developers shall open an issue and reference it in all commits.

For more details, see the UFS_UTILS wiki page: https://github.com/ufs-community/UFS_UTILS/wiki

Changes that support current or future NCEP operations will be given priority for inclusion into the authoritative repository.
