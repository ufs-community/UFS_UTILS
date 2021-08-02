# chgres_cube

# Introduction

The program chgres.F90 creates initial condition files to “coldstart”
the forecast model. The initial conditions are created from either
Global Forecast System (GFS) gridded binary version 2 (GRIB2), NOAA
Environmental Modeling System Input/Output (NEMSIO) data, or Network
Common Data Form (NetCDF) data.

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The chgres_cube program is part of the [NCEPLIBS
UFS_UTILS](https://github.com/NOAA-EMC/UFS_UTILS) project.

## Where to find GFS GRIB2, NEMSIO and NetCDF data

### GRIB2

- 0.25-degree data (last 10 days only) - Use the
  <b>gfs.tHHz.pgrb2.0p25.f000</b> files in subdirectory
  <b>gfs.YYYYMMDD/HH</b> here:
  <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>.

- 0.5-degree data - Use the <b>gfs_4_YYYYMMDD_00HH_000.grb2</b> file,
  under <b>GFS Forecasts 004 (0.5-deg)</b> from the NCDC - Global
  Forecast System
  <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs>. Note:
  Tests were not done with the AVN, MRF or analysis data.

- 1.0-degree data - Use the <b>gfs_3_YYYYMMDD_00HH_000.grb2 file</b>,
  under <b>GFS Forecasts 003 (1-deg)</b> from NCDC - Global Forecast
  System
  <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs>. Note:
  Tests were not done with the AVN, MRF or analysis data.

### NEMSIO

- T1534 gaussian (last 10 days only) - Use the
  <b>gfs.tHHz.atmanl.nemsio</b> (atmospheric fields) and
  <b>gfs.tHHz.sfcanl.nemsio</b> (surface fields) files in subdirectory
  gfs.YYYYMMDD/HH from
  <https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod>.

### NetCDF

- T1534 gaussian (don't have any more details at this time).

# Initializing with GRIB2 data - some caveats

Keep this in mind when using GFS GRIB2 data for model initialization:

- GRIB2 data does not contain the fields needed for the Near Sea
  Surface Temperature (NSST) scheme. See the next section for options
  on running the forecast model in this situation.

- Data is coarse (in vertical and horizontal) compared to the NCEP
  operational GFS . May not provide a good initialization (especially
  for the surface). Recommendations:

 - C96 - use 0.25, 0.5 or 1.0-degree GRIB2 data
 - C192 - use 0.25 or 0.5-degree GRIB2 data
 - C384 - use 0.25-degree GRIB2 data
 - C768 - try the 0.25-degree GRIB2 data. But it may not work well.

- Sea/lake ice thickness and column temperatures are not
  available. So, nominal values of 1.5 m and 265 K are used.

- Soil moisture in the GRIB2 files is created using bilinear
  interpolation and, therefore, may be a mixture of values from
  different soil types. Could result in poor latent/sensible heat
  fluxes.

- Ozone is not available at all isobaric levels. Missing levels are
  set to a nominal value defined in the variable mapping (VARMAP) file
  (1E-07).

- Only tested with GRIB2 data from GFS v14 and v15 (from 12z July 19,
  2017 to current). May not work with older GFS data. Will not work
  with GRIB2 data from other models.

### Near Sea Surface Temperature (NSST) data and GRIB2 initialization

The issue with not having NSST data is important. In GFS we use the
foundation temperature (Tref) and add a diurnal warming/cooling layer
using NSST. This is the surface temperature that is passed to the
atmospheric boundary layer. This is a critical feature, especially
when we are doing Data Assimilation.

When using NEMSIO or NetCDF data to initialize the model, both the
foundation and surface temperature are available and the atmospheric
model should be run using the NSST option as this will properly
account for in the forward run of the model.

In GRIB2 files only the Tsfc is stored and that is set as foundation
temperature as well. So the diurnal heating / cooling is baked into
the initial condition for the extent of the run. This can be critical
if the model is being initialized when the ocean is warm and
initialization is occuring at the peak of the diurnal warming. That
warm ocean will be baked in for the extent of the run and may spawn
off a number of fake hurricanes. The user has two options -- either to
use a spin up cycle to spin up NSST (set <b>nstf_name</b> =
[2,1,0,0,0] in <b>input.nml</b> of the model namelist file. This will
create a diurnal cycle after 24 hours of spin up), or to run the model
without any NSST option ( set <b>nstf_name</b> = [0,0,0,0,0] in
<b>input.nml</b> of the model namelist file. The user will also have
to choose one of the no NSST physics suite options in
<b>input.nml</b>).

Note, that neither of these two options will get rid of the underlying
baked in heating/cooling in the surface temperature fields. For most
cases this may not be an issue, but where it is then the user will
either have to initialize the model with NEMSIO or NetCDF data, or
replace the surface temperature in the GRIB2 fields with independently
obtained foundation temperature.

# chgres_cube namelist options

Namelist variables with “input” in their name refer to data input to
chgres_cube. Namelist variables with “target” in their name refer to
the FV3 horizontal and vertical grid (i.e., the target grid
chgres_cube is mapping to).

When using GRIB2 data as input to chgres_cube, set namelist as
follows:

 - fix_dir_target_grid - Path to the tiled FV3 surface climatological
   files (such as albedo).
 
 - mosaic_file_target_grid - Path and name of the FV3 mosaic file.
 
 - orog_dir_target_grid - directory containing the tiled FV3 orography
   and grid files (NetCDF).
 
 - orog_files_target_grid - names of the six tiled FV3 orography
   files.
 
 - vcoord_file_target_grid - path and name of the model vertical
   coordinate definition file (“global_hyblev.l$LEVS.txt).
 
 - data_dir_input_grid - directory containing the GRIB2 initial
   conditions data
 
 - grib2_file_input_grid - name of the GRIB2 input data file
 
 - varmap_file - path and name of the variable mapping (VARMAP) table.
   See below for details on this table.
 
 - input_type - input data type. Set to ‘grib2’
 
 - cycle_mon/day/hour - month/day/hour of your model initialization
 
 - convert_atm - set to ‘true’ to process the atmospheric fields
 
 - convert_sfc - set to ‘true’ to process the surface fields

When using NEMSIO data as input to chgres_cube, set namelist as follows:

 - fix_dir_target_grid - Path to the tiled FV3 surface climatological
   files (such as albedo).
 
 - mosaic_file_target_grid - Path and name of the FV3 mosaic file.
 
 - orog_dir_target_grid - directory containing the tiled FV3 orography
   and grid files (NetCDF).
 
 - orog_files_target_grid - names of the six tiled FV3 orography
   files.
 
 - vcoord_file_target_grid - path and name of the model vertical
   coordinate definition file (“global_hyblev.l$LEVS.txt).
 
 - data_dir_input_grid - directory containing the NEMSIO input data
 
 - atm_files_input_grid - name of the NEMSIO input atmospheric data
   file
 
 - sfc_files_input_grid - name of the NEMSIO input surface/Near Sea
   Surface Temperature (NSST) data file
 
 - input_type - input data type. Set to ‘gaussian_nemsio’.
 
 - cycle_mon/day/hour - month/day/hour of your model run
 
 - convert_atm - set to ‘true’ to process the atmospheric fields
 
 - convert_sfc - set to ‘true’ to process the surface fields
 
 - convert_nst - set to ‘true’ to process NSST fields
 
 - tracers_input - names of tracer records in input file. For GFDL
   microphysics, set to
   “spfh”,”clwmr”,”o3mr”,”icmr”,”rwmr”,”snmr”,”grle”.
 
 - tracers - names of tracer records in output file expected by model.
   For GFDL microphysics, set to
   “sphum”,”liq_wat”,”o3mr”,”ice_wat”,”rainwat”,”snowwat”,”graupel”.
 
When using NetCDF data as input to chgres_cube, set namelist as follows:

 - fix_dir_target_grid - Path to the tiled FV3 surface climatological
   files (such as albedo).
 
 - mosaic_file_target_grid - Path and name of the FV3 mosaic file.
 
 - orog_dir_target_grid - directory containing the tiled FV3 orography
   and grid files (NetCDF).
 
 - orog_files_target_grid - names of the six tiled FV3 orography
   files.
 
 - vcoord_file_target_grid - path and name of the model vertical
   coordinate definition file (“global_hyblev.l$LEVS.txt).
 
 - data_dir_input_grid - directory containing the NetCDF input data
 
 - atm_files_input_grid - name of the NetCDF input atmospheric data
   file
 
 - sfc_files_input_grid - name of the NetCDF input surface/Near Sea
   Surface Temperature (NSST) data file
 
 - input_type - input data type. Set to ‘gaussian_netcdf’.
 
 - cycle_mon/day/hour - month/day/hour of your model run
 
 - convert_atm - set to ‘true’ to process the atmospheric fields
 
 - convert_sfc - set to ‘true’ to process the surface fields
 
 - convert_nst - set to ‘true’ to process NSST fields
 
 - tracers_input - names of tracer records in input file. For GFDL
   microphysics, set to
   “spfh”,”clwmr”,”o3mr”,”icmr”,”rwmr”,”snmr”,”grle”.
 
 - tracers - names of tracer records in output file expected by model.
   For GFDL microphysics, set to
   “sphum”,”liq_wat”,”o3mr”,”ice_wat”,”rainwat”,”snowwat”,”graupel”.

# Program inputs and outputs

## Inputs

The following four sets of files are located here:
https://ftp.emc.ncep.noaa.gov/EIB/UFS/global/fix/fix_fv3_gmted2010.v20191213/

 - FV3 mosaic file - (NetCDF format)
   - CRES_mosaic.nc

 - FV3 grid files - (NetCDF format)
   - CRES_grid.tile1.nc
   - CRES_grid.tile2.nc
   - CRES_grid.tile3.nc
   - CRES_grid.tile4.nc
   - CRES_grid.tile5.nc
   - CRES_grid.tile6.nc

 - FV3 orography files - (NetCDF format)
   - CRES_oro_data.tile1.nc
   - CRES_oro_data.tile2.nc
   - CRES_oro_data.tile3.nc
   - CRES_oro_data.tile4.nc
   - CRES_oro_data.tile5.nc
   - CRES_oro_data.tile6.nc

 - FV3 surface climatological files - Located under the ./fix_sfc sub-directory. One file for each tile. NetCDF format.
   - CRES.facsf.tileX.nc (fractional coverage for strong/weak zenith angle dependent albedo)
   - CRES.maximum_snow_albedo.tileX.nc (maximum snow albedo)
   - CRES.slope_type.tileX.nc (slope type)
   - CRES.snowfree_albedo.tileX.nc (snow-free albedo)
   - CRES.soil_type.tileX.nc (soil type)
   - CRES.subtrate_temperature.tileX.nc (soil substrate temperature)
   - CRES.vegetation_greenness.tileX.nc (vegetation greenness)
   - CRES.vegetation_type.tileX.nc (vegetation type)

 - FV3 vertical coordinate file. Text file. Located here: https://ftp.emc.ncep.noaa.gov/EIB/UFS/global/fix/fix_am.v20191213/
   - global_hyblev.l$LEVS.txt

 - Input data files. GRIB2, NEMSIO or NetCDF. See above section for how to find this data.

## Outputs

 - Atmospheric “coldstart” files. NetCDF.
   - out.atm.tile1.nc
   - out.atm.tile2.nc
   - out.atm.tile3.nc
   - out.atm.tile4.nc
   - out.atm.tile5.nc
   - out.atm.tile6.nc

 - Surface/Near Sea Surface Temperature (NSST) “coldstart” files. NetCDF
   - out.sfc.tile1.nc
   - out.sfc.tile1.nc
   - out.sfc.tile1.nc
   - out.sfc.tile1.nc
   - out.sfc.tile1.nc
   - out.sfc.tile1.nc

# Running the program stand alone

 - Locate your input files. See above for a list.
 
 - Set the namelist for your experiment. See above for an explanation
   of the namelist entries.
 
 - Link the namelist to Fortran unit number 41:
 <pre> ln -fs your-namelist-file ./fort.41</pre>
 
 - Load any required runtime libraries. For example, you may need to
   load libraries for NetCDF and/or your Fortran compiler.
 
 - Run the program with an MPI task count that is a multiple of six.
   This is an ESMF library requirement when processing a six-tiled
   global grid.

# Variable Mapping (VARMAP) table

The VARMAP table, set in the chgres_cube namelist (variable
varmap_file), controls how chgres_cube handles variables that might be
missing from the GRIB2 files. Since there are so many different
versions of GRIB2 files, it's often uncertain what fields are
available even if you know what source model the data is coming from.
Each file contains the following: (Note, only the GFS physics suite is
currently supported.)

Column 1: Name the code searches for in the table. Do not change.
Some definitions:

 - dzdt - vertical velocity
 - sphum - specific humidity
 - liq_wat - liquid water mixing ratio
 - o3mr - ozone mixing ratio
 - ice_wat - ice water mixing ratio
 - rainwat - rain water mixing ratio
 - snowwat - snow water mixing ratio
 - graupel - graupel mixing ratio
 - vtype - vegetation type
 - sotype - soil type
 - vfrac - plant greenness fraction
 - fricv - friction velocity
 - sfcr - roughness length
 - tprcp - precipitation rate
 - ffmm - surface exchange coefficient for momentum
 - ffhh - surface exchange coefficient for heat
 - f10m - log((sfcr+10)/sfcr)
 - soilw - total volumetric soil moisture
 - soill - liquid volumetric soil moisture
 - soilt - soil column temperature
 - cnwat - plant canopy water content
 - hice - sea/lake ice thickness
 - weasd - snow liquid equivalent
 - snod - physical snow depth

Column 2: Name of the variable in the output “coldstart”
files. Unimplemented.

Column 3: Behavior when the code can't find the variable in the input
file. Options are:

 - "skip": Don't write to the output file.
 - "set_to_fill": Set to user-specified field value (see column 4).
 - "intrp": LnP interpolation to missing levels. No extrapolation allowd.
 - "stop": Force an exception and stop code execution. Use this if you
   absolutely require a field to be present.

Column 4: If column 3 = "set_to_fill", then this value is used to fill
in all points in the input field. These values may be overwritten by
the code before output depending on the variable (especially for
surface variables).

Column 5: Variable type descriptor. Variable names designated as
tracers are used to populate the list of tracers to read from the
GRIB2 file and write to output, so make sure all tracers you wish to
read have an entry. Note that if you wish to add a tracer name that is
not already included in the appropriate VARMAP file, this will require
modification of the chgres_cube code. Valid choices are:

 - “T”: 3-dimensional tracer array
 - “D”: 3-dimensional non-tracer array
 - “S”: 2-dimensional surface array


