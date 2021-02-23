.. _ufs_utils:


****************************
Introduction
****************************

The Unified Forecast Systems (UFS) Utilities repository contains pre-processing programs for the UFS weather model.  These programs set up the model grid and create coldstart initial conditions. The repository is hosted on `Github <https://github.com/NOAA-EMC/UFS_UTILS>`_.  Information on checking out the code and making changes to it is available on the repository `wiki page <https://github.com/NOAA-EMC/UFS_UTILS/wiki>`_.

***********************************
Grid Generation
***********************************

The following programs are used to create a grid.

      * make_hgrid
      * regional_grid_esg
      * make_solo_mosaic
      * orog
      * global_equiv_resol
      * shave
      * filter_topo
      * sfc_climo_gen

The grid generation process is run by these scripts (located under ./ush)

      * fv3gfs_grid_driver.sh  (driver script)
      * fv3gfs_make_grid.sh (creates the geo-referencing for the grid)
      * fv3gfs_make_orog.sh (creates the land-sea mask and terrain)
      * fv3gfs_filter_topo.sh (filters the orography) 
      * sfc_climo_gen.sh (creates climatological surface fields, such as soil type)

***************************************************
Description of each program
***************************************************

chgres_cube
===========

.. include:: chgres_cube.rst

make_hgrid
==========

Introduction
------------

The make_hgrid program computes geo-referencing parameters for global uniform grids.  (Extended Schmidt gnomonic regional grids are created by the regional_esg_grid program.)  The parameters include geographic latitude and longitude, and grid cell area.  See the output data section for a full list of parameters.  Grids are gnomonic such that all great circles are straight lines.  The parameters are computed on the staggered or "supergrid" - which has twice the resolution of the model grid.  The chgres_cube initialization program maps mass fields - such as temperature - at the supergrid centroids, and u/v winds at the face mid-points.  The supergrid is shown here:

.. _figure_reference:

.. figure:: _static/supergrid.png

Code Structure
--------------

Location of source code: ./sorc/fre-nctools.fd/tools/make_hgrid.  Relevant routines:

      * make_hgrid.c - main driver routine
      * create_gnomonic_cubic_grid.c - contains routines for creating a gnomonic grid

Namelist Options
----------------

The program is controlled by these script variables:

      * Global uniform grid
             * res - x/y dimension of one tile.  The "CRES" resolution.  Example: a 96x96 tile would be classified as C96.  It may be converted to physical resolution as follows: resol = (360 degrees / 4*CRES) * 111 km. 

Program inputs and outputs
--------------------------

**Input data:**

None

**Output data:**

Tiled "grid" files (NetCDF) containing geo-referencing records.  File naming convention: CRES_grid.tile#.nc. File records include:

      * x - geographic longitude (degrees)
      * y - geographic latitude (degrees)
      * dx - grid edge 'x' distance (m)
      * dy - grid edge 'y' distance (m)
      * area - grid cell area (m^2)
      * angle_dx - grid vertex 'x' angle with respect to geographic east (degrees)
      * angle_dy - grid vertex 'y' angle with respect to geographic north (degrees)

regional_esg_grid
=================

Introduction
------------

The regional_esg_grid program computes geo-referencing parameters for the Extended Schmidt Gnomonic (ESG) regional grid.  The parameters include geographic latitude and longitude, and grid cell area.  See the output data section for a full list of parameters.  The ESG grid is designed to have nearly homogenous grid spacing.  Like the make_hgrid program, the parameters are computed on the staggered or "supergrid".  For more information on the Extended Schmidt Gnomonic, see: `Purser, et. al <https://dtcenter.org/sites/default/files/events/2020/2-purser-james.pdf>`_.

Code Structure
--------------

Location of source code: ./sorc/grid_tools.fd/regional_esg_grid.fd.  Relevant routines:

      * regional_esg_grid.f90 - Main driver routine.  Reads program namelist.  Writes output file.
      * pseg.f90 - Suite of routines to perform the ESG regional grid mapping.

Namelist Options
----------------

The program is controlled by these script variables:

      * target_lon - center longitude of grid - degrees
      * target_lat - center latitude of grid - degrees
      * idim - dimension of grid in 'i' direction
      * jdim - dimension of grid in 'j' direction
      * delx - grid spacing in degrees in the 'i' direction on the supergrid.  The physical grid spacing is related to delx as follows: 2*delx(circumf_earth / 360 deg).
      * dely - grid spacing in degrees in the 'j' direction on the supergrid.
      * halo - number of rows/cols of the lateral boundary halo

Program Inputs and Outputs
--------------------------

**Input data:**

None

**Output data:**

A tiled "grid" file (NetCDF) containing geo-referencing records for the supergrid.  File naming convention: CRES_grid.tile7.nc.  Here, CRES is the global equivalent resolution as computed by the global_equiv_resol program.  Note: the forecast model assumes regional grids are tile 1.  File records include:

      * x - geographic longitude (degrees)
      * y - geographic latitude (degrees)
      * dx - grid edge 'x' distance (m)
      * dy - grid edge 'y' distance (m)
      * area - grid cell area (m^2)
      * angle_dx - grid vertex 'x' angle with respect to geographic east (degrees)
      * angle_dy - grid vertex 'y' angle with respect to geographic north (degrees)


make_solo_mosaic
================

Introduction
------------

This program creates the "mosaic" file, which contains information about the tiled "grid" files (such as name and tile number).  For global grids, it also defines the orientation of the six tiles.  All output records are listed below under "output data".  There are no runtime-selectable options for this program.

Code structure
--------------

Location of source code ./sorc/fre-nctools.fd/tools/make_solo_mosaic.  Relevant routines:

      * make_solo_mosaic.c - main driver routine.
      * get_contact.c - computes the number of aligned contacts between two tiles.


Program inputs and outputs
--------------------------

**Input data:**  

The tiled "grid" files (CRES_grid.tile#.nc) created by the make_hgrid or regional_esg_grid programs - (NetCDF)

**Output data:** 

The mosaic file - CRES_mosaic.nc (NetCDF).  Contains these records

      * Mosaic - name of mosaic (character)
      * Gridlocation - directory containing the "grid" files (character)
      * Gridfiles - names of each "grid" file (character array)
      * Gridtiles - list of each tile number (character array)
      * Contacts - list of tile contact regions - global grids only (character array)
      * Contact_index - list of contact regions as specified by i/j index - global grids only (character array).

global_equiv_resol
==================

Introduction
------------

This program computes the global equivalent resolution for regional grids.  For example, a global grid with x/y dimensions of 96 would have a global resolution (CRES) of C96.  And the approximate physical resolution would be:

      * Res in km = (360 degrees / 4*CRES) * 111 km = 104 km 

Using the average cell size (in m^2) of the regional grid, the equivalent global resolution is computed according to:

      * CRES = nint( (2*pi*rad_earth)/(4*avg_cell_size) )

There are no runtime-selectable options for this program.

Code structure
--------------

Location of source code:  ./sorc/grid_tools.fd/global_equiv_resol.fd.  Relevant routine:

      * global_equiv_resol.f90 - Contains the entire program.

Program inputs and outputs
--------------------------

**Input data:**  

The regional "grid" file (CRES_grid.tile#.nc) created by the regional_esg_grid program - (NetCDF).  Uses the grid cell area record.

**Output data:**  

Adds the equivalent resolution as an attribute to the "grid" file - CRES_grid.tile#.nc  (NetCDF).

orog
====

Introduction
------------

This program computes the land mask, land fraction, orography and gravity wave drag (GWD) fields on the model grid.  See the output data section for a complete list of fields.

Land-sea mask and land fraction are created from a global 30-arc second University of Maryland land cover (land/non-land flag) dataset.  Land fraction is determined by averaging all 30-second land cover points located within the model grid box.   Points with a land fraction of 50% or more are given a land-sea mask of "land".  Orography and GWD fields are created from two datasets: 1) 30-arc-second `USGS GMTED2010 <https://www.usgs.gov/core-science-systems/eros/coastal-changes-and-impacts/gmted2010?qt-science_support_page_related_con=0#qt-science_support_page_related_con>`_ orography data; 2) for Antarctica, 30-arc-second `RAMP <https://nsidc.org/data/nsidc-0082>`_ terrain data (Radarsat Antarctic Mapping Project).  Fields are determined from all 30-arc second data located within the model grid box.  The orography is simply the average of the 30-arc second values.  It is later filtered by the filter_topo program.  For details on the GWD fields, see:

      * Kim, Y-J and A. Arakawa, 1995: Improvement of orographic gravity wave parameterization using a mesoscale gravity wave model.  J. Atmos. Sci. 52, pp 1875-1902.
      * Lott, F. and M. J. Miller: 1977: A new sub-grid scale orographic drag parameterization: Its formulation and testing, QJRMS, 123, pp 101-127.

Code structure
--------------

The source code is located - ./sorc/orog_mask_tools.fd/orog.fd.  Some important subroutines:

      * MAKEMT2 - computes land fraction, land-sea mask, orography, standard deviation of orography, and convexity.  
      * MAKEPC2 - computes anisotropy (gamma), slope of orography (sigma) and mountain range angle (theta).
      * MAKEOA2 - computes maximum height (elvmax), orographic asymmetry (oa) and length scale (ol).

Program inputs and outputs
--------------------------

**Input data:**

      * The "grid" files (CRES_grid.tile#.nc) containing the geo-reference records for the grid - (NetCDF).  Created by the make_hgrid or regional_esg_grid programs.
      * Global 30-arc-second University of Maryland land cover data.  Used to create the land-sea mask.
             * ./fix/fix_orog/landcover30.fixed (unformatted binary)
      * Global 30-arc-second USGS GMTED2010 orography data.
             * ./fix/fix_orog/gmted2010.30sec.int (unformatted binary)
      * 30-arc-second RAMP Antarctic terrain data (Radarsat Antarctic Mapping Project)
             * ./fix/fix_orog/thirty.second.antarctic.new.bin (unformatted binary)

**Output data:**  

Orography files - one for each tile - oro.CRES.tile#.nc (NetCDF).  Contains these records:

      * geolon - longitude (degrees east)
      * geolat - latitude (degrees north)
      * slmsk - land-sea mask (0 - nonland; 1 - land)
      * land_frac - land fraction (percent)
      * orog_raw - orography (meters)
      * orog_filt - same as orog_raw
      * stddev - standard deviation of orography (meters) 
      * convexity - orographic convexity (unitless)
      * oa[1-4] - orographic asymmetry - four directional components - W/S/SW/NW
      * ol[1-4] - orographic length scale - four directional components - W/S/SW/NW
      * theta - angle of mountain range with respect to east (degrees)
      * gamma - anisotropy (unitless)
      * sigma - slope of orography (unitless)
      * elvmax - maximum height above mean (meters)

filter_topo
===========

Introduction
------------

The FV3 terrain filtering algorithm has several unique properties compared to conventional topography filters. The resulting topography filtered by this algorithm has conserved globally integrated elevations. More importantly, this filter has the following island-preserving properties: 1) No Gibbs ringing at the coastlines where discontinuities occur; 2) the filtered terrain's coastlines strictly match the source terrain's. The detailed implementation of this terrain filtering algorithm will be described in a forthcoming publication by Dr. Shian-Jiann Lin and his group.

Code structure
--------------

Location of source code: ./sorc/grid_tools.fd/filter_topo.fd. The entire program is contained in filter_topo.F90.

Namelist options
----------------

Program execution is controlled via a namelist.  The namelist variables are:

      * topo_file - Name of the orography file (See input data) (character)
      * topo_field - Name of the filtered orography record in the orography file ("orog_filt") (character)
      * mask_field - Name of the land-sea mask record in the orography file ("land_frac") (character)
      * grid_file - The mosaic file (See input data) (character)
      * zero_ocean - Flag to turn on the "island-preserving" property.  Default is true (logical)
      * stretch_fac - Stretching factor.  Equal to "1" for global uniform grids. Not applicable for ESG regional grids (floating point)
      * res - The "CRES" resolution (floating point)
      * grid_type - 0 for a gnomonic grid (integer)
      * regional - True for an ESG regional grid (logical)

Program inputs and outputs
--------------------------

**Input data:**

      * mosaic file - the mosaic file from the make_solo_mosaic program - CRES_mosaic.nc (NetCDF)
      * grid file - the "grid" file from the make_hgrid or regional_esg programs  - CRES_grid.tile#.nc - (NetCDF)
      * orography file - the orography file from the orog program - oro.CRES.tile#.nc (NetCDF)

**Output data:**

      * The filtered orography is written to the "orog_filt" record of the input orography file - oro.CRES.tile#.nc (NetCDF).

Filtering parameters
--------------------

      * n_del2 - Second-order strong filtering coefficient.
      * n_del2_weak - Second-order weak filtering coefficient - used to more finely smooth the topography compared to the strong filter.
      * cd4 - dimensionless coefficient for delta-4 diffusion.
      * peak_fac - overshoot factor for the mountain peak
      * max_slope - maximum allowable terrain slope

shave
=====

Introduction
------------

The "grid" and "orography" files for regional grids are first created with rows and columns extending beyond the halo.  This is required for the topography filtering code to work correctly in the halo region.  After filtering, the shave program removes these extra points from the files.

Code structure
--------------

Location of the source code: ./sorc/grid_tools.fd/shave.fd.  The entire program is contained in shave_nc.F90.

Program control options
-----------------------

The program is controlled by these parameters read from standard input:

      * The i/j dimensions of the compute domain (not including halo) (integer)
      * The number of halo rows/columns (integer)
      * The file name with the extra points (character string)
      * The name of the output file with the extra points removed (character string)

Program inputs and outputs
--------------------------

**Input data:**

      * Model "grid" files (CRES_grid.tile#.nc) created by the make_hgrid or regional_esg_grid programs - (NetCDF)
      * Model orography files (oro.CRES.tile#.nc)  after topography filtering - (NetCDF)

**Output data:** 

With and without the halo.

      * Model "grid" files - CRES.grid.tile#.halo#.nc (NetCDF) 
      * Model orography files - CRES.oro_data.tile#.halo#.nc (NetCDF)


sfc_climo_gen
=============

Introduction
------------

The sfc_climo_gen (surface climatological field generation) program creates surface climatological fields such as soil type, vegetation type and albedo.  Some fields may be time-varying: for example, snow-free albedo is monthly.  But they are static - i.e., they only need to be generated when the grid is created.  The program uses the ESMF library to horizontally interpolate the source data to the model grid.  For regional grids, the program will output files with and without the halo region when the "halo" namelist variable is set.

Code structure
--------------

Location of the source code: ./sorc/sfc_climo_gen.fd.  Brief description of each module:

      * driver.F90 - The main driver routine.
      * interp.F90 - The interpolation driver routine.  Reads the input source data and interpolates it to the model grid.
      * model_grid.F90 - Defines the ESMF grid object for the model grid.
      * output.f90 - Writes the output surface data to a NetCDF file.  For regional grids, will output separate files with and without the halo.
      * program_setup.f90 - Reads the namelist and sets up program execution.
      * search.f90 - Replace undefined values on the model grid with a valid value at a nearby neighbor. Undefined values are typically associated with isolated islands where there is no source data.
      * source_grid.F90 - Reads the grid specifications and land/sea mask for the source data.  Sets up the ESMF grid object for the source grid.
      * utils.f90 - Contains error handling utility routines.

Namelist options
----------------

Program execution is controlled via a namelist.  The namelist variables are:

      * input_facsf_file - path/name of input fractional strong/weak zenith angle albedo data
      * input_substrate_file - path/name of input soil substrate temperature data
      * input_maximum_snow_albedo_file - path/name of input maximum snow albedo data
      * input_snowfree_albedo_file - path/name of input snow-free albedo data
      * input_slope_type_file - path/name of input global slope type data
      * input_soil_type_file - path/name of input soil type data
      * input_vegetation_type_file - path/name of vegetation type data
      * input_vegetation_greenness_file - path/name of monthly vegetation greenness data
      * mosaic_file_mdl - path/name of the model mosaic file
      * orog_dir_mdl - directory containing the model orography files
      * orog_files_mdl - list of model orography files.  For global uniform grids, all six files are listed.
      * halo - number of rows/cols of the lateral boundary halo (regional grids only).  When selected, the program will output files with and without the halo region.
      * maximum_snow_albedo_method - interpolation method for this field.  Bilinear or conservative.  Default is bilinear.
      * snowfree_albedo_method -  interpolation method for this field.  Bilinear or conservative.  Default is bilinear.
      * vegetation_greenness_method -  interpolation method for this field.  Bilinear or conservative.  Default is bilinear.

Program inputs and outputs
--------------------------

**Input data:** 

The global surface climatological data is located in ./fix/fix_sfc_climo.  All NetCDF.

      * Global 1-degree fractional coverage strong/weak zenith angle albedo - facsf.1.0.nc
      * Global 0.05-degree maximum snow albedo - maximum_snow_albedo.0.05.nc
      * Global 2.6 x 1.5-degree soil substrate temperature - substrate_temperature.2.6x1.5.nc
      * Global 0.05-degree four component monthly snow-free albedo - snowfree_albedo.4comp.0.05.nc
      * Global 1.0-degree categorical slope type - slope_type.1.0.nc
      * Global 0.05-degree categorical STATSGO soil type - soil_type.statsgo.0.05.nc
      * Global 0.05-degree categorical IGBP vegetation type - vegetation_type.igbp.0.05.nc
      * Global 0.144-degree monthly vegetation greenness in percent - vegetation_greenness.0.144.nc
      * Model mosaic file - CRES_mosaic.nc (NetCDF)
      * Model orography files including halo - CRES_oro_data.tile#.halo#.nc (NetCDF)
      * Model grid files including halo - CRES_grid.tile#.halo#.nc (NetCDF)

**Output files:** 

All files with and without halo (all NetCDF).

      * Fractional coverage strong/weak zenith angle albedo - CRES_facsf.tile#.halo#.nc
      * Maximum snow albedo - CRES_maximum_snow_albedo.tile#.halo#nc
      * Soil substrate temperature - CRES_substrate_temperature.tile#.halo#.nc
      * Snow free albedo - CRES_snowfree_albedo.tile#.halo#.nc
      * Slope type - CRES_slope_type.tile#.halo#.nc
      * Soil type - CRES_soil_type.tile#.halo#.nc
      * Vegetation type - CRES_vegetation_type.tile#.halo#.nc
      * Vegetation greenness - CRES_vegetation_greenness.tile#.halo#.nc
